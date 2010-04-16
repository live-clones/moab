#include "iMeshP.h"
#include "iMesh_MOAB.hpp"
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/CN.hpp"
#include "moab/MeshTopoUtil.hpp"
#include "FileOptions.hpp"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"

#define IS_BUILDING_MB
#include "Internals.hpp"
#undef IS_BUILDING_MB

#include <assert.h>
#include <sstream>

#ifdef USE_MPI    
#include "moab_mpi.h"
#endif

using namespace moab;


/********************* Error Handling **************************/

#define FIXME printf("Warning: function has incomplete implementation: %s\n", __func__ )



/******** Type-safe casting between MOAB and ITAPS types *********/

#ifndef TEMPLATE_FUNC_SPECIALIZATION
// if no template specializtion, disable some type checking
template <typename T, typename S> inline
T itaps_cast( S handle )
{ 
  assert(sizeof(S) >= sizeof(T)); 
  return reinterpret_cast<T>(handle);
}
#else

// basic template method : only works to cast to equivalent types (no-op)
template <typename T, typename S> inline
T itaps_cast( S h )
{ return h; }
// verify size and do reinterpret cast
template <typename T> inline T itaps_cast_internal_( EntityHandle h )
{
  assert(sizeof(T) >= sizeof(EntityHandle));
  return reinterpret_cast<T>(h);
}
// verify size and do reinterpret cast
template <typename T> inline EntityHandle* itaps_cast_ptr_( T* h )
{
  assert(sizeof(T) == sizeof(EntityHandle));
  return reinterpret_cast<EntityHandle*>(h);
}
// verify size and do reinterpret cast
template <typename T> inline const EntityHandle* itaps_cast_const_ptr_( const T* h )
{
  assert(sizeof(T) == sizeof(EntityHandle));
  return reinterpret_cast<const EntityHandle*>(h);
}
// verify set-type handle before cast
template <typename T> inline T itaps_set_cast_( EntityHandle h )
{
  assert(TYPE_FROM_HANDLE(h) == MBENTITYSET);
  return itaps_cast_internal_<T>(h);
}

// define conversion routines between itaps handle and EntityHandle types
#define DECLARE_ALLOWED_ITAPS_CONVERSION( ITAPS_HANDLE_TYPE ) \
  template <> inline \
  ITAPS_HANDLE_TYPE \
  itaps_cast<ITAPS_HANDLE_TYPE,EntityHandle>( EntityHandle h ) \
  { return itaps_cast_internal_<ITAPS_HANDLE_TYPE>(h); } \
  \
  template <> inline \
  EntityHandle \
  itaps_cast<EntityHandle,ITAPS_HANDLE_TYPE>( ITAPS_HANDLE_TYPE handle ) \
  { return reinterpret_cast<EntityHandle>(handle); } \
  \
  template <> inline \
  EntityHandle* \
  itaps_cast<EntityHandle*,ITAPS_HANDLE_TYPE*>( ITAPS_HANDLE_TYPE* ptr ) \
  { return itaps_cast_ptr_(ptr); } \
  \
  template <> inline \
  const EntityHandle* \
  itaps_cast<const EntityHandle*,const ITAPS_HANDLE_TYPE*>( const ITAPS_HANDLE_TYPE* ptr ) \
  { return itaps_cast_const_ptr_(ptr); }


// define conversion routines between itaps handle and EntityHandle types
// but limit to EntityHandle for MBENTITYSET type.
#define DECLARE_ALLOWED_ITAPS_SET_CONVERSION( ITAPS_HANDLE_TYPE ) \
  template <> inline \
  ITAPS_HANDLE_TYPE \
  itaps_cast<ITAPS_HANDLE_TYPE,EntityHandle>( EntityHandle h ) \
  { return itaps_set_cast_<ITAPS_HANDLE_TYPE>(h); } \
  \
  template <> inline \
  EntityHandle \
  itaps_cast<EntityHandle,ITAPS_HANDLE_TYPE>( ITAPS_HANDLE_TYPE handle ) \
  { return reinterpret_cast<EntityHandle>(handle); } \
  \
  template <> inline \
  EntityHandle* \
  itaps_cast<EntityHandle*,ITAPS_HANDLE_TYPE*>( ITAPS_HANDLE_TYPE* ptr ) \
  { return itaps_cast_ptr_(ptr); } \
  \
  template <> inline \
  const EntityHandle* \
  itaps_cast<const EntityHandle*,const ITAPS_HANDLE_TYPE*>( const ITAPS_HANDLE_TYPE* ptr ) \
  { return itaps_cast_const_ptr_(ptr); }

DECLARE_ALLOWED_ITAPS_SET_CONVERSION( iMeshP_PartitionHandle )
//DECLARE_ALLOWED_ITAPS_SET_CONVERSION( iMeshP_PartHandle )
DECLARE_ALLOWED_ITAPS_SET_CONVERSION( iBase_EntitySetHandle )
DECLARE_ALLOWED_ITAPS_CONVERSION( iBase_EntityHandle )


template <> inline
Tag itaps_cast<Tag,iBase_TagHandle>( iBase_TagHandle h )
  { return reinterpret_cast<Tag>(h); }
template <> inline
iBase_TagHandle itaps_cast<iBase_TagHandle,Tag>( Tag h )
  { return reinterpret_cast<iBase_TagHandle>(h); }


#endif

#define PCOMM ParallelComm::get_pcomm( MBI, itaps_cast<EntityHandle>(partition_handle) )

// Need a different function name for Tag because (currently)
// both Tag and iBase_EntityHandle are void**.
iBase_TagHandle itaps_tag_cast( Tag t )
{ 
  assert(sizeof(iBase_TagHandle) >= sizeof(Tag)); 
  return reinterpret_cast<iBase_TagHandle>(t);
}

/********************* ITAPS arrays **************************/

// Handle returning Range in ITAPS array (do ALLOCATE_ARRAY and copy).
#define RANGE_TO_ITAPS_ARRAY( RANGE, NAME ) do { \
  ALLOC_CHECK_ARRAY_NOFAIL( NAME, (RANGE).size() ); \
  std::copy( (RANGE).begin(), (RANGE).end(), itaps_cast<EntityHandle*>(*(NAME)) ); \
  } while (false)


static inline ErrorCode get_entities( Interface* iface,
                                        EntityHandle set,
                                        int type, int topology, 
                                        Range& entities )
{
  if (topology != iMesh_ALL_TOPOLOGIES)
    return iface->get_entities_by_type( set, mb_topology_table[topology], entities );
  else if (type != iBase_ALL_TYPES)
    return iface->get_entities_by_dimension( set, type, entities );
  else
    return iface->get_entities_by_handle( set, entities );
}

static inline ErrorCode remove_not_owned( ParallelComm* pcomm, Range& ents )
{
  ErrorCode rval;
  
  std::vector<unsigned char> pstatus(ents.size());
  rval = pcomm->get_moab()->tag_get_data(pcomm->pstatus_tag(), ents, &pstatus[0]);
  if (MB_SUCCESS != rval)
    return rval;
    
  Range::iterator i = ents.begin();
  std::vector<unsigned char>::const_iterator j;
  for (j = pstatus.begin(); j != pstatus.end(); ++j) {
    if (*j & PSTATUS_NOT_OWNED)
      i = ents.erase( i );
    else
      ++i;
  }
  
  return MB_SUCCESS;
}

static inline ErrorCode count_owned( ParallelComm* pcomm, const Range& ents, int& n )
{
  ErrorCode rval;
  n = 0;
  
  std::vector<unsigned char> pstatus(ents.size());
  rval = pcomm->get_moab()->tag_get_data(pcomm->pstatus_tag(), ents, &pstatus[0]);
  if (MB_SUCCESS != rval)
    return rval;
    
  std::vector<unsigned char>::const_iterator j;
  for (j = pstatus.begin(); j != pstatus.end(); ++j)
    if (!(*j & PSTATUS_NOT_OWNED))
      ++n;
  
  return MB_SUCCESS;
}

static void set_intersection_query( iMesh_Instance instance,
                                    iMeshP_PartHandle set1,
                                    iBase_EntitySetHandle set2,
                                    int type,
                                    int topo,
                                    Range& result,
                                    int* err )
{
  ErrorCode rval;
  
  if (!set1) {
    rval = get_entities( MBI, itaps_cast<EntityHandle>(set2), type, topo, result );
    CHKERR(rval,"Invalid Part handle");
  }
  else if (!set2) {
    rval = get_entities( MBI, itaps_cast<EntityHandle>(set1), type, topo, result );
    CHKERR(rval,"Invalid set handle");
  }
  else {
    Range r1, r2;
    rval = get_entities( MBI, itaps_cast<EntityHandle>(set1), type, topo, r1 );
    CHKERR(rval,"Invalid Part handle");
    rval = get_entities( MBI, itaps_cast<EntityHandle>(set2), type, topo, r2 );
    CHKERR(rval,"Invalid set handle");
    result.merge( intersect( r1, r2) );
  }
  
  RETURN (iBase_SUCCESS);
}
  


/********************* iMeshP API **************************/

#ifdef __cplusplus
extern "C" {
#endif

void iMeshP_createPartitionAll( iMesh_Instance instance,
                        /*in*/  MPI_Comm communicator,
                        /*out*/ iMeshP_PartitionHandle *partition_handle,
                                int *err )
{
  *partition_handle = 0;

  Tag prtn_tag;
  ErrorCode rval = MBI->tag_create( PARALLEL_PARITIONING_TAG_NAME, 
                                      sizeof(int), 
                                      MB_TAG_SPARSE,
                                      MB_TYPE_INTEGER, 
                                      prtn_tag, 
                                      0, 
                                      true ); CHKERR(rval,"tag creation failed");
  
  EntityHandle handle;
  rval = MBI->create_meshset( MESHSET_SET, handle ); CHKERR(rval,"set creation failed");
  ParallelComm* pcomm = ParallelComm::get_pcomm( MBI, handle, &communicator );
  if (!pcomm) {
    MBI->delete_entities( &handle, 1 );
    RETURN(iBase_FAILURE);
  }
  
  *partition_handle = itaps_cast<iMeshP_PartitionHandle>(handle);
  RETURN (iBase_SUCCESS);
}

void iMeshP_destroyPartitionAll( iMesh_Instance instance,
                                 iMeshP_PartitionHandle partition_handle,
                                 int *err)
{
  ParallelComm* pcomm = PCOMM;
  if (pcomm)
    delete pcomm;
  EntityHandle handle = itaps_cast<EntityHandle>(partition_handle);
  ErrorCode rval = MBI->delete_entities( &handle, 1 ); CHKERR(rval,"entity deletion failed");
  RETURN (iBase_SUCCESS);
}

void iMeshP_getPartIdFromPartHandle( iMesh_Instance instance,
                                     const iMeshP_PartitionHandle partition_handle,
                                     const iMeshP_PartHandle part_handle,
                                     iMeshP_Part *part_id,
                                     int *err )
{
  int junk1 = 1, junk2;
  iMeshP_getPartIdsFromPartHandlesArr( instance, partition_handle, &part_handle, 1, 
                                       &part_id, &junk1, &junk2, err );
}

void iMeshP_getPartHandleFromPartId( iMesh_Instance instance,
                                     const iMeshP_PartitionHandle partition_handle,
                                     iMeshP_Part part_id,
                                     iMeshP_PartHandle *part_handle,
                                     int *err )
{
  int junk1 = 1, junk2;
  iMeshP_getPartHandlesFromPartsIdsArr( instance, partition_handle, &part_id, 1, 
                                        &part_handle, &junk1, &junk2, err );
}

void iMeshP_getPartIdsFromPartHandlesArr(
            iMesh_Instance instance,
            const iMeshP_PartitionHandle partition_handle,
            const iMeshP_PartHandle *part_handles,
            const int part_handles_size,
            iMeshP_Part **part_ids,
            int *part_ids_allocated,
            int *part_ids_size,
            int *err)
{
  ErrorCode rval;
  ParallelComm* pcomm = PCOMM;
  ALLOC_CHECK_ARRAY( part_ids, part_handles_size );
  for (int i = 0; i < part_handles_size; ++i) {
    int id;
    rval = pcomm->get_part_id( itaps_cast<EntityHandle>(part_handles[i]), id );
    (*part_ids)[i] = id;
    CHKERR(rval,"error getting part id");
  }
  KEEP_ARRAY(part_ids);
  RETURN(iBase_SUCCESS);
}

void iMeshP_getPartHandlesFromPartsIdsArr(
            iMesh_Instance instance,
            const iMeshP_PartitionHandle partition_handle,
            const iMeshP_Part *part_ids,
            const int part_ids_size,
            iMeshP_PartHandle **part_handles,
            int *part_handles_allocated,
            int *part_handles_size,
            int *err)
{
  ErrorCode rval;
  ParallelComm* pcomm = PCOMM;
  ALLOC_CHECK_ARRAY( part_handles, part_ids_size );
  for (int i = 0; i < part_ids_size; ++i) {
    EntityHandle handle;
    rval = pcomm->get_part_handle( part_ids[i], handle );
    CHKERR(rval,"error getting part handle");
    (*part_handles)[i] = itaps_cast<iMeshP_PartHandle>(handle);
  }
  KEEP_ARRAY(part_handles);
  RETURN(iBase_SUCCESS);
}

void iMeshP_getPartitionComm( iMesh_Instance instance,
                              iMeshP_PartitionHandle partition_handle,
                              MPI_Comm* communicator_out,
                              int* err )
{
  ParallelComm* pcomm = PCOMM;
  if (!pcomm)
    RETURN (iBase_FAILURE);
  *communicator_out = pcomm->proc_config().proc_comm();
  RETURN (iBase_SUCCESS);
}

void iMeshP_syncPartitionAll( iMesh_Instance instance,
                              iMeshP_PartitionHandle partition_handle,
                              int* err )
{
  ParallelComm* pcomm = PCOMM;
  if (!pcomm)
    ERROR (iBase_FAILURE,"No PComm");
  ErrorCode rval = pcomm->collective_sync_partition();
  CHKERR(rval,"collective sync failed");
  RETURN(iBase_SUCCESS);
}

void iMeshP_getNumPartitions( iMesh_Instance instance,
                              int* num_partitions_out,
                              int* err )
{
  std::vector<ParallelComm*> pcomms;
  ErrorCode rval = ParallelComm::get_all_pcomm( MBI, pcomms );
  CHKERR(rval,"Internal error retreiving PComms");
  
  std::vector<ParallelComm*>::iterator i;
  *num_partitions_out = 0;
  for (i = pcomms.begin(); i != pcomms.end(); ++i)
    if ((*i)->get_partitioning())
      (*num_partitions_out)++;

  RETURN (iBase_SUCCESS );
}

void iMeshP_getPartitions( iMesh_Instance instance,
                           iMeshP_PartitionHandle **partition_handle,
                           int *partition_handle_allocated, 
                           int *partition_handle_size, 
                           int *err )
{
  std::vector<ParallelComm*> pcomms;
  ErrorCode rval = ParallelComm::get_all_pcomm( MBI, pcomms );
  CHKERR(rval,"Internal error retreiving PComms");
  
  std::vector<ParallelComm*>::iterator i;
  int count = 0;
  for (i = pcomms.begin(); i != pcomms.end(); ++i)
    if ((*i)->get_partitioning())
      ++count;
  ALLOC_CHECK_ARRAY_NOFAIL( partition_handle, count );
  
  *partition_handle_size = 0;
  for (i = pcomms.begin(); i != pcomms.end(); ++i)
    if ((*i)->get_partitioning())
      (*partition_handle)[(*partition_handle_size)++] 
        = itaps_cast<iMeshP_PartitionHandle>((*i)->get_partitioning());

  RETURN (iBase_SUCCESS );
}

void iMeshP_getNumGlobalParts( iMesh_Instance instance,
                               const iMeshP_PartitionHandle partition_handle,
                               int *num_global_part, 
                               int *err )
{
  ParallelComm* pcomm = PCOMM;
  if (!pcomm) 
    ERROR (iBase_FAILURE,"No PComm");
  
  ErrorCode rval = pcomm->get_global_part_count( *num_global_part );
  CHKERR (rval,"PComm::get_global_part_count failed");
  RETURN(iBase_SUCCESS);
}

void iMeshP_getNumLocalParts(iMesh_Instance instance,
                          const iMeshP_PartitionHandle partition_handle,
                          int *num_local_part, 
                          int *err)
{
  ParallelComm* pcomm = PCOMM;
  if (!pcomm) 
    ERROR (iBase_FAILURE,"No PComm");
  
  *num_local_part = pcomm->partition_sets().size();
  RETURN (iBase_SUCCESS);
}

void iMeshP_getLocalParts( iMesh_Instance instance,
                           const iMeshP_PartitionHandle partition_handle,
                           iMeshP_PartHandle **part_handles,
                           int *part_handles_allocated,
                           int *part_handles_size,
                           int *err )
{
  ParallelComm* pcomm = PCOMM;
  if (!pcomm)
    ERROR (iBase_FAILURE,"No PComm");
  
  RANGE_TO_ITAPS_ARRAY( pcomm->partition_sets(), part_handles );
  RETURN (iBase_SUCCESS);
}

void iMeshP_getRankOfPart( iMesh_Instance instance,
                           const iMeshP_PartitionHandle partition_handle,
                           const iMeshP_Part part_id,
                           int *rank,
                           int *err )
{
  int junk1 = 1, junk2 = 1;
  iMeshP_getRankOfPartArr( instance, partition_handle, &part_id,
                           1, &rank, &junk1, &junk2, err );
}

void iMeshP_getRankOfPartArr( iMesh_Instance instance,
                              const iMeshP_PartitionHandle partition_handle,
                              const iMeshP_Part *part_ids,
                              const int part_ids_size,
                              int **rank, 
                              int *rank_allocated, 
                              int *rank_size,
                              int *err )
{
  ParallelComm* pcomm = PCOMM;
  if (!pcomm)
    ERROR (iBase_FAILURE,"No PComm");
  
  ALLOC_CHECK_ARRAY( rank, part_ids_size );
  ErrorCode rval = MB_SUCCESS;
  for (int i = 0; i < part_ids_size; ++i) {
    rval = pcomm->get_part_owner( part_ids[i], (*rank)[i] );
    CHKERR(rval,"PComm::get_part_owner failed");
  }
  KEEP_ARRAY(rank);
  RETURN(iBase_SUCCESS);
}

void iMeshP_getNumOfTypeAll( iMesh_Instance instance,
                             const iMeshP_PartitionHandle partition_handle,
                             const iBase_EntitySetHandle entity_set_handle,
                             const int entity_type, 
                             int *num_type, 
                             int *err )
{
  ParallelComm* pcomm = PCOMM;
  if (!pcomm)
    ERROR (iBase_FAILURE,"No PComm");
  
  Range entities;
  ErrorCode rval = get_entities( MBI,
                                   itaps_cast<EntityHandle>(entity_set_handle),
                                   entity_type,
                                   iMesh_ALL_TOPOLOGIES,
                                   entities );
  int count = 0;
  if (MB_SUCCESS == rval)
    rval = count_owned( pcomm, entities, count );
  
  int vals[2] = { count, rval }, sums[2];
  int ierr = MPI_Allreduce( vals, sums, 2, MPI_INT, MPI_SUM, pcomm->proc_config().proc_comm() );
  assert(iBase_SUCCESS == 0);
  if (ierr || sums[1])
    RETURN (iBase_FAILURE);
  
  *num_type = sums[0];
  RETURN (iBase_SUCCESS);
}

void iMeshP_getNumOfTopoAll( iMesh_Instance instance,
                             const iMeshP_PartitionHandle partition_handle,
                             const iBase_EntitySetHandle entity_set_handle,
                             const int entity_topology, 
                             int *num_topo, 
                             int *err )
{
  ParallelComm* pcomm = PCOMM;
  if (!pcomm)
    ERROR (iBase_FAILURE,"No PComm");
  
  Range entities;
  ErrorCode rval = get_entities( MBI,
                                   itaps_cast<EntityHandle>(entity_set_handle),
                                   iBase_ALL_TYPES,
                                   entity_topology,
                                   entities );
  int count = 0;
  if (MB_SUCCESS == rval)
    rval = count_owned( pcomm, entities, count );
  
  int vals[2] = { count, rval }, sums[2];
  int ierr = MPI_Allreduce( vals, sums, 2, MPI_INT, MPI_SUM, pcomm->proc_config().proc_comm() );
  assert(iBase_SUCCESS == 0);
  if (ierr || sums[1])
    RETURN (iBase_FAILURE);
  
  *num_topo = sums[0];
  RETURN (iBase_SUCCESS);
}

void iMeshP_createPart( iMesh_Instance instance,
                        iMeshP_PartitionHandle partition_handle,
                        iMeshP_PartHandle *part_handle,
                        int *err )
{
  ParallelComm* pcomm = PCOMM;
  if (!pcomm)
    ERROR (iBase_FAILURE,"No PComm");
  
  EntityHandle h;
  ErrorCode rval = pcomm->create_part( h );
  CHKERR(rval,"Part creation failed");
  *part_handle = itaps_cast<iMeshP_PartHandle>(h);
}

void iMeshP_destroyPart( iMesh_Instance instance,
                         iMeshP_PartitionHandle partition_handle,
                         iMeshP_PartHandle part_handle,
                         int *err )
{
  ParallelComm* pcomm = PCOMM;
  if (!pcomm)
    ERROR (iBase_FAILURE,"No PComm");
  
  ErrorCode rval = pcomm->destroy_part( itaps_cast<EntityHandle>(part_handle) );
  CHKERR(rval,"Part destruction failed");
  RETURN(iBase_SUCCESS);
}

void iMeshP_getNumPartNbors( iMesh_Instance instance,
                             iMeshP_PartitionHandle partition_handle,
                             iMeshP_PartHandle part_handle,
                             int entity_type,
                             int *num_part_nbors,
                             int *err )
{
  int junk1 = 1, junk2 = 1;
  iMeshP_getNumPartNborsArr( instance, partition_handle,
                             &part_handle, 1, entity_type,
                             &num_part_nbors, &junk1, &junk2,
                             err );
}

void iMeshP_getNumPartNborsArr( iMesh_Instance instance,
                                const iMeshP_PartitionHandle partition_handle,
                                const iMeshP_PartHandle *part_handles,
                                const int part_handles_size,
                                int entity_type,
                                int **num_part_nbors,
                                int *num_part_nbors_allocated,
                                int *num_part_nbors_size,
                                int *err )
{
  ParallelComm* pcomm = PCOMM;
  if (!pcomm) 
    ERROR (iBase_FAILURE,"No PComm");

  ALLOC_CHECK_ARRAY( num_part_nbors, part_handles_size );
  
  int n, neighbors[MAX_SHARING_PROCS];
  ErrorCode rval;
  for (int i = 0; i < part_handles_size; ++i) {
    EntityHandle h = itaps_cast<EntityHandle>(part_handles[i]);
    rval = pcomm->get_part_neighbor_ids( h, neighbors, n ); 
    CHKERR(rval,"error getting neighbor ids");
    (*num_part_nbors)[i] = n;
  }
  
  KEEP_ARRAY(num_part_nbors);
  RETURN(iBase_SUCCESS);
}


void iMeshP_getPartNbors( iMesh_Instance instance,
                          const iMeshP_PartitionHandle partition_handle,
                          const iMeshP_PartHandle part_handle,
                          int entity_type,
                          int *num_part_nbors,
                          iMeshP_Part **nbor_part_ids,
                          int *nbor_part_ids_allocated,
                          int *nbor_part_ids_size,
                          int *err )
{
  int junk1 = 1, junk2 = 1;
  iMeshP_getPartNborsArr( instance, partition_handle, 
                          &part_handle, 1, entity_type, 
                          &num_part_nbors, &junk1, &junk2,
                          nbor_part_ids, nbor_part_ids_allocated, 
                          nbor_part_ids_size, err );
}

void iMeshP_getPartNborsArr( iMesh_Instance instance,
                             const iMeshP_PartitionHandle partition_handle,
                             const iMeshP_PartHandle *part_handles,
                             const int part_handles_size,
                             int entity_type,
                             int **num_part_nbors,
                             int *num_part_nbors_allocated,
                             int *num_part_nbors_size,
                             iMeshP_Part **nbor_part_ids,
                             int *nbor_part_ids_allocated,
                             int *nbor_part_ids_size,
                             int *err ) 
{
  ParallelComm* pcomm = PCOMM;
  if (!pcomm) 
    ERROR (iBase_FAILURE,"No PComm");

  ALLOC_CHECK_ARRAY( num_part_nbors, part_handles_size );
  
  std::vector<int> all_neighbors;
  int n, pnbor[MAX_SHARING_PROCS];
  ErrorCode rval;
  for (int i = 0; i < part_handles_size; ++i) {
    EntityHandle h = itaps_cast<EntityHandle>(part_handles[i]);
    rval = pcomm->get_part_neighbor_ids( h, pnbor, n ); 
    CHKERR(rval,"error getting neighbor ids");
    (*num_part_nbors)[i] = n;
    std::copy( pnbor, pnbor+n, std::back_inserter(all_neighbors) );
  }
  
  ALLOC_CHECK_ARRAY_NOFAIL( nbor_part_ids, all_neighbors.size() );
  memcpy( *nbor_part_ids, &all_neighbors[0], sizeof(int)*all_neighbors.size() );
  
  KEEP_ARRAY(num_part_nbors);
  RETURN(iBase_SUCCESS);
}

static ErrorCode get_boundary_entities( ParallelComm* pcomm,
                                   EntityHandle part_handle,
                                   int entity_type,
                                   int entity_topology,
                                   int adj_part_id,
                                   Range& entities_out)
{
  int* adj_part_id_ptr = (adj_part_id == iMeshP_ALL_PARTS) ? 0 : &adj_part_id;
  
  Range iface_sets;
  ErrorCode rval = pcomm->get_interface_sets( 
                         itaps_cast<EntityHandle>(part_handle),
                         iface_sets, adj_part_id_ptr ); 
  if (MB_SUCCESS != rval)
    return rval;
  
  for (Range::iterator i = iface_sets.begin(); i != iface_sets.end(); ++i) {
    rval = get_entities( pcomm->get_moab(), *i, entity_type, entity_topology, entities_out );
    if (MB_SUCCESS != rval)
      return rval;
  }
  
  return MB_SUCCESS;
}

void iMeshP_getNumPartBdryEnts( iMesh_Instance instance,
                                const iMeshP_PartitionHandle partition_handle,
                                const iMeshP_PartHandle part_handle, 
                                const int entity_type, 
                                const int entity_topology, 
                                const iMeshP_Part target_part_id, 
                                int *num_entities, 
                                int *err )
{
  Range entities;
  ErrorCode rval = get_boundary_entities( PCOMM,
                                            itaps_cast<EntityHandle>(part_handle),
                                            entity_type,
                                            entity_topology,
                                            target_part_id,
                                            entities );
  CHKERR(rval,"failed to get boundary entities");
  *num_entities = entities.size();
  RETURN(iBase_SUCCESS);
}

void iMeshP_getPartBdryEnts( iMesh_Instance instance,
                             const iMeshP_PartitionHandle partition_handle,
                             const iMeshP_PartHandle part_handle,
                             const int entity_type,
                             const int entity_topology,
                             const iMeshP_Part target_part_id,
                             iBase_EntityHandle **entity_handles,
                             int *entity_handles_allocated,
                             int *entity_handles_size,
                             int *err )
{
  Range entities;
  ErrorCode rval = get_boundary_entities( PCOMM,
                                            itaps_cast<EntityHandle>(part_handle),
                                            entity_type,
                                            entity_topology,
                                            target_part_id,
                                            entities );
  CHKERR(rval,"failed to get boundary entities");
  RANGE_TO_ITAPS_ARRAY( entities, entity_handles );
  RETURN(iBase_SUCCESS);
}

void iMeshP_initPartBdryEntIter( iMesh_Instance instance,
                                 const iMeshP_PartitionHandle partition_handle,
                                 const iMeshP_PartHandle part_handle,
                                 const int entity_type,
                                 const int entity_topology,
                                 const iMeshP_Part nbor_part_id,
                                 iMesh_EntityIterator* entity_iterator,
                                 int* err )
{
  Range entities;
  ErrorCode rval = get_boundary_entities( PCOMM,
                                            itaps_cast<EntityHandle>(part_handle),
                                            entity_type,
                                            entity_topology,
                                            nbor_part_id,
                                            entities ); 
  CHKERR(rval,"error getting boundary");
  *entity_iterator = create_itaps_iterator( entities );
  RETURN( entity_iterator ? iBase_SUCCESS : iBase_FAILURE );
}

void iMeshP_initPartBdryEntArrIter( iMesh_Instance instance,
                                    const iMeshP_PartitionHandle partition_handle,
                                    const iMeshP_PartHandle part_handle,
                                    const int entity_type,
                                    const int entity_topology,
                                    const int array_size,
                                    const iMeshP_Part nbor_part_id,
                                    iMesh_EntityArrIterator* entity_iterator,
                                    int* err )
{
  Range entities;
  ErrorCode rval = get_boundary_entities( PCOMM,
                                            itaps_cast<EntityHandle>(part_handle),
                                            entity_type,
                                            entity_topology,
                                            nbor_part_id,
                                            entities ); 
  CHKERR(rval,"error getting boundary");
  *entity_iterator = (iMesh_EntityArrIterator) create_itaps_iterator( entities, array_size );
  RETURN( entity_iterator ? iBase_SUCCESS : iBase_FAILURE );
}



void iMeshP_getNumOfType( iMesh_Instance instance,
                          const iMeshP_PartitionHandle ,
                          const iMeshP_PartHandle part_handle,
                          const iBase_EntitySetHandle entity_set_handle,
                          const int entity_type,
                          int *num_type,
                          int *err )
{
  Range r;
  set_intersection_query( instance, part_handle, entity_set_handle, 
                          entity_type, iMesh_ALL_TOPOLOGIES, r, err );
  *num_type = r.size();
}

void iMeshP_getNumOfTopo( iMesh_Instance instance,
                          const iMeshP_PartitionHandle partition_handle,
                          const iMeshP_PartHandle part_handle,
                          const iBase_EntitySetHandle entity_set_handle,
                          const int entity_topology,
                          int *num_topo,
                          int *err )
{
  Range r;
  set_intersection_query( instance, part_handle, entity_set_handle, 
                          iBase_ALL_TYPES, entity_topology, r, err );
  *num_topo = r.size();
}

void iMeshP_getAdjEntIndices(iMesh_Instance instance,
                             iMeshP_PartitionHandle partition,
                             iMeshP_PartHandle part,
                             iBase_EntitySetHandle entity_set_handle,
                             int entity_type_requestor,
                             int entity_topology_requestor,
                             int entity_type_requested,
                             iBase_EntityHandle** entity_handles,
                             int* entity_handles_allocated,
                             int* entity_handles_size,
                             iBase_EntityHandle** adj_entity_handles,
                             int* adj_entity_handles_allocated,
                             int* adj_entity_handles_size,
                             int** adj_entity_indices,
                             int* adj_entity_indices_allocated,
                             int* adj_entity_indices_size,
                             int** offset,
                             int* offset_allocated,
                             int* offset_size,
                             int *err)
{
  const int allocated_entity_handles = (*entity_handles_allocated == 0);
  const int allocated_indices = (*adj_entity_indices_allocated == 0);
  const int allocated_offset = (*offset_allocated == 0);

  // get source entities
  iMeshP_getEntities( instance, 
                      partition, part,
                      entity_set_handle,
                      entity_type_requestor, 
                      entity_topology_requestor,
                      entity_handles,
                      entity_handles_allocated,
                      entity_handles_size,
                      err );
  if (iBase_SUCCESS != *err)
    return;

  // get adjacencies
  iBase_EntityHandle* all_adj_handles = 0;
  int size = 0, alloc = 0;
  iMesh_getEntArrAdj( instance,
                      *entity_handles, *entity_handles_size,
                      entity_type_requested,
                      &all_adj_handles, &alloc, &size,
                      offset, offset_allocated, offset_size,
                      err );
  if (*err != iBase_SUCCESS) {
    if (allocated_entity_handles) {
      free( *entity_handles );
      *entity_handles = 0;
      *entity_handles_allocated = 0;
    }
    return;
  }

  // allocate or check size of adj_entity_indices
  *adj_entity_indices_size = size;
  if (allocated_indices) {
    MALLOC(*adj_entity_indices, sizeof(iBase_EntityHandle)*size, int*);
//    *adj_entity_indices = (int*)malloc(sizeof(iBase_EntityHandle)*size);
    if (!*adj_entity_indices) 
      *err = iBase_MEMORY_ALLOCATION_FAILED;
    else
      *adj_entity_indices_allocated = size;
  }
  else if (*adj_entity_indices_allocated < size) {
    *err = iBase_BAD_ARRAY_DIMENSION;
  }
  if (iBase_SUCCESS != *err) {
    free( all_adj_handles );
    if (allocated_entity_handles) {
      free( *entity_handles );
      *entity_handles = 0;
      *entity_handles_allocated = 0;
    }
    if (allocated_offset) {
      free( *offset );
      *offset = 0;
      *offset_allocated = 0;
    }
    return;
  }

  // Now create an array of unique sorted handles from all_adj_handles.
  // We need to create a copy because we still need all_adj_handles.  We
  // will eventually need to copy the resulting unique list into 
  // adj_entity_handles, so if adj_entity_handles is already allocated and
  // of sufficient size, use it rather than allocating another temporary.
  iBase_EntityHandle* unique_adj = 0;
  if (*adj_entity_handles_allocated >= size) {
    unique_adj = *adj_entity_handles;
  }
  else {
    MALLOC(unique_adj, sizeof(iBase_EntityHandle) * size, iBase_EntityHandle*);
  }
  std::copy( all_adj_handles, all_adj_handles+size, unique_adj );
  std::sort( unique_adj, unique_adj + size );
  *adj_entity_handles_size = std::unique( unique_adj, unique_adj + size ) - unique_adj;

  // If we created a temporary array for unique_adj rather than using
  // already allocated space in adj_entity_handles, allocate adj_entity_handles
  // and copy the unique handle list into it
  if (*adj_entity_handles != unique_adj) {
    if (!*adj_entity_handles_allocated) {
      MALLOC(*adj_entity_handles, 
             sizeof(iBase_EntityHandle) * *adj_entity_handles_size, iBase_EntityHandle*);
      if (!*adj_entity_handles)
        *err = iBase_MEMORY_ALLOCATION_FAILED;
      else
        *adj_entity_handles_allocated = *adj_entity_handles_size;
    }
    else if (*adj_entity_handles_allocated < *adj_entity_handles_size) 
      *err = iBase_BAD_ARRAY_DIMENSION;
    if (iBase_SUCCESS != *err) {
      free( unique_adj );
      free( all_adj_handles );
      if (allocated_entity_handles) {
        free( *entity_handles );
        *entity_handles = 0;
        *entity_handles_allocated = 0;
      }
      if (allocated_offset) {
        free( *offset );
        *offset = 0;
        *offset_allocated = 0;
      }
      if (allocated_indices) {
        free( *adj_entity_indices );
        *adj_entity_indices = 0;
        *adj_entity_indices_allocated = 0;
      }
      return;
    }

    std::copy( unique_adj, unique_adj + *adj_entity_handles_size, *adj_entity_handles );
    free( unique_adj );
    unique_adj = *adj_entity_handles;
  }

  // convert from adjacency list to indices into unique_adj
  for (int i = 0; i < *adj_entity_indices_size; ++i)
    (*adj_entity_indices)[i] = std::lower_bound( unique_adj, 
      unique_adj + *adj_entity_handles_size, all_adj_handles[i] ) - unique_adj;
  free( all_adj_handles );
}

void iMeshP_getEntities( iMesh_Instance instance,
                         const iMeshP_PartitionHandle ,
                         const iMeshP_PartHandle part_handle,
                         const iBase_EntitySetHandle entity_set_handle,
                         const int entity_type,
                         const int entity_topology,
                         iBase_EntityHandle** entity_handles,
                         int* entity_handles_allocated,
                         int* entity_handles_size,
                         int *err )
{
  Range r;
  set_intersection_query( instance, part_handle, entity_set_handle, 
                          entity_type, entity_topology, r, err );
  if (iBase_SUCCESS != *err)
    return;
  
  RANGE_TO_ITAPS_ARRAY( r, entity_handles );
  RETURN(iBase_SUCCESS);
}

void iMeshP_getAdjEntities( iMesh_Instance instance,
                            const iMeshP_PartitionHandle partition_handle,
                            const iMeshP_PartHandle part_handle,
                            const iBase_EntitySetHandle entity_set_handle,
                            const int entity_type_requestor,
                            const int entity_topology_requestor,
                            const int entity_type_requested,
                            iBase_EntityHandle** adj_entity_handles,
                            int* adj_entity_handles_allocated,
                            int* adj_entity_handles_size,
                            int** offset,
                            int* offset_allocated,
                            int* offset_size,
                            int** in_entity_set,
                            int* in_entity_set_allocated,
                            int* in_entity_set_size,
                            int *err )
{
  ErrorCode rval;
  Range r;
  set_intersection_query( instance, part_handle, entity_set_handle,
                           entity_type_requestor, entity_topology_requestor,
                           r, err );
  if (iBase_SUCCESS != *err)
    return;
  
    // count adjacencies
  std::vector<EntityHandle> tmp_storage;
  int num_adj = 0;
  int num_conn;
  const EntityHandle* conn_ptr;
  for (Range::iterator i = r.begin(); i != r.end(); ++i)  {
    if (entity_type_requested || TYPE_FROM_HANDLE(*i) == MBPOLYHEDRON) {
      tmp_storage.clear();
      rval = MBI->get_adjacencies( &*i, 1, entity_type_requested, false, tmp_storage );
      CHKERR(rval,"get_adjacencies failed");
      num_adj += tmp_storage.size();
    }
    else {
      rval = MBI->get_connectivity( *i, conn_ptr, num_conn, false, &tmp_storage );
      CHKERR(rval,"get_connectivity failed");
      num_adj += num_conn;
    }
  }
  
    // get adjacencies
  ALLOC_CHECK_ARRAY( adj_entity_handles, num_adj );
  ALLOC_CHECK_ARRAY( offset, r.size() );
  int arr_pos = 0;
  int* offset_iter = *offset;
  for (Range::iterator i = r.begin(); i != r.end(); ++i)  {
    *offset_iter = arr_pos; 
    ++offset_iter;

    tmp_storage.clear();
    rval = MBI->get_adjacencies( &*i, 1, entity_type_requested, false, tmp_storage );
    CHKERR(rval,"get_adjacencies failed");
    for (std::vector<EntityHandle>::iterator j = tmp_storage.begin(); j != tmp_storage.end(); ++j) {
      (*adj_entity_handles)[arr_pos] = itaps_cast<iBase_EntityHandle>(*j);
      ++arr_pos;
    }
  }

    // get in_entity_set
  iMesh_isEntArrContained( instance, 
                           entity_set_handle, 
                           *adj_entity_handles,
                           *adj_entity_handles_size,
                           in_entity_set,
                           in_entity_set_allocated,
                           in_entity_set_size,
                           err );
  
  if (iBase_SUCCESS == *err) {
    KEEP_ARRAY(adj_entity_handles);
    KEEP_ARRAY(offset);
  }
}

void iMeshP_initEntIter( iMesh_Instance instance,
                         const iMeshP_PartitionHandle partition_handle,
                         const iMeshP_PartHandle part_handle,
                         const iBase_EntitySetHandle entity_set_handle,
                         const int requested_entity_type,
                         const int requested_entity_topology,
                         iMesh_EntityIterator* entity_iterator,
                         int *err )
{
  Range r;
  set_intersection_query( instance, part_handle, entity_set_handle,
                           requested_entity_type, requested_entity_topology,
                           r, err );
  if (iBase_SUCCESS != *err)
    return;
  
  *entity_iterator = create_itaps_iterator( r );
  RETURN (iBase_SUCCESS);
}

void iMeshP_initEntArrIter( iMesh_Instance instance,
                            const iMeshP_PartitionHandle partition_handle,
                            const iMeshP_PartHandle part_handle,
                            const iBase_EntitySetHandle entity_set_handle,
                            const int requested_entity_type,
                            const int requested_entity_topology,
                            const int requested_array_size,
                            iMesh_EntityArrIterator* entArr_iterator,
                            int *err )
{
  Range r;
  set_intersection_query( instance, part_handle, entity_set_handle,
                           requested_entity_type, requested_entity_topology,
                           r, err );
  if (iBase_SUCCESS != *err)
    return;
  
  *entArr_iterator = (iMesh_EntityArrIterator)create_itaps_iterator( r, requested_array_size );
  RETURN (iBase_SUCCESS);
}

void iMeshP_getEntOwnerPart( iMesh_Instance instance,
                             const iMeshP_PartitionHandle partition_handle,
                             const iBase_EntityHandle entity_handle,
                             iMeshP_Part *part_id,
                             int* err )
{ 
  int junk1 = 1, junk2 = 1;
  iMeshP_getEntOwnerPartArr( instance, partition_handle, &entity_handle, 1,
                             &part_id, &junk1, &junk2, err );
}

void iMeshP_getEntOwnerPartArr( iMesh_Instance instance,
                                const iMeshP_PartitionHandle partition_handle,
                                const iBase_EntityHandle *entity_handles,
                                const int entity_handles_size,
                                iMeshP_Part **part_ids,
                                int *part_ids_allocated,
                                int *part_ids_size,
                                int* err )
{
  ParallelComm* pcomm = PCOMM;
  if (!pcomm) 
    ERROR (iBase_FAILURE,"No PComm");
  
  int id;
  ALLOC_CHECK_ARRAY( part_ids, entity_handles_size );
  ErrorCode rval = MB_SUCCESS;
  for (int i = 0; i < entity_handles_size; ++i) {
    EntityHandle h = itaps_cast<EntityHandle>(entity_handles[i]);
    rval = pcomm->get_owning_part( h, id );
    (*part_ids)[i] = id;
    CHKERR(rval,"Failet get part owner");
  }
  KEEP_ARRAY(part_ids);
  RETURN(iBase_SUCCESS);
}
  
void iMeshP_isEntOwner( iMesh_Instance instance,
                        const iMeshP_PartitionHandle partition_handle,
                        const iMeshP_PartHandle part_handle,
                        const iBase_EntityHandle entity_handle,
                        int* is_owner,
                        int *err )
{
  int junk1 = 1, junk2 = 1;
  iMeshP_isEntOwnerArr( instance, partition_handle,
                        part_handle, &entity_handle, 1,
                        &is_owner, &junk1, &junk2,
                        err );
}

void iMeshP_isEntOwnerArr( iMesh_Instance instance,
                           const iMeshP_PartitionHandle partition_handle,
                           const iMeshP_PartHandle part_handle,
                           const iBase_EntityHandle *entity_handles,
                           const int entity_handles_size,
                           int** is_owner,
                           int* is_owner_allocated,
                           int* is_owner_size,
                           int *err )
{
  ErrorCode rval;
  ParallelComm* pcomm = PCOMM;
  if (!pcomm) 
    ERROR (iBase_FAILURE,"No PComm");
  
  int id;
  rval = pcomm->get_part_id( itaps_cast<EntityHandle>(part_handle), id );
  CHKERR(rval,"error getting part id");
  
  ALLOC_CHECK_ARRAY( is_owner, entity_handles_size );
  *is_owner_size = entity_handles_size;
  
  int owner;
  for (int i = 0; i < entity_handles_size; ++i) {
    rval = pcomm->get_owner( itaps_cast<EntityHandle>(entity_handles[i]), owner );
    CHKERR(rval,"error getting owner");
    (*is_owner)[i] = (owner == id);
  }
  
  KEEP_ARRAY(is_owner);
  RETURN(iBase_SUCCESS);
}

void iMeshP_getEntStatus(iMesh_Instance instance,
                         /*in*/ const iMeshP_PartitionHandle partition_handle,
                         /*in*/ const iMeshP_PartHandle part_handle, 
                         /*in*/ const iBase_EntityHandle entity_handle, 
                         /*out*/ int* par_status, // Values=INTERNAL,BOUNDARY,GHOST
                         int *err) 
{
  int junk1 = 1, junk2 = 1;
  iMeshP_getEntStatusArr( instance, partition_handle,
                          part_handle, &entity_handle, 1,
                          &par_status, &junk1, &junk2,
                          err );
}

void iMeshP_getEntStatusArr(iMesh_Instance instance,
                            /*in*/ const iMeshP_PartitionHandle partition_handle,
                            /*in*/ const iMeshP_PartHandle part_handle, 
                            /*in*/ const iBase_EntityHandle *entity_handles, 
                            /*in*/ const int entity_handles_size, 
                            /*inout*/ int** par_status, // Values=INTERNAL,BOUNDARY,GHOST
                            /*inout*/ int* par_status_allocated, 
                            /*inout*/ int* par_status_size, 
                            int *err) 
{
  ParallelComm* pcomm = PCOMM;
  if (!pcomm) 
    ERROR (iBase_FAILURE,"No PComm");

  std::vector<unsigned char> pstatus(entity_handles_size);
  ErrorCode result = MBI->tag_get_data(pcomm->pstatus_tag(), 
                                         itaps_cast<const EntityHandle*>(entity_handles), 
                                         entity_handles_size,
                                         &pstatus[0]); 
  CHKERR(result,"error getting pstatus_tag");

  ALLOC_CHECK_ARRAY( par_status, entity_handles_size );
  for (int i = 0; i < entity_handles_size; i++) {
    if (!pstatus[i]) 
      (*par_status)[i] = iMeshP_INTERNAL;
    else if (pstatus[i] & PSTATUS_GHOST) 
      (*par_status)[i] = iMeshP_GHOST;
    else if (pstatus[i] & PSTATUS_INTERFACE)
      (*par_status)[i] = iMeshP_BOUNDARY;
  }

  KEEP_ARRAY(par_status);
  RETURN(iBase_SUCCESS);
}

void iMeshP_getNumCopies( iMesh_Instance instance,
                          const iMeshP_PartitionHandle partition_handle,
                          const iBase_EntityHandle entity_handle,
                          int *num_copies_ent,
                          int *err )
{
  ParallelComm* pcomm = PCOMM;
  if (!pcomm) 
    ERROR (iBase_FAILURE,"No PComm");

  int ids[MAX_SHARING_PROCS];
  ErrorCode rval = pcomm->get_sharing_parts( 
                              itaps_cast<EntityHandle>(entity_handle),
                              ids, *num_copies_ent );
  CHKERR(rval,"ParallelComm::get_sharing_parts failed");
  RETURN(iBase_SUCCESS);
}

void iMeshP_getCopyParts( iMesh_Instance instance,
                          const iMeshP_PartitionHandle partition_handle,
                          const iBase_EntityHandle entity_handle,
                          iMeshP_Part **part_ids,
                          int *part_ids_allocated,
                          int *part_ids_size,
                          int *err )
{
  ParallelComm* pcomm = PCOMM;
  if (!pcomm) 
    ERROR (iBase_FAILURE,"No PComm");

  int ids[MAX_SHARING_PROCS], num_ids;
  ErrorCode rval = pcomm->get_sharing_parts( 
                              itaps_cast<EntityHandle>(entity_handle),
                              ids, num_ids ); 
  CHKERR(rval,"ParallelComm::get_sharing_parts failed");
  ALLOC_CHECK_ARRAY_NOFAIL( part_ids, num_ids );
  std::copy( ids, ids+num_ids, *part_ids );
  RETURN (iBase_SUCCESS);
}



void iMeshP_getCopies( iMesh_Instance instance,
                       const iMeshP_PartitionHandle partition_handle, 
                       const iBase_EntityHandle entity_handle, 
                       iMeshP_Part **part_ids, 
                       int *part_ids_allocated, 
                       int *part_ids_size, 
                       iBase_EntityHandle **copies_entity_handles, 
                       int *copies_entity_handles_allocated, 
                       int *copies_entity_handles_size, 
                       int *err )
{
  ParallelComm* pcomm = PCOMM;
  if (!pcomm) 
    ERROR (iBase_FAILURE,"No PComm");

  int ids[MAX_SHARING_PROCS], num_ids;
  EntityHandle handles[MAX_SHARING_PROCS];
  ErrorCode rval = pcomm->get_sharing_parts( 
                              itaps_cast<EntityHandle>(entity_handle),
                              ids, num_ids, handles ); 
  CHKERR(rval,"ParallelComm::get_sharing_parts failed");
  ALLOC_CHECK_ARRAY_NOFAIL( part_ids, num_ids );
  ALLOC_CHECK_ARRAY_NOFAIL( copies_entity_handles, num_ids );
  for (int i = 0; i < num_ids; ++i) {
    (*part_ids)[i] = ids[i];
    (*copies_entity_handles)[i] = itaps_cast<iBase_EntityHandle>(handles[i]);
  }
  RETURN (iBase_SUCCESS);
}

void iMeshP_getCopyOnPart( iMesh_Instance instance,
                           const iMeshP_PartitionHandle partition_handle,
                           const iBase_EntityHandle entity_handle,
                           const iMeshP_Part part_id,
                           iBase_EntityHandle* copy_entity_handle,
                           int *err )
{
  ParallelComm* pcomm = PCOMM;
  if (!pcomm) 
    ERROR (iBase_FAILURE,"No PComm");

  int ids[MAX_SHARING_PROCS], num_ids;
  EntityHandle handles[MAX_SHARING_PROCS];
  ErrorCode rval = pcomm->get_sharing_parts( 
                              itaps_cast<EntityHandle>(entity_handle),
                              ids, num_ids, handles ); 
  CHKERR(rval,"ParallelComm::get_sharing_parts failed");
  int idx = std::find( ids, ids+num_ids, part_id ) - ids;
  if (idx == num_ids)
    RETURN (iBase_FAILURE);
  
  *copy_entity_handle = itaps_cast<iBase_EntityHandle>(handles[idx]);
  RETURN (iBase_SUCCESS);
}


void iMeshP_getOwnerCopy( iMesh_Instance instance,
                          const iMeshP_PartitionHandle partition_handle,
                          const iBase_EntityHandle entity_handle,
                          iMeshP_Part *owner_part_id,
                          iBase_EntityHandle *owner_entity_handle,
                          int *err )
{ 
  ParallelComm* pcomm = PCOMM;
  if (!pcomm) 
    ERROR (iBase_FAILURE,"No PComm");

  int id;
  EntityHandle h;
  ErrorCode rval = pcomm->get_owning_part( 
                          itaps_cast<EntityHandle>(entity_handle),
                          id, &h );  
  CHKERR(rval,"Failed to get owner");
  *owner_part_id = id;
  *owner_entity_handle = itaps_cast<iBase_EntityHandle>(h);
  RETURN(iBase_SUCCESS);
}

void iMeshP_waitForAnyRequest( iMesh_Instance instance,
                               iMeshP_PartitionHandle partition_handle,
                               iMeshP_RequestHandle *requests,
                               int request_size,
                               int* index,
                               int* err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_waitForAllRequests( iMesh_Instance instance,
                                iMeshP_PartitionHandle partition_handle,
                                iMeshP_RequestHandle *requests,
                                int requests_size,
                                int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_waitForRequestEnt(
            iMesh_Instance instance,
            const iMeshP_PartitionHandle partition,
            iMeshP_RequestHandle request,
            iBase_EntityHandle **out_entities,
            int *out_entities_allocated,
            int *out_entities_size,
            int *err)
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_testRequest( iMesh_Instance instance,
                  const iMeshP_PartitionHandle partition_handle,
                  iMeshP_RequestHandle req,
                  int *complete,
                  int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_pollForRequests( iMesh_Instance instance,
                             iMeshP_PartitionHandle partition_handle,
                             iMeshP_RequestHandle **requests_completed,
                             int *requests_allocated,
                             int *requests_size,
                             int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_exchEntArrToPartsAll( iMesh_Instance instance,
                                  const iMeshP_PartitionHandle partition_handle,
                                  const iBase_EntityHandle *entity_handles,
                                  const int entity_handles_size,
                                  const iMeshP_Part *target_part_ids,
                                  int command_code,
                                  int update_ghost,
                                  iMeshP_RequestHandle *request,
                                  int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_migrateEntity( iMesh_Instance instance,
                           const iMeshP_PartitionHandle partition_handle,
                           const iMeshP_PartHandle part_handle,
                           const iBase_EntityHandle local_entity_handle,
                           iMeshP_RequestHandle *request,
                           int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_updateVtxCoords( iMesh_Instance instance,
                             const iMeshP_PartitionHandle partition_handle,
                             const iBase_EntityHandle local_vertex_handle,
                             int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_replaceOnPartBdry( iMesh_Instance instance,
                               const iMeshP_PartitionHandle partition_handle,
                               const iBase_EntityHandle *old_entities,
                               const int old_entities_size,
                               const iBase_EntityHandle *new_entities,
                               const int new_entities_size,
                               const int *offset,
                               const int offset_size,
                               int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_addGhostOf( iMesh_Instance instance,
                        const iMeshP_PartitionHandle partition_handle,
                        const iMeshP_Part target_part_id,
                        const iBase_EntityHandle entity_to_copy,
                        iMeshP_RequestHandle *request,
                        int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_rmvGhostOf( iMesh_Instance instance,
                        const iMeshP_PartitionHandle partition_handle,
                        const iMeshP_Part target_part_id,
                        const iBase_EntityHandle copy_to_purge,
                        int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_syncMeshAll( iMesh_Instance instance,
                         const iMeshP_PartitionHandle partition_handle,
                         int *err )
{
  FIXME; // for now we only sync vertex coordinates
         // need to update ParallelComm::update_shared_mesh to fix this
  
  ParallelComm* pcomm = PCOMM;
  ErrorCode rval = pcomm->update_shared_mesh();
  CHKERR(rval,"update failed");
  RETURN (iBase_SUCCESS);
}
		            
void iMeshP_pushTags( iMesh_Instance instance,
                      const iMeshP_PartitionHandle partition_handle,
                      iBase_TagHandle source_tag,
                      iBase_TagHandle dest_tag,
                      int entity_type,
                      int entity_topo,
                      int *err )
{
  ParallelComm* pcomm = PCOMM;
  DimensionPair types;
  if (entity_topo != iMesh_ALL_TOPOLOGIES)
    types.first = types.second = mb_topology_table[entity_topo];
  else if (entity_type != iBase_ALL_TYPES)
    types = CN::TypeDimensionMap[entity_type];
  else { 
    types.first = MBVERTEX; 
    types.second = MBENTITYSET; 
    --types.second; 
  }
  
  std::vector<Tag> src_tags(1, itaps_cast<Tag>(source_tag));
  std::vector<Tag> dst_tags(1, itaps_cast<Tag>(dest_tag));
  
  ErrorCode rval;
  Range entities;
  for (EntityType t = types.first; t <= types.second; ++t) {
    rval = MBI->get_entities_by_type_and_tag( 0, t, &src_tags[0], 0, 1, 
                                              entities, Interface::UNION );
    CHKERR(rval,"error getting entities to push");
  }
  
  rval = pcomm->exchange_tags( src_tags, dst_tags, entities );
  CHKERR(rval,"tag data communication failed");
  RETURN (iBase_SUCCESS);
}

void iMeshP_pushTagsEnt( iMesh_Instance instance,
                         const iMeshP_PartitionHandle partition_handle,
                         iBase_TagHandle source_tag,
                         iBase_TagHandle dest_tag,
                         const iBase_EntityHandle *entities,
                         int entities_size,
                         int *err )
{

  Range range;
  const EntityHandle* ents = itaps_cast<const EntityHandle*>(entities);
  std::copy( ents, ents+entities_size, range_inserter(range) );
  
  std::vector<Tag> src_tags(1, itaps_cast<Tag>(source_tag));
  std::vector<Tag> dst_tags(1, itaps_cast<Tag>(dest_tag));
  ParallelComm* pcomm = PCOMM;
  ErrorCode rval = pcomm->exchange_tags( src_tags, dst_tags, range );
  CHKERR(rval,"tag data communication failed");
  RETURN (iBase_SUCCESS);
}

void iMeshP_iPushTags( iMesh_Instance instance,
                       const iMeshP_PartitionHandle partition_handle,
                       iBase_TagHandle source_tag,
                       iBase_TagHandle dest_tag,
                       int entity_type,
                       int entity_topo,
                       iMeshP_RequestHandle *req,
                       int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_iPushTagsEnt( iMesh_Instance instance,
                          const iMeshP_PartitionHandle partition_handle,
                          iBase_TagHandle source_tag,
                          iBase_TagHandle dest_tag,
                          const iBase_EntityHandle *entities,
                          int entities_size,
                          iMeshP_RequestHandle *req,
                          int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_createGhostEntsAll( iMesh_Instance instance,
                                iMeshP_PartitionHandle partition_handle,
                                int ghost_dim,
                                int bridge_dim,
                                int num_layers,
                                int include_copies,
                                int *err )
{
  if (include_copies) {
    FIXME; RETURN(iBase_NOT_SUPPORTED);
  }
  
  ParallelComm* pcomm = PCOMM;
  ErrorCode rval;
  rval = pcomm->exchange_ghost_cells( ghost_dim, bridge_dim, num_layers, true );
  CHKERR(rval,"ghost exchange failed");
  RETURN(iBase_SUCCESS);
}

void iMeshP_deleteGhostEnts( iMesh_Instance instance,
                             iMeshP_PartitionHandle partition_handle,
                             int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_ghostEntInfo( const iMesh_Instance instance,
                          const iMeshP_PartitionHandle partition_handle,
                          int *ghost_rules_allocated,
                          int *ghost_rules_size,
                          int **ghost_dim,
                          int **bridge_dim,
                          int **num_layers,
                          int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

static void append_option( std::string& opt,
                           const char* option,
                           const char* default_value = 0 )
{
  std::string::size_type i;

    // construct a std::string containing at least the leading
    // separator.
  char sep = opt.empty() ? ';' : opt[0];
  
    // make sure that the separator is ok (not contained in the option)
  const char separators[] = ";:,.!@$#%&*^|/\"'\\~%";
  if (strchr(option,sep) || (default_value && (sep == '=' || strchr(default_value,sep)))) {
      // need a new separator.  
    int c, e = sizeof(separators)-1;
    for (c = 0; c < e; ++c) 
      if (!strchr(opt.c_str(),separators[c]) &&
          !strchr(option,separators[c]) &&
          (!default_value || !strchr(default_value,separators[c])))
        break;
    if (c == e) {
      opt.clear();
      return;
    }
    
    i = 0;
    while (std::string::npos != (i = opt.find(sep,i))) 
      opt[i] = separators[c]; 
    sep = separators[c];
  }
  
    // search for the required option
  std::string search(&sep,1);
  search += option;
  const std::string::size_type sl = search.length();
  i = opt.find( search );
  while (i != std::string::npos) {
    std::string::size_type end = i + sl;
    if (end == opt.size() || opt[end] == sep || opt[end] == '=')
      break;
    i = end;
  }
  
    // if string doesn't already contain the option, append it.
  if (i == std::string::npos) {
    opt += search;
    if (default_value) {
      opt += "=";
      opt += default_value;
    }
  }
}

void iMeshP_loadAll( iMesh_Instance instance,
                  const iMeshP_PartitionHandle partition,
                  const iBase_EntitySetHandle entity_set_handle,
                  const char *name,
                  const char *options,
                  int *err,
                  int name_len,
                  int options_len )
{
  ErrorCode rval;

    // create partition set if user didn't give us one.
  EntityHandle partitioning;
  if (partition) {
    partitioning = itaps_cast<EntityHandle>(partition);
  }
  else {
    rval = MBI->create_meshset( MESHSET_SET, partitioning );
    CHKERR(rval,"failed to create meshset");
  }
  
    // get ParallelComm for partition
  MPI_Comm default_comm = MPI_COMM_WORLD;
  ParallelComm* pcomm = ParallelComm::get_pcomm( MBI, partitioning, &default_comm );
  if (!pcomm) {
    RETURN (iBase_FAILURE);
  }

    // add necessary values to options string
  std::string opt( options, options_len );
  append_option( opt, "PARALLEL" );
  append_option( opt, "PARTITION_DISTRIBUTE" );
  append_option( opt, "PARALLEL_RESOLVE_SHARED_ENTS" );
  std::ostringstream id;
  id << pcomm->get_id();
  append_option( opt, "PCOMM", id.str().c_str() );
  
    // load the file
  iMesh_load( instance, entity_set_handle, name, opt.c_str(), err, name_len, opt.length() );
  if (*err) return;
  
  rval = pcomm->collective_sync_partition();
  CHKERR(rval,"collective sync failed");
  RETURN(iBase_SUCCESS);
}

void iMeshP_saveAll( iMesh_Instance instance,
                  const iMeshP_PartitionHandle partition,
                  const iBase_EntitySetHandle entity_set_handle,
                  const char *name,
                  const char *options,
                  int *err,
                  const int name_len,
                  int options_len )
{
  EntityHandle set;
  set = entity_set_handle ? itaps_cast<EntityHandle>(entity_set_handle)
                          : itaps_cast<EntityHandle>(partition);
  iMesh_save( instance, itaps_cast<iBase_EntitySetHandle>(set), 
              name, options, err, name_len, options_len );

}




//  Map from processes to parts:  
//  Given a partition handle and a process rank,
//  return the part handles owned by the process.
//  COMMUNICATION:  None++.
  void iMeshP_getPartsOnRank(iMesh_Instance instance,
                             const iMeshP_PartitionHandle partition_handle,
                             /*in*/    const int rank,
                             /*inout*/ iMeshP_PartHandle **part_handles, 
                             /*inout*/ int *part_handles_allocated, 
                             /*out*/   int *part_handles_size, 
                             int *err) 
  {
    EntityHandle p = itaps_cast<EntityHandle>(partition_handle);
    ParallelComm *pc = ParallelComm::get_pcomm(MBI, p);
    if (!pc) RETURN(iBase_ERROR_MAP[MB_FAILURE]);

    Range part_sets;
  
    ALLOC_CHECK_ARRAY_NOFAIL( part_handles, pc->partition_sets().size() );
    Range::iterator rit;
    int i;
    for (i = 0, rit = pc->partition_sets().begin(); 
         rit != pc->partition_sets().end(); rit++, i++)
      (*part_handles)[i] = itaps_cast<iMeshP_PartHandle>(*rit);
  
    RETURN(iBase_SUCCESS);
  }
    
  void iMeshP_getPartsArrOnRank(iMesh_Instance instance,
                                const iMeshP_PartitionHandle partition_handle,
                                /*in*/    const int *rank,
                                /*in*/    const int rank_size,
                                /*inout*/ iMeshP_PartHandle **part_handles, 
                                /*inout*/ int *part_handles_allocated, 
                                /*out*/   int *part_handles_size, 
                                int *err) 
  {
    EntityHandle p = itaps_cast<EntityHandle>(partition_handle);
    ParallelComm *pc = ParallelComm::get_pcomm(MBI, p);
    if (!pc) RETURN(iBase_ERROR_MAP[MB_FAILURE]);

    if (rank[0] != (int)pc->proc_config().proc_rank() || rank_size > 1) {
      RETURN(iBase_ERROR_MAP[MB_NOT_IMPLEMENTED]);
    }
  
    iMeshP_getPartsOnRank(instance, partition_handle, rank[0],
                          part_handles, part_handles_allocated, part_handles_size,
                          err);
  }

#ifdef __cplusplus
} // extern "C"
#endif
