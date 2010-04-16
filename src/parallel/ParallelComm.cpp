#include "moab/Interface.hpp"
#include "moab/ParallelComm.hpp"
#include "moab/WriteUtilIface.hpp"
#include "moab/ReadUtilIface.hpp"
#include "SequenceManager.hpp"
#include "EntitySequence.hpp"
#include "TagServer.hpp"
#include "MBTagConventions.hpp"
#include "moab/Skinner.hpp"
#include "MBParallelConventions.h"
#include "moab/Core.hpp"
#include "Error.hpp"
#include "ElementSequence.hpp"
#include "moab/CN.hpp"
#include "moab/RangeMap.hpp"
#include "moab/MeshTopoUtil.hpp"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <functional>
#include <numeric>

#define MIN(a,b) (a < b ? a : b)
#define MAX(a,b) (a > b ? a : b)
const bool debug = false;

#include <math.h>
#include <assert.h>



extern "C" 
{
#include "types.h"
#include "minmax.h"
#include "gs.h"
#include "errmem.h"
#include "sort.h"
#include "tuple_list.h"
}

#ifdef USE_MPI
#include "moab_mpi.h"
#endif
#undef DEBUG_MPE
//#define DEBUG_MPE 1
#ifdef DEBUG_MPE
#include "mpe.h"
int IFACE_START, IFACE_END;
int GHOST_START, GHOST_END;
int SHAREDV_START, SHAREDV_END;
int RESOLVE_START, RESOLVE_END;
int ENTITIES_START, ENTITIES_END;
int RHANDLES_START, RHANDLES_END;
#endif

namespace moab {

const unsigned int ParallelComm::INITIAL_BUFF_SIZE = 1024;

const int MAX_BCAST_SIZE = (1<<28);

#undef DEBUG_COMM
//#define DEBUG_COMM 1
#undef DEBUG_PACKING
#undef DEBUG_MSGS
//#define DEBUG_MSGS 1
#undef DEBUG_BARRIER
#define DEBUG_BARRIER 1
#ifdef DEBUG_MSGS
std::vector<ParallelComm::Buffer*> msgs;
#endif
#ifdef DEBUG_PACKING
unsigned int __PACK_num = 0, __UNPACK_num = 0, __PACK_count = 0, __UNPACK_count = 0;
std::string __PACK_string, __UNPACK_string;

#define PC(n, m) {\
          if (__PACK_num == (unsigned int)n && __PACK_string == m) __PACK_count++;\
          else {\
            if (__PACK_count > 1) std::cerr << " (" << __PACK_count << "x)";\
            __PACK_count = 1; __PACK_string = m; __PACK_num = n;\
            std::cerr << std::endl << "PACK: " << n << m;\
          }}
#define UPC(n, m) {\
          if (__UNPACK_num == (unsigned int)n && __UNPACK_string == m) __UNPACK_count++;\
          else {\
            if (__UNPACK_count > 1) std::cerr << "(" << __UNPACK_count << "x)";\
            __UNPACK_count = 1; __UNPACK_string = m; __UNPACK_num = n;\
            std::cerr << std::endl << "UNPACK: " << n << m;\
          }}
#else
#define PC(n, m)
#define UPC(n, m)
#endif

template <typename T> static inline
void UNPACK( unsigned char*& buff, T* val, size_t count )
{
  memcpy( val, buff, count*sizeof(T) );
  buff += count*sizeof(T);
}

template <typename T> static inline
void PACK( unsigned char*& buff, const T* val, size_t count )
{
  memcpy( buff, val, count*sizeof(T) );
  buff += count*sizeof(T);
}

static inline
void PACK_INTS( unsigned char*& buff, const int* int_val, size_t num )
  { 
    unsigned long _extra;
    ALIGN_BY_SIZE<int>(buff, _extra);
    PACK( buff, int_val, num ); PC(num, " ints"); }

static inline
void PACK_INT( unsigned char*& buff, int int_val )
  { PACK_INTS( buff, &int_val, 1 ); }

static inline
void PACK_DBL( unsigned char*& buff, const double* dbl_val, size_t num )
  { 
    unsigned long extra;
    ALIGN_BUFFER(buff, extra);
    PACK( buff, dbl_val, num ); PC(num, " doubles"); }

static inline
void PACK_EH( unsigned char*& buff, const EntityHandle* eh_val, size_t num )
  { 
    unsigned long extra;
    ALIGN_BUFFER(buff, extra);
    PACK( buff, eh_val, num ); PC(num, " handles"); }

static inline
void PACK_CHAR_64( unsigned char*& buff, const char* str )
{
  strncpy( reinterpret_cast<char*>(buff), str, 64 );
  buff += 64;
  PC(64, " chars");
}

static inline
void PACK_VOID( unsigned char*& buff, const void* val, size_t num )
{
  PACK( buff, reinterpret_cast<const unsigned char*>(val), num );
  PC(num, " void");
}

static inline
void PACK_BYTES( unsigned char*& buff, const void* val, int num )
  { PACK_INT(buff, num); PACK_VOID(buff, val, num); }

static inline
void PACK_RANGE( unsigned char*& buff, const Range& rng )
{
  PACK_INT( buff, rng.psize() );
  Range::const_pair_iterator cit;
  for (cit = rng.const_pair_begin(); cit != rng.const_pair_end(); ++cit) {
    EntityHandle eh[2] = { cit->first, cit->second };
    PACK_EH(buff, eh, 2);
  }
  PC(rng.psize(), "-subranged range");
}

static inline
void UNPACK_INTS( unsigned char*& buff, int* int_val, size_t num )
  { 
    unsigned long _extra;
    ALIGN_BY_SIZE<int>(buff, _extra);
    UNPACK(buff, int_val, num); UPC(num, " ints"); 
}

static inline
void UNPACK_INT( unsigned char*& buff, int& int_val )
  { UNPACK_INTS( buff, &int_val, 1 ); }

static inline
void UNPACK_DBL( unsigned char*& buff, double* dbl_val, size_t num )
  { 
    long _extra;
    ALIGN_BUFFER(buff, _extra);
    UNPACK(buff, dbl_val, num); UPC(num, " doubles"); }

static inline
void UNPACK_EH( unsigned char*& buff, EntityHandle* eh_val, size_t num )
  { 
    long _extra;
    ALIGN_BUFFER(buff, _extra);
    UNPACK(buff, eh_val, num); UPC(num, " handles"); }

static inline
void UNPACK_CHAR_64( unsigned char*& buff, char* char_val )
{
  memcpy( buff, char_val, 64 );
  buff += 64;
  UPC(64, " chars");
}

static inline
void UNPACK_VOID( unsigned char*& buff, void* val, size_t num )
{
  UNPACK(buff, reinterpret_cast<unsigned char*>(val), num);
  UPC(num, " void");
}

static inline
void UNPACK_TYPE( unsigned char*& buff, EntityType& type )
{
  int int_type = MBMAXTYPE;
  UNPACK_INT(buff, int_type);
  type = static_cast<EntityType>(int_type);
  assert(type >= MBVERTEX && type <= MBMAXTYPE);
}

static inline
void UNPACK_RANGE( unsigned char*& buff, Range& rng )
{
  int num_subs;
  EntityHandle eh[2];
  UNPACK_INT( buff, num_subs );
  for (int i = 0; i < num_subs; ++i) {
    UPC(num_subs, "-subranged range"); 
    UNPACK_EH(buff, eh, 2); 
    rng.insert(eh[0], eh[1]);
  }
}    

enum MBMessageTag {MB_MESG_ANY=MPI_ANY_TAG, 
                   MB_MESG_ENTS_ACK,
                   MB_MESG_ENTS_SIZE,
                   MB_MESG_ENTS_LARGE,
                   MB_MESG_REMOTEH_ACK,
                   MB_MESG_REMOTEH_SIZE,
                   MB_MESG_REMOTEH_LARGE,
                   MB_MESG_TAGS_ACK,
                   MB_MESG_TAGS_SIZE,
                   MB_MESG_TAGS_LARGE};
    
static inline size_t RANGE_SIZE(const Range& rng)
  { return 2*sizeof(EntityHandle)*rng.psize()+sizeof(int); }

static inline void PRINT_DEBUG_ISEND(int from, int to, unsigned char *buff,
                                     int tag, int size) 
{
#ifdef DEBUG_COMM
  std::cerr << "Isend, " << from << "->" << to
            << ", buffer ptr = " << (void*)buff << ", tag=" << tag 
            << ", size=" << size << std::endl; std::cerr.flush();
#endif
}

static inline void PRINT_DEBUG_IRECV(int to, int from, unsigned char *buff, int size,
                                     int tag, int incoming) 
{
#ifdef DEBUG_COMM
  std::cerr << "Irecv, " << to << "<-" << from << ", buffer ptr=" << (void*)buff
            << ", size=" << size << ", tag=" << tag
            << (tag < MB_MESG_REMOTEH_ACK ? ", incoming1=" : 
                (tag < MB_MESG_TAGS_ACK ? ", incoming2=" : ", incoming="))
            << incoming << std::endl; std::cerr.flush();
#endif
}

static inline void PRINT_DEBUG_RECD(MPI_Status status) 
{
#ifdef DEBUG_COMM
  int this_count;
  int success = MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &this_count);
  if (MPI_SUCCESS != success) this_count = -1;
  std::cerr << "Received from " << status.MPI_SOURCE
            << ", count = " << this_count 
            << ", tag = " << status.MPI_TAG
            << std::endl; std::cerr.flush();
#endif    
}

static inline void PRINT_DEBUG_WAITANY(std::vector<MPI_Request> &reqs, int tag, int proc) 
{
#ifdef DEBUG_COMM
  std::cerr << "Waitany, p=" << proc
            << (tag < MB_MESG_REMOTEH_ACK ? ", recv_ent_reqs = " : 
                (tag < MB_MESG_TAGS_ACK ? ", recv_remoteh_reqs = " : ", recv_tag_reqs = "));
  for (unsigned int i = 0; i < reqs.size(); i++) std::cerr << " " << reqs[i];
  std::cerr << std::endl; std::cerr.flush();
#endif
}


#define RR(a) if (MB_SUCCESS != result) {\
          dynamic_cast<Core*>(mbImpl)->get_error_handler()->set_last_error(a);\
          return result;}

#define RRA(a) if (MB_SUCCESS != result) {\
      std::string tmp_str; mbImpl->get_last_error(tmp_str);\
      tmp_str.append("\n"); tmp_str.append(a);\
      dynamic_cast<Core*>(mbImpl)->get_error_handler()->set_last_error(tmp_str); \
      return result;}

#define RRAI(i, a) if (MB_SUCCESS != result) {                \
      std::string tmp_str; i->get_last_error(tmp_str);\
      tmp_str.append("\n"); tmp_str.append(a);\
      dynamic_cast<Core*>(i)->get_error_handler()->set_last_error(tmp_str); \
      return result;}

/** Name of tag used to store ParallelComm Index on mesh paritioning sets */
const char* PARTITIONING_PCOMM_TAG_NAME = "__PRTN_PCOMM";
 
/** \brief Tag storing parallel communication objects
 *
 * This tag stores pointers to ParallelComm communication
 * objects; one of these is allocated for each different
 * communicator used to read mesh.  ParallelComm stores
 * partition and interface sets corresponding to its parallel mesh.
 * By default, a parallel read uses the first ParallelComm object
 * on the interface instance; if instantiated with one, ReadParallel
 * adds this object to the interface instance too.
 *
 * Tag type: opaque
 * Tag size: MAX_SHARING_PROCS*sizeof(ParallelComm*)
 */
#define PARALLEL_COMM_TAG_NAME "__PARALLEL_COMM"


ParallelComm::ParallelComm(Interface *impl, MPI_Comm comm, int* id ) 
        : mbImpl(impl), procConfig(comm),
          sharedpTag(0), sharedpsTag(0), 
          sharedhTag(0), sharedhsTag(0), pstatusTag(0), ifaceSetsTag(0),
          partitionTag(0), globalPartCount(-1), partitioningSet(0)
{
  initialize();
  
  if (id)
    *id = pcommID;
}

ParallelComm::ParallelComm(Interface *impl,
                               std::vector<unsigned char> &tmp_buff, 
                               MPI_Comm comm,
                               int* id) 
    : mbImpl(impl), procConfig(comm),
      sharedpTag(0), sharedpsTag(0), 
      sharedhTag(0), sharedhsTag(0), pstatusTag(0), ifaceSetsTag(0),
      partitionTag(0), globalPartCount(-1), partitioningSet(0)
{
  initialize();
  
  if (id)
    *id = pcommID;
}

ParallelComm::~ParallelComm() 
{
  remove_pcomm(this);
  delete_all_buffers();
}

void ParallelComm::initialize() 
{
  tagServer = dynamic_cast<Core*>(mbImpl)->tag_server();
  sequenceManager = dynamic_cast<Core*>(mbImpl)->sequence_manager();
  
    // initialize MPI, if necessary
  int flag = 1;
  int retval = MPI_Initialized(&flag);
  if (MPI_SUCCESS != retval || !flag) {
    int argc = 0;
    char **argv = NULL;
    
      // mpi not initialized yet - initialize here
    retval = MPI_Init(&argc, &argv);
  }

    // reserve space for vectors
  buffProcs.reserve(MAX_SHARING_PROCS);
  localOwnedBuffs.reserve(MAX_SHARING_PROCS);
  remoteOwnedBuffs.reserve(MAX_SHARING_PROCS);

  pcommID = add_pcomm(this);
}

int ParallelComm::add_pcomm(ParallelComm *pc) 
{
    // add this pcomm to instance tag
  std::vector<ParallelComm *> pc_array(MAX_SHARING_PROCS, 
                                         (ParallelComm*)NULL);
  Tag pc_tag = pcomm_tag(mbImpl, true);
  assert(0 != pc_tag);
  
  ErrorCode result = mbImpl->tag_get_data(pc_tag, 0, 0, (void*)&pc_array[0]);
  if (MB_SUCCESS != result && MB_TAG_NOT_FOUND != result) 
    return -1;
  int index = 0;
  while (index < MAX_SHARING_PROCS && pc_array[index]) index++;
  if (index == MAX_SHARING_PROCS) {
    index = -1;
    assert(false);
  }
  else {
    pc_array[index] = pc;
    mbImpl->tag_set_data(pc_tag, 0, 0, (void*)&pc_array[0]);
  }
  return index;
}

void ParallelComm::remove_pcomm(ParallelComm *pc) 
{
    // remove this pcomm from instance tag
  std::vector<ParallelComm *> pc_array(MAX_SHARING_PROCS);
  Tag pc_tag = pcomm_tag(mbImpl, true);
  
  ErrorCode result = mbImpl->tag_get_data(pc_tag, 0, 0, (void*)&pc_array[0]);
  std::vector<ParallelComm*>::iterator pc_it = 
    std::find(pc_array.begin(), pc_array.end(), pc);
  assert(MB_SUCCESS == result && 
         pc_it != pc_array.end());
  *pc_it = NULL;
  mbImpl->tag_set_data(pc_tag, 0, 0, (void*)&pc_array[0]);
}

//! assign a global id space, for largest-dimension or all entities (and
//! in either case for vertices too)
ErrorCode ParallelComm::assign_global_ids(EntityHandle this_set,
                                              const int dimension, 
                                              const int start_id,
                                              const bool largest_dim_only,
                                              const bool parallel) 
{
  Range entities[4];
  int local_num_elements[4];
  ErrorCode result;
  std::vector<unsigned char> pstatus;
  for (int dim = 0; dim <= dimension; dim++) {
    if (dim == 0 || !largest_dim_only || dim == dimension) {
      result = mbImpl->get_entities_by_dimension(this_set, dim, entities[dim]); 
      RRA("Failed to get vertices in assign_global_ids.");
    }

      // need to filter out non-locally-owned entities!!!
    pstatus.resize(entities[dim].size());
    result = mbImpl->tag_get_data(pstatus_tag(), entities[dim], &pstatus[0]);
    RRA("Failed to get pstatus in assign_global_ids.");
    
    Range dum_range;
    Range::iterator rit;
    unsigned int i;
    for (rit = entities[dim].begin(), i = 0; rit != entities[dim].end(); rit++, i++)
      if (pstatus[i] & PSTATUS_NOT_OWNED)
        dum_range.insert(*rit);
    entities[dim] = subtract( entities[dim], dum_range);
    
    local_num_elements[dim] = entities[dim].size();
  }
  
    // communicate numbers
  std::vector<int> num_elements(procConfig.proc_size()*4);
#ifdef USE_MPI
  if (procConfig.proc_size() > 1 && parallel) {
    int retval = MPI_Allgather(local_num_elements, 4, MPI_INT,
                               &num_elements[0], 4, 
                               MPI_INT, procConfig.proc_comm());
    if (0 != retval) return MB_FAILURE;
  }
  else
#endif
    for (int dim = 0; dim < 4; dim++) num_elements[dim] = local_num_elements[dim];
  
    // my entities start at one greater than total_elems[d]
  int total_elems[4] = {start_id, start_id, start_id, start_id};
  
  for (unsigned int proc = 0; proc < procConfig.proc_rank(); proc++) {
    for (int dim = 0; dim < 4; dim++) total_elems[dim] += num_elements[4*proc + dim];
  }
  
    //assign global ids now
  Tag gid_tag;
  int zero = 0;
  result = mbImpl->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int), 
                              MB_TAG_DENSE, MB_TYPE_INTEGER, gid_tag,
                              &zero, true);
  if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) return result;
  
  for (int dim = 0; dim < 4; dim++) {
    if (entities[dim].empty()) continue;
    num_elements.resize(entities[dim].size());
    int i = 0;
    for (Range::iterator rit = entities[dim].begin(); rit != entities[dim].end(); rit++)
      num_elements[i++] = total_elems[dim]++;
    
    result = mbImpl->tag_set_data(gid_tag, entities[dim], &num_elements[0]); 
    RRA("Failed to set global id tag in assign_global_ids.");
  }
  
  return MB_SUCCESS;
}

int ParallelComm::get_buffers(int to_proc, bool *is_new) 
{
  int ind = -1;
  std::vector<unsigned int>::iterator vit = 
    std::find(buffProcs.begin(), buffProcs.end(), to_proc);
  if (vit == buffProcs.end()) {
    ind = buffProcs.size();
    buffProcs.push_back((unsigned int)to_proc);
    localOwnedBuffs.push_back(new Buffer(INITIAL_BUFF_SIZE));
    remoteOwnedBuffs.push_back(new Buffer(INITIAL_BUFF_SIZE));
    if (is_new) *is_new = true;
  }
  else {
    ind = vit - buffProcs.begin();
    if (is_new) *is_new = false;
  }
  assert(ind < MAX_SHARING_PROCS);
  return ind;
}

ErrorCode ParallelComm::broadcast_entities( const int from_proc,
                                                Range &entities,
                                                const bool adjacencies,
                                                const bool tags) 
{
#ifndef USE_MPI
  return MB_FAILURE;
#else
  
  ErrorCode result = MB_SUCCESS;
  int success;
  int buff_size;

  Buffer buff(INITIAL_BUFF_SIZE);
  buff.reset_ptr(sizeof(int));
  if ((int)procConfig.proc_rank() == from_proc) {
    result = add_verts(entities);
    RRA("Failed to add adj vertices.");

    buff.reset_ptr(sizeof(int));
    result = pack_buffer( entities, adjacencies, tags, 
                          false, -1, &buff); 
    RRA("Failed to compute buffer size in broadcast_entities.");
    buff.set_stored_size();
    buff_size = buff.buff_ptr - buff.mem_ptr;
  }

  success = MPI_Bcast( &buff_size, 1, MPI_INT, from_proc, procConfig.proc_comm() );
  if (MPI_SUCCESS != success) {
    result = MB_FAILURE;
    RRA("MPI_Bcast of buffer size failed.");
  }
  
  if (!buff_size) // no data
    return MB_SUCCESS;

  if ((int)procConfig.proc_rank() != from_proc) 
    buff.reserve(buff_size);

  size_t offset = 0;
  while (buff_size) {
    int size = std::min( buff_size, MAX_BCAST_SIZE );
    success = MPI_Bcast(buff.mem_ptr+offset, size, MPI_UNSIGNED_CHAR, from_proc, procConfig.proc_comm() );
    if (MPI_SUCCESS != success) {
      result = MB_FAILURE;
      RRA("MPI_Bcast of buffer failed.");
    }
    
    offset += size;
    buff_size -= size;
  }

  if ((int)procConfig.proc_rank() != from_proc) {
    std::vector<std::vector<EntityHandle> > dum1a, dum1b;
    std::vector<std::vector<int> > dum1p;
    std::vector<EntityHandle> dum2;
    std::vector<unsigned int> dum3;
    buff.reset_ptr(sizeof(int));
    result = unpack_buffer(&buff, false, from_proc, -1, 
                           dum1a, dum1b, dum1p, dum2, dum2, dum3, entities);
    RRA("Failed to unpack buffer in broadcast_entities.");
  }

  return MB_SUCCESS;
#endif
}

ErrorCode ParallelComm::pack_buffer(Range &orig_ents, 
                                        const bool adjacencies,
                                        const bool tags,
                                        const bool store_remote_handles,
                                        const int to_proc,
                                        Buffer *buff) 
{
    // pack the buffer with the entity ranges, adjacencies, and tags sections
    // 
    // Note: new entities used in subsequent connectivity lists, sets, or tags, 
    //   are referred to as (MBMAXTYPE + index), where index is into vector 
    //   of new entities, 0-based
  ErrorCode result;

  Range set_range;
  std::vector<Range> set_ranges;
  std::vector<Tag> all_tags;
  std::vector<Range> tag_ranges;
  std::vector<int> set_sizes;
  std::vector<unsigned int> options_vec;

  Range::const_iterator rit;

    // entities
  result = pack_entities(orig_ents, buff,
                         store_remote_handles, to_proc, false); 
  RRA("Packing entities failed.");
  
    // sets
  result = pack_sets(orig_ents, buff,
                     store_remote_handles, to_proc); 
  RRA("Packing sets (count) failed.");

    // tags
  Range final_ents;
  if (tags) {
    result = get_tag_send_list(orig_ents, all_tags, tag_ranges );
    RRA("Failed to get tagged entities.");
    result = pack_tags(orig_ents, all_tags, all_tags, tag_ranges, 
                       buff, store_remote_handles, to_proc);
    RRA("Packing tags (count) failed.");
  }

  return result;
}
 
ErrorCode ParallelComm::unpack_buffer(Buffer *buff,
                                      const bool store_remote_handles,
                                          const int from_proc,
                                          const int ind,
                                          std::vector<std::vector<EntityHandle> > &L1hloc,
                                          std::vector<std::vector<EntityHandle> > &L1hrem,
                                          std::vector<std::vector<int> > &L1p,
                                          std::vector<EntityHandle> &L2hloc, 
                                          std::vector<EntityHandle> &L2hrem,
                                          std::vector<unsigned int> &L2p,
                                          Range &new_ents) 
{
#ifdef DEBUG_PACKING
  std::cerr << "Unpack buffer: starting at " << buff->get_size() << std::endl;
#endif  
    ErrorCode result;
    result = unpack_entities(buff, store_remote_handles,
                             ind, false, L1hloc, L1hrem, L1p, L2hloc, L2hrem, L2p, new_ents);
  RRA("Unpacking entities failed.");

  result = unpack_sets(buff, new_ents, store_remote_handles, from_proc);
  RRA("Unpacking sets failed.");

  result = unpack_tags(buff, new_ents, store_remote_handles, from_proc);
  RRA("Unpacking tags failed.");

#ifdef DEBUG_PACKING
  std::cerr << std::endl;
#endif
  
  return MB_SUCCESS;
}

int ParallelComm::num_subranges(const Range &this_range)
{
    // ok, have all the ranges we'll pack; count the subranges
  int num_sub_ranges = 0;
  for (Range::const_pair_iterator pit = this_range.const_pair_begin(); 
       pit != this_range.const_pair_end(); pit++)
    num_sub_ranges++;

  return num_sub_ranges;
}

int ParallelComm::estimate_ents_buffer_size(Range &entities,
                                            const bool store_remote_handles,
                                            unsigned int &num_aligns) 
{
  int buff_size = 0;
  num_aligns = 0;
  std::vector<EntityHandle> dum_connect_vec;
  const EntityHandle *connect;
  int num_connect;

  int num_verts = entities.num_of_type(MBVERTEX);
    // # verts + coords + handles
  buff_size += 2*sizeof(int) + 3*sizeof(double)*num_verts;
  num_aligns++;
  if (store_remote_handles) {
    buff_size += sizeof(EntityHandle)*num_verts;
    num_aligns++;
  }

    // do a rough count by looking at first entity of each type
  for (EntityType t = MBEDGE; t < MBENTITYSET; t++) {
    const Range::iterator rit = entities.lower_bound(t);
    if (TYPE_FROM_HANDLE(*rit) != t) continue;
    
    ErrorCode result = mbImpl->get_connectivity(*rit, connect, num_connect, 
                                                false, &dum_connect_vec);
    RRA("Failed to get connectivity to estimate buffer size.");

      // number, type, nodes per entity
    buff_size += 3*sizeof(int);
    int num_ents = entities.num_of_type(t);
      // connectivity, handle for each ent
    buff_size += (num_connect+1)*sizeof(EntityHandle)*num_ents;
    num_aligns++;
  }

      // extra entity type at end, passed as int
  buff_size += sizeof(int);

  return buff_size;
}

int ParallelComm::estimate_sets_buffer_size(Range &entities,
                                            const bool store_remote_handles,
                                            unsigned int &num_aligns) 
{
  num_aligns = 0;
  
    // number of sets
  int buff_size = sizeof(int);
  
    // do a rough count by looking at first entity of each type
  Range::iterator rit = entities.lower_bound(MBENTITYSET);
  ErrorCode result;
  
  for (; rit != entities.end(); rit++) {
    unsigned int options;
    result = mbImpl->get_meshset_options(*rit, options);
    RRA("Failed to get meshset options.");

    buff_size += sizeof(int);
    
    Range set_range;
    if (options & MESHSET_SET) {
        // range-based set; count the subranges
      result = mbImpl->get_entities_by_handle(*rit, set_range);
      RRA("Failed to get set entities.");

        // set range
      int rsize = RANGE_SIZE(set_range);
      buff_size += rsize;
      num_aligns += rsize;
    }
    else if (options & MESHSET_ORDERED) {
        // just get the number of entities in the set
      int num_ents;
      result = mbImpl->get_number_entities_by_handle(*rit, num_ents);
      RRA("Failed to get number entities in ordered set.");

        // set vec
      buff_size += sizeof(EntityHandle) * num_ents + sizeof(int);
      num_aligns++;
    }

      // get numbers of parents/children
    int num_par, num_ch;
    result = mbImpl->num_child_meshsets(*rit, &num_ch);
    RRA("Failed to get num children.");

    result = mbImpl->num_parent_meshsets(*rit, &num_par);
    RRA("Failed to get num parents.");

    buff_size += (num_ch + num_par) * sizeof(EntityHandle) + 2*sizeof(int);
    num_aligns += 2;
  }

  return buff_size;
}

ErrorCode ParallelComm::pack_entities(Range &entities,
                                          Buffer *buff,
                                          const bool store_remote_handles,
                                          const int to_proc,
                                          const bool is_iface,
                                          std::vector<std::set<unsigned int> > *entprocs,
                                          Range *allsent) 
{
    // packed information:
    // 1. # entities = E
    // 2. for e in E
    //   a. # procs sharing e, incl. sender and receiver = P
    //   b. for p in P (procs sharing e)
    //   c. for p in P (handle for e on p) (Note1)
    // 3. vertex/entity info

    // get an estimate of the buffer size & pre-allocate buffer size
  unsigned int num_aligns;
  unsigned int buff_size = estimate_ents_buffer_size(entities, 
                                                     store_remote_handles, 
                                                     num_aligns);
  buff->check_space(buff_size, num_aligns);
  
  WriteUtilIface *wu;
  ErrorCode result = mbImpl->query_interface(std::string("WriteUtilIface"), 
                                               reinterpret_cast<void**>(&wu));
  RRA("Couldn't get WriteUtilIface.");

  unsigned int num_ents;

    // first pack procs/handles sharing this ent, not including this dest but including
    // others (with zero handles)
  if (store_remote_handles) {

      // buff space is at least proc+handle for each entity; use avg of 4 other procs
      // to estimate buff size, but check later
    buff->check_space(sizeof(int) + (5*sizeof(int) + sizeof(EntityHandle))*entities.size(),
                      entities.size());

      // 1. # entities = E
    PACK_INT(buff->buff_ptr, entities.size());
  
    Range::iterator rit;
  
      // pre-fetch sharedp and pstatus
    std::vector<int> sharedp_vals(entities.size());
    result = mbImpl->tag_get_data(sharedp_tag(), entities, &sharedp_vals[0]);
    RRA("Failed to get sharedp_tag.");
    std::vector<char> pstatus_vals(entities.size());
    result = mbImpl->tag_get_data(pstatus_tag(), entities, &pstatus_vals[0]);
    RRA("Failed to get sharedp_tag.");
  
    unsigned int i;
    std::vector<std::set<unsigned int> >::iterator vit;
    int tmp_procs[MAX_SHARING_PROCS];
    EntityHandle tmp_handles[MAX_SHARING_PROCS];
    std::set<unsigned int> dumprocs;

      // 2. for e in E
    for (rit = entities.begin(), i = 0; 
         rit != entities.end(); rit++, i++) {
      result = build_sharedhps_list(*rit, pstatus_vals[i], sharedp_vals[i],
                                    (entprocs ? (*entprocs)[allsent->index(*rit)] : dumprocs),
                                    num_ents, tmp_procs, tmp_handles);
      RRA("Failed to build sharedhps.");

        // now pack them
      buff->check_space((num_ents+1)*sizeof(int) + 
                       num_ents*sizeof(EntityHandle));
      PACK_INT(buff->buff_ptr, num_ents);
      PACK_INTS(buff->buff_ptr, tmp_procs, num_ents);
      PACK_EH(buff->buff_ptr, tmp_handles, num_ents);

#ifndef NDEBUG
        // check for duplicates in proc list
      unsigned int dp = 0;
      for (; dp < MAX_SHARING_PROCS && -1 != tmp_procs[dp]; dp++)
        dumprocs.insert(tmp_procs[dp]);
      assert(dumprocs.size() == dp);
      dumprocs.clear();
#endif      
    }
  }
  
    // pack vertices
  Range these_ents = entities.subset_by_type(MBVERTEX);
  num_ents = these_ents.size();

  if (num_ents) {
    buff_size = 2*sizeof(int) + 3*num_ents*sizeof(double);
    buff->check_space(buff_size);

    // type, # ents
    PACK_INT(buff->buff_ptr, ((int) MBVERTEX));
    PACK_INT(buff->buff_ptr, ((int) num_ents));

    unsigned long extra;
    ALIGN_BUFFER(buff->buff_ptr, extra);
    result = mbImpl->get_coords(these_ents, (double*)buff->buff_ptr);
    RRA("Couldn't get vertex coordinates.");
    PC(3*num_ents, " doubles");

    buff->buff_ptr += 3 * num_ents * sizeof(double);

#ifdef DEBUG_PACKING
  std::cerr << "Packed " << these_ents.size() << " ents of type " 
            << CN::EntityTypeName(TYPE_FROM_HANDLE(*these_ents.begin())) << std::endl;
#endif      
  }

    // now entities; go through range, packing by type and equal # verts per element
  Range::iterator start_rit = entities.find(*these_ents.rbegin());
  start_rit++;
  int last_nodes = -1;
  EntityType last_type = MBMAXTYPE;
  these_ents.clear();
  Range::iterator end_rit = start_rit;
  EntitySequence *seq;
  ElementSequence *eseq;
  
  while (start_rit != entities.end() || !these_ents.empty()) {
      // cases:
      // A: !end, last_type == MBMAXTYPE, seq: save contig sequence in these_ents
      // B: !end, last type & nodes same, seq: save contig sequence in these_ents
      // C: !end, last type & nodes different: pack these_ents, then save contig sequence in these_ents
      // D: end: pack these_ents

      // find the sequence holding current start entity, if we're not at end
    eseq = NULL;
    if (start_rit != entities.end()) {
      result = sequenceManager->find(*start_rit, seq);
      RRA("Couldn't find entity sequence.");
      if (NULL == seq) return MB_FAILURE;
      eseq = dynamic_cast<ElementSequence*>(seq);
    }

      // pack the last batch if at end or next one is different
    if (!these_ents.empty() &&
        (!eseq || eseq->type() != last_type ||
         last_nodes != (int) eseq->nodes_per_element())) {
      result = pack_entity_seq(last_nodes, store_remote_handles,
                               to_proc, these_ents, entities, buff);
      RRA("Failed to pack entities from a sequence.");
      these_ents.clear();
    }

    if (eseq) {
        // continuation of current range, just save these entities
        // get position in entities list one past end of this sequence
      end_rit = entities.lower_bound(start_rit, entities.end(), eseq->end_handle()+1);

        // put these entities in the range
      std::copy(start_rit, end_rit, range_inserter(these_ents));

      last_type = eseq->type();
      last_nodes = eseq->nodes_per_element();
    }
    else if (start_rit != entities.end() &&
             TYPE_FROM_HANDLE(*start_rit) == MBENTITYSET)
      break;

    start_rit = end_rit;
  }

    // pack MBMAXTYPE to indicate end of ranges
  buff->check_space(sizeof(int));
  PACK_INT(buff->buff_ptr, ((int)MBMAXTYPE));

  buff->set_stored_size();

  if (22 == procConfig.proc_rank() || 23 == procConfig.proc_rank())
    print_buffer(buff->mem_ptr, MB_MESG_ENTS_LARGE, to_proc, true);

  return MB_SUCCESS;
}

ErrorCode ParallelComm::build_sharedhps_list(const EntityHandle entity,
                                                 const unsigned char pstatus,
                                                 const int sharedp, 
                                                 const std::set<unsigned int> &entprocs,
                                                 unsigned int &num_ents,
                                                 int *tmp_procs,
                                                 EntityHandle *tmp_handles)
{
  num_ents = 0;
  unsigned char pstat;
  ErrorCode result = get_sharing_data(entity, tmp_procs, tmp_handles,
                                        pstat, num_ents);
  RRA("Failed in get_sharing_data.");
  assert(pstat == pstatus);
  
    // build shared proc/handle lists
    // start with multi-shared, since if it is the owner will be first
  if (pstatus & PSTATUS_MULTISHARED) {
  }
  else if (pstatus & PSTATUS_NOT_OWNED) {
      // if not multishared and not owned, other sharing proc is owner, put that
      // one first
    assert("If not owned, I should be shared too" &&
           pstatus & PSTATUS_SHARED &&
           num_ents == 1);
    tmp_procs[1] = procConfig.proc_rank();
    tmp_handles[1] = entity;
    num_ents = 2;
  }
  else if (pstatus & PSTATUS_SHARED) {
      // if not multishared and owned, I'm owner
    assert("shared and owned, should be only 1 sharing proc" &&
           1 == num_ents);
    tmp_procs[1] = tmp_procs[0];
    tmp_procs[0] = procConfig.proc_rank();
    tmp_handles[1] = tmp_handles[0];
    tmp_handles[0] = entity;
    num_ents = 2;
  }
  else {
      // not shared yet, just add owner (me)
    tmp_procs[0] = procConfig.proc_rank();
    tmp_handles[0] = entity;
    num_ents = 1;
  }

#ifndef NDEBUG
  int tmp_ps = num_ents;
#endif
  
    // now add others, with zero handle for now
  for (std::set<unsigned int>::iterator sit = entprocs.begin();
       sit != entprocs.end(); sit++) {
    assert("these procs shouldn't already be in the shared list" &&
           std::find(tmp_procs, tmp_procs+tmp_ps, *sit) == tmp_procs+tmp_ps);
    tmp_procs[num_ents] = *sit;
    tmp_handles[num_ents] = 0;
    num_ents++;
  }

    // put -1 after procs and 0 after handles
  if (MAX_SHARING_PROCS > num_ents) {
    tmp_procs[num_ents] = -1;
    tmp_handles[num_ents] = 0;
  }
  
  return MB_SUCCESS;
}

ErrorCode ParallelComm::pack_entity_seq(const int nodes_per_entity,
                                            const bool store_remote_handles,
                                            const int to_proc,
                                            Range &these_ents,
                                            Range &entities,
                                            Buffer *buff) 
{
  int tmp_space = 3*sizeof(int) + nodes_per_entity*these_ents.size()*sizeof(EntityHandle);
  buff->check_space(tmp_space);
  
    // pack the entity type
  PACK_INT(buff->buff_ptr, ((int)TYPE_FROM_HANDLE(*these_ents.begin())));

    // pack # ents
  PACK_INT(buff->buff_ptr, these_ents.size());
      
    // pack the nodes per entity
  PACK_INT(buff->buff_ptr, nodes_per_entity);
      
    // pack the connectivity
  const EntityHandle *connect;
  int num_connect;
  unsigned long extra;
  ALIGN_BUFFER(buff->buff_ptr, extra);
  std::vector<EntityHandle> dum_connect;
  EntityHandle *start_vec = (EntityHandle*)buff->buff_ptr;
  ErrorCode result = MB_SUCCESS;
  for (Range::const_iterator rit = these_ents.begin(); rit != these_ents.end(); rit++) {
    result = mbImpl->get_connectivity(*rit, connect, num_connect, false,
                                      &dum_connect);
    RRA("Failed to get connectivity.");
    assert(num_connect == nodes_per_entity);
    PACK_EH(buff->buff_ptr, connect, num_connect);
  }

    // substitute destination handles
  result = get_remote_handles(store_remote_handles, start_vec, start_vec,
                              nodes_per_entity*these_ents.size(), to_proc,
                              entities);
  RRA("Trouble getting remote handles when packing entities.");

#ifdef DEBUG_PACKING
  std::cerr << std::endl << "Packed " << these_ents.size() << " ents of type " 
            << CN::EntityTypeName(TYPE_FROM_HANDLE(*these_ents.begin())) << std::endl;
#endif      

  return result;
}


ErrorCode ParallelComm::get_remote_handles(const bool store_remote_handles,
                                               EntityHandle *from_vec, 
                                               EntityHandle *to_vec_tmp,
                                               int num_ents, int to_proc,
                                               const Range &new_ents) 
{
    // NOTE: THIS IMPLEMENTATION IS JUST LIKE THE RANGE-BASED VERSION, NO REUSE
    // AT THIS TIME, SO IF YOU FIX A BUG IN THIS VERSION, IT MAY BE IN THE
    // OTHER VERSION TOO!!!
  if (0 == num_ents) return MB_SUCCESS;
  
    // use a local destination ptr in case we're doing an in-place copy
  std::vector<EntityHandle> tmp_vector;
  EntityHandle *to_vec = to_vec_tmp;
  if (to_vec == from_vec) {
    tmp_vector.resize(num_ents);
    to_vec = &tmp_vector[0];
  }

  if (!store_remote_handles) {
    int err;
      // in this case, substitute position in new_ents list
    for (int i = 0; i < num_ents; i++) {
      int ind = new_ents.index(from_vec[i]);
      to_vec[i] = CREATE_HANDLE(MBMAXTYPE, ind, err);
      assert(to_vec[i] != 0 && !err && -1 != ind);
    }
  }
  else {
    Tag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
    ErrorCode result = get_shared_proc_tags(sharedp_tag, sharedps_tag, 
                                              sharedh_tag, sharedhs_tag, pstatus_tag);
  
      // get single-proc destination handles and shared procs
    std::vector<int> sharing_procs(num_ents);
    result = mbImpl->tag_get_data(sharedh_tag, from_vec, num_ents,
                                  to_vec);
    RRA("Failed to get shared handle tag for remote_handles.");
    result = mbImpl->tag_get_data(sharedp_tag, from_vec, num_ents, &sharing_procs[0]);
    RRA("Failed to get sharing proc tag in remote_handles.");
    for (int j = 0; j < num_ents; j++) {
      if (to_vec[j] && sharing_procs[j] != to_proc)
        to_vec[j] = 0;
    }
    
    EntityHandle tmp_handles[MAX_SHARING_PROCS];
    int tmp_procs[MAX_SHARING_PROCS];
    int i;
      // go through results, and for 0-valued ones, look for multiple shared proc
    EntityHandle *tmp_eh;
    for (tmp_eh = to_vec, i = 0; i < num_ents; i++) {
      if (!to_vec[i]) {
        result = mbImpl->tag_get_data(sharedps_tag, from_vec+i, 1, tmp_procs);
        if (MB_SUCCESS == result) {
          for (int j = 0; j < MAX_SHARING_PROCS; j++) {
            if (-1 == tmp_procs[j]) break;
            else if (tmp_procs[j] == to_proc) {
              result = mbImpl->tag_get_data(sharedhs_tag, from_vec+i, 1, tmp_handles);
              RRA("Trouble getting sharedhs tag.");
              to_vec[i] = tmp_handles[j];
              assert(to_vec[i]);
              break;
            }
          }
        }
        if (!to_vec[i]) {
          int j = new_ents.index(from_vec[i]);
          if (-1 == j) {
            result = MB_FAILURE;
            std::cout << "Failed to find new entity in send list, proc " 
                      << procConfig.proc_rank() << std::endl;
            for (int j = 0; j <= num_ents; j++) 
              std::cout << j << ": " << from_vec[j] << " " << to_vec[j] 
                        << std::endl;
            RRA("Failed to find new entity in send list.");
          }
          int err;
          to_vec[i] = CREATE_HANDLE(MBMAXTYPE, j, err);
          if (err) {
            result = MB_FAILURE;
            RRA("Failed to create handle in remote_handles.");
          }
        }
      }
    }
  }
  
    // memcpy over results if from_vec and to_vec are the same
  if (to_vec_tmp == from_vec) 
    memcpy(from_vec, to_vec, num_ents * sizeof(EntityHandle));
  
  return MB_SUCCESS;
}

ErrorCode ParallelComm::get_remote_handles(const bool store_remote_handles,
                                               const Range &from_range, 
                                               EntityHandle *to_vec,
                                               int to_proc,
                                               const Range &new_ents) 
{
    // NOTE: THIS IMPLEMENTATION IS JUST LIKE THE VECTOR-BASED VERSION, NO REUSE
    // AT THIS TIME, SO IF YOU FIX A BUG IN THIS VERSION, IT MAY BE IN THE
    // OTHER VERSION TOO!!!
  if (from_range.empty()) return MB_SUCCESS;
  
  if (!store_remote_handles) {
    int err;
      // in this case, substitute position in new_ents list
    Range::iterator rit;
    unsigned int i;
    for (rit = from_range.begin(), i = 0; rit != from_range.end(); rit++, i++) {
      int ind = new_ents.index(*rit);
      to_vec[i] = CREATE_HANDLE(MBMAXTYPE, ind, err);
      assert(to_vec[i] != 0 && !err && -1 != ind);
    }
  }
  else {
    Tag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
    ErrorCode result = get_shared_proc_tags(sharedp_tag, sharedps_tag, 
                                              sharedh_tag, sharedhs_tag, pstatus_tag);
  
      // get single-proc destination handles and shared procs
    std::vector<int> sharing_procs(from_range.size());
    result = mbImpl->tag_get_data(sharedh_tag, from_range, to_vec);
    RRA("Failed to get shared handle tag for remote_handles.");
    result = mbImpl->tag_get_data(sharedp_tag, from_range, &sharing_procs[0]);
    RRA("Failed to get sharing proc tag in remote_handles.");
    for (unsigned int j = 0; j < from_range.size(); j++) {
      if (to_vec[j] && sharing_procs[j] != to_proc)
        to_vec[j] = 0;
    }
    
    EntityHandle tmp_handles[MAX_SHARING_PROCS];
    int tmp_procs[MAX_SHARING_PROCS];
      // go through results, and for 0-valued ones, look for multiple shared proc
    Range::iterator rit;
    unsigned int i;
    for (rit = from_range.begin(), i = 0; rit != from_range.end(); rit++, i++) {
      if (!to_vec[i]) {
        result = mbImpl->tag_get_data(sharedhs_tag, &(*rit), 1, tmp_handles);
        if (MB_SUCCESS == result) {
          result = mbImpl->tag_get_data(sharedps_tag, &(*rit), 1, tmp_procs);
          RRA("Trouble getting sharedps tag.");
          for (int j = 0; j < MAX_SHARING_PROCS; j++)
            if (tmp_procs[j] == to_proc) {
              to_vec[i] = tmp_handles[j];
              break;
            }
        }
      
        if (!to_vec[i]) {
          int j = new_ents.index(*rit);
          if (-1 == j) {
            result = MB_FAILURE;
            RRA("Failed to find new entity in send list.");
          }
          int err;
          to_vec[i] = CREATE_HANDLE(MBMAXTYPE, j, err);
          if (err) {
            result = MB_FAILURE;
            RRA("Failed to create handle in remote_handles.");
          }
        }
      }
    }
  }
  
  return MB_SUCCESS;
}

ErrorCode ParallelComm::get_remote_handles(const bool store_remote_handles,
                                               const Range &from_range, 
                                               Range &to_range,
                                               int to_proc,
                                               const Range &new_ents) 
{
  std::vector<EntityHandle> to_vector(from_range.size());

  ErrorCode result =
    get_remote_handles(store_remote_handles, from_range, &to_vector[0],
                       to_proc, new_ents);
  RRA("Trouble getting remote handles.");
  std::copy(to_vector.begin(), to_vector.end(), range_inserter(to_range));
  return result;
}

ErrorCode ParallelComm::unpack_entities(Buffer *buff,
                                        const bool store_remote_handles,
                                        const int from_ind,
                                        const bool is_iface,
                                        std::vector<std::vector<EntityHandle> > &L1hloc,
                                        std::vector<std::vector<EntityHandle> > &L1hrem,
                                        std::vector<std::vector<int> > &L1p,
                                        std::vector<EntityHandle> &L2hloc, 
                                        std::vector<EntityHandle> &L2hrem,
                                        std::vector<unsigned int> &L2p,
                                        Range &new_ents) 
{
    // general algorithm:
    // - unpack # entities
    // - save start of remote handle info, then scan forward to entity definition data
    // - for all vertices or entities w/ same # verts:
    //   . get entity type, num ents, and (if !vert) # verts 
    //   . for each ent:
    //      o get # procs/handles in remote handle info
    //      o if # procs/handles > 2, check for already-created entity:
    //        x get index of owner proc (1st in proc list), resize L1 list if nec
    //        x look for already-arrived entity in L2 by owner handle
    //      o if no existing entity:
    //        x if iface, look for existing entity with same connect & type
    //        x if none found, create vertex or element
    //        x if !iface & multi-shared, save on L2
    //        x if !iface, put new entity on new_ents list
    //      o update proc/handle, pstatus tags, adjusting to put owner first if iface
    //      o if !iface, save new handle on L1 for all sharing procs

    // lists of handles/procs to return to sending/other procs
    // L1hloc[p], L1hrem[p]: handle pairs [h, h'], where h is the local proc handle
    //         and h' is either the remote proc handle (if that is known) or
    //         the owner proc handle (otherwise);
    // L1p[p]: indicates whether h is remote handle (= -1) or owner (rank of owner)
    // L2hloc, L2hrem: local/remote handles for entities shared by > 2 procs;
    //         remote handles are on owning proc
    // L2p: owning procs for handles in L2hrem

  ErrorCode result;
  bool done = false;
  ReadUtilIface *ru = NULL;

  result = mbImpl->query_interface(std::string("ReadUtilIface"), 
                                   reinterpret_cast<void**>(&ru));
  RRA("Failed to get ReadUtilIface.");

#ifdef DEBUG_PACKING
  std::cerr << "Unpack entities: starting at " << buff->get_size() << std::endl;
#endif  

    // procs the sending proc is telling me I'll be receiving from
  std::set<unsigned int> comm_procs;

    // 1. # entities = E
  int num_ents;
  unsigned char *buff_save = buff->buff_ptr;
  unsigned char *buff_orig = buff->buff_ptr;
  int i, j;
  unsigned long extra;

  if (store_remote_handles) {
    UNPACK_INT(buff->buff_ptr, num_ents);

      // save place where remote handle info starts, then scan forward to ents
    for (i = 0; i < num_ents; i++) {
      UNPACK_INT(buff->buff_ptr, j);
      if (j < 0) {
        std::cout << "Should be non-negative # proc/handles.";
        return MB_FAILURE;
      }
      
//      buff_ptr += j * (sizeof(int)+sizeof(EntityHandle));
      buff->buff_ptr += j * sizeof(int);
      ALIGN_BUFFER(buff->buff_ptr, extra);
      buff->buff_ptr += j*sizeof(EntityHandle);
    }
  }

  std::vector<EntityHandle> msg_ents;
  
  while (!done) {
    EntityType this_type = MBMAXTYPE;
    UNPACK_TYPE(buff->buff_ptr, this_type);
    assert(this_type != MBENTITYSET);

      // MBMAXTYPE signifies end of entities data
    if (MBMAXTYPE == this_type) break;

    assert(!is_iface || this_type != MBVERTEX);
    
      // get the number of ents
    int num_ents2, verts_per_entity;
    UNPACK_INT(buff->buff_ptr, num_ents2);

      // unpack the nodes per entity
    if (MBVERTEX != this_type && num_ents2) {
      UNPACK_INT(buff->buff_ptr, verts_per_entity);
    }
      
    std::vector<int> ps(MAX_SHARING_PROCS, -1);
    std::vector<EntityHandle> hs(MAX_SHARING_PROCS, 0);
    for (int e = 0; e < num_ents2; e++) {
        // check for existing entity, otherwise make new one
      EntityHandle new_h = 0;

      EntityHandle *connect;
      double *coords;
      int num_ps = -1;

        //=======================================
        // unpack all the data at once, to make sure the buffer pointers
        // are tracked correctly
        //=======================================
      if (store_remote_handles) {
          // pointers to other procs/handles
        UNPACK_INT(buff_save, num_ps);
        if (0 >= num_ps) {
          std::cout << "Shouldn't ever be fewer than 1 procs here." << std::endl;
          return MB_FAILURE;
        }
        
        UNPACK_INTS(buff_save, &ps[0], num_ps);
        UNPACK_EH(buff_save, &hs[0], num_ps);
      }

      ALIGN_BUFFER(buff->buff_ptr, extra);
      if (MBVERTEX == this_type) {
        coords = (double*) buff->buff_ptr;
        buff->buff_ptr += 3*sizeof(double);
      }
      else {
        connect = (EntityHandle*) buff->buff_ptr;
        buff->buff_ptr += verts_per_entity * sizeof(EntityHandle);

          // update connectivity to local handles
        result = get_local_handles(connect, verts_per_entity, msg_ents);
        RRA("Couldn't get local handles.");
      }

        //=======================================
        // now, process that data; begin by finding an identical 
        // entity, if there is one
        //=======================================
      if (store_remote_handles) {
        result = find_existing_entity(is_iface, ps[0], hs[0], num_ps, 
                                      connect, verts_per_entity,
                                      this_type,
                                      L2hloc, L2hrem, L2p,
                                      new_h);
        RRA("Trouble getting existing entity.");
      }

        //=======================================
        // if we didn't find one, we'll have to create one
        //=======================================
      bool created_here = false;
      if (!new_h && !is_iface) {
        
        if (MBVERTEX == this_type) {
            // create a vertex
          result = mbImpl->create_vertex(coords, new_h);
          RRA("Couldn't make new vertex.");
        }
        else {
            // create the element
          result = mbImpl->create_element(this_type, connect, verts_per_entity, new_h);
          RRA("Couldn't make new vertex.");

            // update adjacencies
          result = ru->update_adjacencies(new_h, 1, 
                                          verts_per_entity, connect);
          RRA("Failed to update adjacencies.");
        }

          // should have a new handle now
        assert(new_h);
        
          // if a new multi-shared entity, save owner for subsequent lookup in L2 lists
        if (store_remote_handles && !is_iface && num_ps > 2) {
          L2hrem.push_back(hs[0]);
          L2hloc.push_back(new_h);
          L2p.push_back(ps[0]);
        }

        created_here = true;
      }

        //=======================================
        // take care of sharing data
        //=======================================

        // need to save entities found in order, for interpretation of
        // later parts of this message
      if (!is_iface) {
        assert(new_h);
        msg_ents.push_back(new_h);
      }

      if (created_here) new_ents.insert(new_h);

      if (new_h && store_remote_handles) {
        
          // update sharing data and pstatus, adjusting order if iface
        result = update_remote_data(new_h, &ps[0], &hs[0], num_ps, 
                                    (is_iface ? PSTATUS_INTERFACE :
                                     (created_here ? (PSTATUS_GHOST | PSTATUS_NOT_OWNED) : 0)));
        RRA("");

          // need to send this new handle to all sharing procs
        if (!is_iface) {
          for (j = 0; j < num_ps; j++) {
            if (ps[j] == (int)procConfig.proc_rank()) continue;
            int idx = get_buffers(ps[j]);
            if (idx == (int)L1hloc.size()) {
              L1hloc.resize(idx+1);
              L1hrem.resize(idx+1);
              L1p.resize(idx+1);
            }
            
              // don't bother adding if it's already in the list
            std::vector<EntityHandle>::iterator vit = 
                std::find(L1hloc[idx].begin(), L1hloc[idx].end(), new_h);
            if (vit != L1hloc[idx].end()) {
                // if it's in the list but remote handle isn't known but we know
                // it, replace in the list
              if (L1p[idx][vit-L1hloc[idx].begin()] != -1 && hs[j]) {
                L1hrem[idx][vit-L1hloc[idx].begin()] = hs[j];
                L1p[idx][vit-L1hloc[idx].begin()] = -1;
              }
              else continue;
            }
            else {
              if (!hs[j]) {
                assert(-1 != ps[0] && num_ps > 2);
                L1p[idx].push_back(ps[0]);
                L1hrem[idx].push_back(hs[0]);
              }
              else {
#ifndef NDEBUG
//                assert("either this remote handle isn't in the remote list, or it's for another proc" &&
//                       (std::find(L1hrem[idx].begin(), L1hrem[idx].end(), hs[j]) == 
//                        L1hrem[idx].end() ||
//                        L1p[idx][std::find(L1hrem[idx].begin(), L1hrem[idx].end(), hs[j]) - 
//                                 L1hrem[idx].begin()] != -1));
                if (!is_iface && 22 == procConfig.proc_rank() && 
                    "either this remote handle isn't in the remote list, or it's for another proc" &&
                    (std::find(L1hrem[idx].begin(), L1hrem[idx].end(), hs[j]) != 
                     L1hrem[idx].end() &&
                     L1p[idx][std::find(L1hrem[idx].begin(), L1hrem[idx].end(), hs[j]) - 
                              L1hrem[idx].begin()] == -1)) {
                  print_buffer(buff_orig, MB_MESG_ENTS_LARGE, buffProcs[from_ind], false);
                  
                  std::cerr << "L1hrem[k], L1hloc[k], L1p[k]:" << std::endl;
                  for (unsigned int k = 0; k < L1hrem.size(); k++) {
                    std::cerr << procConfig.proc_rank() << ": k=" << k << ": " << L1hrem[idx][k] << ", " << L1hloc[idx][k] << ", " 
                              << L1p[idx][k] << std::endl;
                    std::cerr << procConfig.proc_rank() << ": ps, hs = ";
                    for (unsigned int l = 0; l < MAX_SHARING_PROCS && -1 != ps[l]; l++)
                      std::cerr << ps[l] << ", " << hs[l];
                    std::cerr << std::endl;
                  }
                }
#endif
                L1p[idx].push_back(-1);
                L1hrem[idx].push_back(hs[j]);
              }
              L1hloc[idx].push_back(new_h);
            }
          }
        }

        assert("Shouldn't be here for non-shared entities" &&
               -1 != num_ps);
        std::fill(&ps[0], &ps[num_ps], -1);
        std::fill(&hs[0], &hs[num_ps], 0);
      }
    }
    
    
#ifdef DEBUG_PACKING
      std::cerr << "Unpacked " << num_ents << " ents of type " 
                << CN::EntityTypeName(TYPE_FROM_HANDLE(this_type)) << std::endl;
#endif      

  }

#ifdef DEBUG_PACKING
  std::cerr << std::endl << "Done unpacking entities." << std::endl;
#endif

  return MB_SUCCESS;
}

ErrorCode ParallelComm::print_buffer(unsigned char *buff_ptr, 
                                         int mesg_tag, 
                                         int from_proc, bool sent) 
{
  std::cerr << procConfig.proc_rank();
  if (sent) std::cerr << " sent";
  else std::cerr << " received";
  std::cerr << " message type " << mesg_tag 
            << " to/from proc " << from_proc << "; contents:" << std::endl;

  int msg_length, num_ents;
  unsigned long extra;
  unsigned char *orig_ptr = buff_ptr;
  UNPACK_INT(buff_ptr, msg_length);
  std::cerr << msg_length << " bytes..." << std::endl;

  if (MB_MESG_ENTS_SIZE == mesg_tag || MB_MESG_ENTS_LARGE == mesg_tag) {

      // 1. # entities = E
    int i, j, k;
    std::vector<int> ps;
    std::vector<EntityHandle> hs;

    UNPACK_INT(buff_ptr, num_ents);
    std::cerr << num_ents << " entities..." << std::endl;

      // save place where remote handle info starts, then scan forward to ents
    for (i = 0; i < num_ents; i++) {
      UNPACK_INT(buff_ptr, j);
      if (0 > j) return MB_FAILURE;
      ps.resize(j);
      hs.resize(j);
      std::cerr << "Entity " << i << ", # procs = " << j << std::endl;
      UNPACK_INTS(buff_ptr, &ps[0], j);
      UNPACK_EH(buff_ptr, &hs[0], j);
      std::cerr << "   Procs: ";
      for (k = 0; k < j; k++) std::cerr << ps[k] << " ";
      std::cerr << std::endl;
      std::cerr << "   Handles: ";
      for (k = 0; k < j; k++) std::cerr << hs[k] << " ";
      std::cerr << std::endl;

      if (buff_ptr-orig_ptr > msg_length) {
        std::cerr << "End of buffer..." << std::endl;
        std::cerr.flush();
        return MB_FAILURE;
      }
    }
  
    while (true) {
      EntityType this_type = MBMAXTYPE;
      UNPACK_TYPE(buff_ptr, this_type);
      assert(this_type != MBENTITYSET);

        // MBMAXTYPE signifies end of entities data
      if (MBMAXTYPE == this_type) break;

        // get the number of ents
      int num_ents2, verts_per_entity;
      UNPACK_INT(buff_ptr, num_ents2);

        // unpack the nodes per entity
      if (MBVERTEX != this_type && num_ents2) {
        UNPACK_INT(buff_ptr, verts_per_entity);
      }

      std::cerr << "Type: " << CN::EntityTypeName(this_type)
                << "; num_ents = " << num_ents2;
      if (MBVERTEX != this_type) std::cerr << "; verts_per_ent = " << verts_per_entity;
      std::cerr << std::endl;
      if (num_ents2 < 0 || num_ents2 > msg_length) {
        std::cerr << "Wrong number of entities, returning." << std::endl;
        return MB_FAILURE;
      }
    
      unsigned long extra;
      for (int e = 0; e < num_ents2; e++) {
          // check for existing entity, otherwise make new one
        EntityHandle *connect;
        double *coords;

        ALIGN_BUFFER(buff_ptr, extra);
        if (MBVERTEX == this_type) {
          coords = (double*) buff_ptr;
          buff_ptr += 3*sizeof(double);
          std::cerr << "xyz = " << coords[0] << ", " << coords[1] << ", " 
                    << coords[2] << std::endl;
        }
        else {
          connect = (EntityHandle*) buff_ptr;
          buff_ptr += verts_per_entity * sizeof(EntityHandle);

            // update connectivity to local handles
          std::cerr << "Connectivity: ";
          for (k = 0; k < verts_per_entity; k++) std::cerr << connect[k] << " ";
          std::cerr << std::endl;
        }

        if (buff_ptr-orig_ptr > msg_length) {
          std::cerr << "End of buffer..." << std::endl;
          std::cerr.flush();
          return MB_FAILURE;
        }
      }
    }
  }
  
  else if (MB_MESG_REMOTEH_SIZE == mesg_tag || MB_MESG_REMOTEH_LARGE == mesg_tag) {
    UNPACK_INT(buff_ptr, num_ents);
    std::cerr << num_ents << " entities..." << std::endl;
    if (0 > num_ents || num_ents > msg_length) {
      std::cerr << "Wrong number of entities, returning." << std::endl;
      return MB_FAILURE;
    }
    std::vector<EntityHandle> L1hloc(num_ents), L1hrem(num_ents);
    std::vector<int> L1p(num_ents);
    UNPACK_INTS(buff_ptr, &L1p[0], num_ents);
    UNPACK_EH(buff_ptr, &L1hrem[0], num_ents);
    UNPACK_EH(buff_ptr, &L1hloc[0], num_ents);
    std::cerr << num_ents << " Entity pairs; hremote/hlocal/proc: " << std::endl;
    for (int i = 0; i < num_ents; i++) {
      EntityType etype = TYPE_FROM_HANDLE(L1hloc[i]);
      std::cerr << CN::EntityTypeName(etype) << ID_FROM_HANDLE(L1hrem[i])  << ", " 
                << CN::EntityTypeName(etype) << ID_FROM_HANDLE(L1hloc[i])  << ", " 
                << L1p[i] << std::endl;
    }

    if (buff_ptr-orig_ptr > msg_length) {
      std::cerr << "End of buffer..." << std::endl;
      std::cerr.flush();
      return MB_FAILURE;
    }

  }
  else if (mesg_tag == MB_MESG_TAGS_SIZE || mesg_tag == MB_MESG_TAGS_LARGE) {
    int num_tags, dum1, num_ents, data_type, tag_size;
    UNPACK_INT(buff_ptr, num_tags);
    std::cerr << "Number of tags = " << num_tags << std::endl;
    for (int i = 0; i < num_tags; i++) {
      std::cerr << "Tag " << i << ":" << std::endl;
      UNPACK_INT(buff_ptr, tag_size);
      UNPACK_INT(buff_ptr, dum1);
      UNPACK_INT(buff_ptr, data_type);
      std::cerr << "Tag size, type, data type = " << tag_size << ", " 
                << dum1 << ", " << data_type << std::endl;
      UNPACK_INT(buff_ptr, dum1);
      std::cerr << "Default value size = " << dum1 << std::endl;
      buff_ptr += dum1;
      UNPACK_INT(buff_ptr, dum1);
      std::cerr << "Tag name = " << (char*) buff_ptr << std::endl;
      buff_ptr += dum1;
      UNPACK_INT(buff_ptr, num_ents);
      std::cerr << "Number of ents = " << num_ents << std::endl;
      ALIGN_BUFFER(buff_ptr, extra);
      unsigned char *tmp_buff = buff_ptr;
      buff_ptr += num_ents*sizeof(EntityHandle);
      int tot_length = 0;
      for (int i = 0; i < num_ents; i++) {
        EntityType etype = TYPE_FROM_HANDLE(*((EntityHandle*)tmp_buff));
        std::cerr << CN::EntityTypeName(etype) << " " 
                  << ID_FROM_HANDLE(*((EntityHandle*)tmp_buff))
                  << ", tag = ";
        if (tag_size == MB_VARIABLE_LENGTH) {
          UNPACK_INT(buff_ptr, dum1);
          tot_length += dum1;
          std::cerr << "(variable, length = " << dum1 << ")" << std::endl;
        }
        else if (data_type == MB_TYPE_DOUBLE) std::cerr << *((double*)buff_ptr) << std::endl;
        else if (data_type == MB_TYPE_INTEGER) std::cerr << *((int*)buff_ptr) << std::endl;
        else if (data_type == MB_TYPE_OPAQUE) std::cerr << "(opaque)" << std::endl;
        else if (data_type == MB_TYPE_HANDLE) 
          std::cerr <<  (EntityHandle)*buff_ptr << std::endl;
        else if (data_type == MB_TYPE_BIT) std::cerr << "(bit)" << std::endl;
        tmp_buff += sizeof(EntityHandle);
        buff_ptr += tag_size;
      }

      if (tag_size == MB_VARIABLE_LENGTH) buff_ptr += tot_length;
    }
  }
  else {
    return MB_FAILURE;
  }

  std::cerr.flush();
  
  return MB_SUCCESS;
}

ErrorCode ParallelComm::list_entities(const EntityHandle *ents, int num_ents) 
{
  if (NULL == ents && 0 == num_ents) {
    sharedEnts.print("Shared entities:\n");
    return MB_SUCCESS;
  }
  
  else if (NULL == ents && 0 != num_ents) {
    return list_entities(sharedEnts);
  }
    
  unsigned char pstat;
  EntityHandle tmp_handles[MAX_SHARING_PROCS];
  int tmp_procs[MAX_SHARING_PROCS];
  unsigned int num_ps;
  ErrorCode result;

  for (int i = 0; i < num_ents; i++) {
    result = mbImpl->list_entities(ents+i, 1);

    result = get_sharing_data(ents[i], tmp_procs, tmp_handles, pstat, num_ps);
    RRA("Failed to get sharing data.");

    std::cout << "Pstatus: ";
    if (!num_ps)
      std::cout << "local " << std::endl;
    else {
      if (pstat & PSTATUS_NOT_OWNED) std::cout << "NOT_OWNED; ";
      if (pstat & PSTATUS_SHARED) std::cout << "SHARED; ";
      if (pstat & PSTATUS_MULTISHARED) std::cout << "MULTISHARED; ";
      if (pstat & PSTATUS_INTERFACE) std::cout << "INTERFACE; ";
      if (pstat & PSTATUS_GHOST) std::cout << "GHOST; ";
      std::cout << std::endl;
      for (unsigned int j = 0; j < num_ps; j++) {
        std::cout << "  proc " << tmp_procs[j] << " id (handle) " 
                  << mbImpl->id_from_handle(tmp_handles[j]) 
                  << "(" << tmp_handles[j] << ")" << std::endl;
      }
    }
    std::cout << std::endl;
  }

  return MB_SUCCESS;
}
  
ErrorCode ParallelComm::list_entities(const Range &ents) 
{
  for (Range::iterator rit = ents.begin(); rit != ents.end(); rit++)
    list_entities(&(*rit), 1);
  return MB_SUCCESS;
}

ErrorCode ParallelComm::update_remote_data(Range &local_range,
                                               Range &remote_range,
                                               int other_proc,
                                               const unsigned char add_pstat) 
{
  Range::iterator rit, rit2;
  ErrorCode result = MB_SUCCESS;

    // for each pair of local/remote handles:
  for (rit = local_range.begin(), rit2 = remote_range.begin(); 
       rit != local_range.end(); rit++, rit2++) {

    result = update_remote_data(*rit, &other_proc, &(*rit2), 1, add_pstat);
    RRA(" ");
  }

  return result;
}
  
ErrorCode ParallelComm::update_remote_data(const EntityHandle new_h,
                                               const int *ps,
                                               const EntityHandle *hs,
                                               const int num_ps,
                                               const unsigned char add_pstat) 
{
  EntityHandle tag_hs[MAX_SHARING_PROCS];
  int tag_ps[MAX_SHARING_PROCS];
  unsigned char pstat;
    // get initial sharing data; tag_ps and tag_hs get terminated with -1 and 0
    // in this function, so no need to initialize
  unsigned int num_exist;
  ErrorCode result = get_sharing_data(new_h, tag_ps, tag_hs, pstat, num_exist);
  RRA("");
  
#ifndef NDEBUG
  {
      // check for duplicates in proc list
    std::set<unsigned int> dumprocs;
    unsigned int dp = 0;
    for (; (int) dp < num_ps && -1 != ps[dp]; dp++)
      dumprocs.insert(ps[dp]);
    assert(dp == dumprocs.size());
  }
#endif      

    // add any new sharing data
  bool changed = false;
  int idx;
  if (!num_exist) {
      // just take what caller passed
    memcpy(tag_ps, ps, num_ps*sizeof(int));
    memcpy(tag_hs, hs, num_ps*sizeof(EntityHandle));
    num_exist = num_ps;
      // if it's only one, hopefully I'm not there yet...
    assert("I shouldn't be the only proc there." &&
           (1 != num_exist || ps[0] != (int)procConfig.proc_rank()));
    changed = true;
  }
  else {
    for (int i = 0; i < num_ps; i++) {
      idx = std::find(tag_ps, tag_ps+num_exist, ps[i]) - tag_ps;
      if (idx == (int) num_exist) {
          // if there's only 1 sharing proc, and it's not me, then
          // we'll end up with 3; add me to the front
        if (!i && num_ps == 1 && num_exist == 1 &&
            ps[0] != (int)procConfig.proc_rank()) {
          int j = 1;
            // if I own this entity, put me at front, otherwise after first
          if (!(pstat & PSTATUS_NOT_OWNED)) {
            tag_ps[1] = tag_ps[0];
            tag_hs[1] = tag_hs[0];
            j = 0;
          }
          tag_ps[j] = procConfig.proc_rank();
          tag_hs[j] = new_h;
          num_exist++;
        }
        
        tag_ps[num_exist] = ps[i];
        tag_hs[num_exist] = hs[i];
        num_exist++;
        changed = true;
      }
      else if (0 == tag_hs[idx]) {
        tag_hs[idx] = hs[i];
        changed = true;
      }
      else if (0 != hs[i]) {
        assert(hs[i] == tag_hs[idx]);
      }
    }
  }
  
    // adjust for interface layer if necessary
  if (add_pstat & PSTATUS_INTERFACE) {
    idx = std::min_element(tag_ps, tag_ps+num_exist) - tag_ps;
    if (idx) {
      int tag_proc = tag_ps[idx];
      tag_ps[idx] = tag_ps[0];
      tag_ps[0] = tag_proc;
      EntityHandle tag_h = tag_hs[idx];
      tag_hs[idx] = tag_hs[0];
      tag_hs[0] = tag_h;
      changed = true;
      if (tag_ps[0] != (int)procConfig.proc_rank()) pstat |= PSTATUS_NOT_OWNED;
    }
  }
    
  if (!changed) return MB_SUCCESS;
  
  assert("interface entities should have > 1 proc" &&
         (!(add_pstat & PSTATUS_INTERFACE) || num_exist > 1));
  assert("ghost entities should have > 1 proc" &&
         (!(add_pstat & PSTATUS_GHOST) || num_exist > 1));
  
    // if it's multi-shared and we created the entity in this unpack,
    // local handle probably isn't in handle list yet
  if (num_exist > 2) {
    idx = std::find(tag_ps, tag_ps+num_exist, procConfig.proc_rank()) - tag_ps;
    assert(idx < (int) num_exist);
    if (!tag_hs[idx])
      tag_hs[idx] = new_h;
  }
      
  int tag_p;
  EntityHandle tag_h;

    // reset single shared proc/handle if was shared and moving to multi-shared
  if (num_exist > 2 && !(pstat & PSTATUS_MULTISHARED) &&
      (pstat & PSTATUS_SHARED)) {
      // must remove sharedp/h first, which really means set to default value
    tag_p = -1;
    result = mbImpl->tag_set_data(sharedp_tag(), &new_h, 1, &tag_p);
    RRA("Couldn't set sharedp tag.");
    tag_h = 0;
    result = mbImpl->tag_set_data(sharedh_tag(), &new_h, 1, &tag_h);
    RRA("Couldn't set sharedh tag.");
  }

    // update pstat
  pstat |= add_pstat;
  
    // set sharing tags
  if (num_exist > 2) {
    std::fill(tag_ps+num_exist, tag_ps+MAX_SHARING_PROCS, -1);
    std::fill(tag_hs+num_exist, tag_hs+MAX_SHARING_PROCS, 0);
    result = mbImpl->tag_set_data(sharedps_tag(), &new_h, 1, tag_ps);
    RRA("Couldn't set sharedps tag.");
    result = mbImpl->tag_set_data(sharedhs_tag(), &new_h, 1, tag_hs);
    RRA("Couldn't set sharedhs tag.");
    pstat |= (PSTATUS_MULTISHARED | PSTATUS_SHARED);

#ifndef NDEBUG
    {
        // check for duplicates in proc list
      std::set<unsigned int> dumprocs;
      unsigned int dp = 0;
      for (; dp < num_exist && -1 != tag_ps[dp]; dp++)
        dumprocs.insert(tag_ps[dp]);
      assert(dp == dumprocs.size());
    }
#endif      
  }
  else if (num_exist == 2 || num_exist == 1) {
    if (tag_ps[0] == (int) procConfig.proc_rank()) {
      assert(2 == num_exist && tag_ps[1] != (int) procConfig.proc_rank());
      tag_ps[0] = tag_ps[1];
      tag_hs[0] = tag_hs[1];
    }
    assert(tag_ps[0] != -1 && tag_hs[0] != 0);
    result = mbImpl->tag_set_data(sharedp_tag(), &new_h, 1, tag_ps);
    RRA("Couldn't set sharedp tag.");
    result = mbImpl->tag_set_data(sharedh_tag(), &new_h, 1, tag_hs);
    RRA("Couldn't set sharedh tag.");
    pstat |= PSTATUS_SHARED;
  }

    // now set new pstatus
  result = mbImpl->tag_set_data(pstatus_tag(), &new_h, 1, &pstat);
  RRA("Couldn't set pstatus tag.");

  if (pstat & PSTATUS_SHARED) sharedEnts.insert(new_h);
  
  return MB_SUCCESS;
}

ErrorCode ParallelComm::get_sharing_data(const Range &entities,
                                             std::set<int> &procs,
                                             int operation)
{
    // get the union or intersection of sharing data for multiple entities

  ErrorCode result;
  int sp2[MAX_SHARING_PROCS];
  int num_ps;
  unsigned char pstat;
  std::set<int> tmp_procs;
  procs.clear();
  
  for (Range::const_iterator rit = entities.begin(); rit != entities.end(); rit++) {
        
      // get sharing procs
    result = get_sharing_data(*rit, sp2, NULL, pstat, num_ps);
    RRA("Problem getting sharing data in get_sharing_data.");
    if (!(pstat & PSTATUS_SHARED) && Interface::INTERSECT == operation) {
      procs.clear();
      return MB_SUCCESS;
    }
        
    if (rit == entities.begin()) {
      std::copy(sp2, sp2+num_ps, std::inserter(procs, procs.begin()));
    }
    else {
      std::sort(sp2, sp2+num_ps);
      tmp_procs.clear();
      if (Interface::UNION == operation) 
        std::set_union(procs.begin(), procs.end(), 
                       sp2, sp2+num_ps, std::inserter(tmp_procs, tmp_procs.end()));
      else if (Interface::INTERSECT == operation)
        std::set_intersection(procs.begin(), procs.end(), 
                              sp2, sp2+num_ps, std::inserter(tmp_procs, tmp_procs.end()));
      else {
        assert("Unknown operation." && false);
        return MB_FAILURE;
      }
      procs.swap(tmp_procs);
    }
    if (Interface::INTERSECT == operation && procs.empty()) 
      return MB_SUCCESS;
  }

  return MB_SUCCESS;
}
  
ErrorCode ParallelComm::get_sharing_data(const EntityHandle entity,
                                             int *ps, 
                                             EntityHandle *hs,
                                             unsigned char &pstat,
                                             unsigned int &num_ps)
{
  ErrorCode result = mbImpl->tag_get_data(pstatus_tag(), &entity, 1, &pstat);
  RRA("Couldn't get pstatus tag.");
  if (pstat & PSTATUS_MULTISHARED) {
    result = mbImpl->tag_get_data(sharedps_tag(), &entity, 1, ps);
    RRA("Couldn't get sharedps tag.");
    if (hs) {
      result = mbImpl->tag_get_data(sharedhs_tag(), &entity, 1, hs);
      RRA("Couldn't get sharedhs tag.");
    }
    num_ps = std::find(ps, ps+MAX_SHARING_PROCS, -1) - ps;
  }
  else if (pstat & PSTATUS_SHARED) {
    result = mbImpl->tag_get_data(sharedp_tag(), &entity, 1, ps);
    RRA("Couldn't get sharedp tag.");
    if (hs) {
      result = mbImpl->tag_get_data(sharedh_tag(), &entity, 1, hs);
      RRA("Couldn't get sharedh tag.");
      hs[1] = 0;
    }
      // initialize past end of data
    ps[1] = -1;
    num_ps = 1;
  }
  else {
    ps[0] = -1;
    if (hs) hs[0] = 0;
    num_ps = 0;
  }

  assert(0 <= num_ps && MAX_SHARING_PROCS >= num_ps);
  
  return MB_SUCCESS;
}
  
ErrorCode ParallelComm::find_existing_entity(const bool is_iface,
                                                 const int owner_p,
                                                 const EntityHandle owner_h,
                                                 const int num_ps,
                                                 const EntityHandle *connect,
                                                 const int num_connect,
                                                 const EntityType this_type,
                                                 std::vector<EntityHandle> &L2hloc,
                                                 std::vector<EntityHandle> &L2hrem,
                                                 std::vector<unsigned int> &L2p,
                                                 EntityHandle &new_h) 
{
  new_h = 0;
  if (!is_iface && num_ps > 2) {
    for (unsigned int i = 0; i < L2hrem.size(); i++) {
      if (L2hrem[i] == owner_h && owner_p == (int) L2p[i]) {
        new_h = L2hloc[i];
        return MB_SUCCESS;
      }
    }        
  }

    // if we got here and it's a vertex, we don't need to look further
  if (MBVERTEX == this_type || !connect || !num_connect) return MB_SUCCESS;
  
  Range tmp_range;
  ErrorCode result = mbImpl->get_adjacencies(connect, num_connect, 
                                               CN::Dimension(this_type), false, 
                                               tmp_range);
  RRA("Problem getting existing entity.");
  if (!tmp_range.empty()) {
      // found a corresponding entity - return target
    new_h = *tmp_range.begin();
  }  
  else {
    new_h = 0;
  }

  return MB_SUCCESS;
}

ErrorCode ParallelComm::get_local_handles(const Range &remote_handles,
                                              Range &local_handles,
                                              const Range &new_ents) 
{
  std::vector<EntityHandle> rh_vec;
  rh_vec.reserve(remote_handles.size());
  std::copy(remote_handles.begin(), remote_handles.end(), std::back_inserter(rh_vec));
  ErrorCode result = get_local_handles(&rh_vec[0], remote_handles.size(), new_ents);
  std::copy(rh_vec.begin(), rh_vec.end(), range_inserter(local_handles));
  return result;
}
  
ErrorCode ParallelComm::get_local_handles(EntityHandle *from_vec, 
                                              int num_ents,
                                              const Range &new_ents) 
{
  std::vector<EntityHandle> tmp_ents;
  std::copy(new_ents.begin(), new_ents.end(), std::back_inserter(tmp_ents));
  return get_local_handles(from_vec, num_ents, tmp_ents);
}

ErrorCode ParallelComm::get_local_handles(EntityHandle *from_vec,
                                              int num_ents,
                                              const std::vector<EntityHandle> &new_ents) 
{
  for (int i = 0; i < num_ents; i++) {
    if (TYPE_FROM_HANDLE(from_vec[i]) == MBMAXTYPE) {
      assert(ID_FROM_HANDLE(from_vec[i]) < (int) new_ents.size());
      from_vec[i] = new_ents[ID_FROM_HANDLE(from_vec[i])];
    }
  }
  
  return MB_SUCCESS;
}

template <typename T> void
insert_in_array( T* array, size_t array_size, size_t location, T value )
{
  assert( location+1 < array_size );
  for (size_t i = array_size-1; i > location; --i)
    array[i] = array[i-1];
  array[location] = value;
}

ErrorCode ParallelComm::pack_range_map(Range &key_range, EntityHandle val_start,
                                           HandleMap &handle_map) 
{
  for (Range::const_pair_iterator key_it = key_range.const_pair_begin(); 
       key_it != key_range.const_pair_end(); key_it++) {
    int tmp_num = (*key_it).second - (*key_it).first + 1;
    handle_map.insert((*key_it).first, val_start, tmp_num);
    val_start += tmp_num;
  }

  return MB_SUCCESS;
}

ErrorCode ParallelComm::pack_sets(Range &entities,
                                      Buffer *buff,
                                      const bool store_remote_handles,
                                      const int to_proc)
{
    // SETS:
    // . #sets
    // . for each set:
    //   - options[#sets] (unsigned int)
    //   - if (unordered) set range 
    //   - else if ordered
    //     . #ents in set
    //     . handles[#ents]
    //   - #parents
    //   - if (#parents) handles[#parents]
    //   - #children
    //   - if (#children) handles[#children]
  
    // now the sets; assume any sets the application wants to pass are in the entities list
  ErrorCode result;
  Range all_sets = entities.subset_by_type(MBENTITYSET);

  unsigned int num_aligns;
  int buff_size = estimate_sets_buffer_size(all_sets, store_remote_handles, num_aligns);
  buff->check_space(buff_size);

    // number of sets
  PACK_INT(buff->buff_ptr, all_sets.size());

    // options for all sets
  std::vector<unsigned int> options(all_sets.size());
  Range::iterator rit;
  std::vector<EntityHandle> members;
  int i;
  for (rit = all_sets.begin(), i = 0; rit != all_sets.end(); rit++, i++) {
      result = mbImpl->get_meshset_options(*rit, options[i]);
      RRA("Failed to get meshset options.");
  }
  buff->check_space(all_sets.size()*sizeof(unsigned int));
  PACK_VOID(buff->buff_ptr, &options[0], all_sets.size()*sizeof(unsigned int));
  
    // vectors/ranges
  for (rit = all_sets.begin(), i = 0; rit != all_sets.end(); rit++, i++) {
      Range set_range;
      if (options[i] & MESHSET_SET) {
        Range set_range;
        result = mbImpl->get_entities_by_handle(*rit, set_range);
        RRA("Failed to get set entities.");

        buff_size = RANGE_SIZE(set_range);
        buff->check_space(buff_size);
        PACK_RANGE(buff->buff_ptr, set_range);
      }
      else if (options[i] & MESHSET_ORDERED) {
        members.clear();
        result = mbImpl->get_entities_by_handle(*rit, members);
        RRA("Failed to get entities in ordered set.");
        
        buff->check_space(members.size()*sizeof(EntityHandle)+sizeof(int));
        PACK_INT(buff->buff_ptr, members.size());
          // check size, since calling pack aligns the buffer
        if (members.size())
          PACK_EH(buff->buff_ptr, &members[0], members.size());
      }
  }
    // pack numbers of parents/children
  unsigned int tot_pch = 0;
  int num_pch;
  buff->check_space(2*all_sets.size()*sizeof(int));
  for (rit = all_sets.begin(), i = 0; rit != all_sets.end(); rit++, i++) {
      // pack parents
    result = mbImpl->num_parent_meshsets(*rit, &num_pch);
    RRA("Failed to get num parents.");
    PACK_INT(buff->buff_ptr, num_pch);
    tot_pch += num_pch;
    result = mbImpl->num_child_meshsets(*rit, &num_pch);
    RRA("Failed to get num children.");
    PACK_INT(buff->buff_ptr, num_pch);
    tot_pch += num_pch;
  }

    // now pack actual parents/children
  members.clear();
  members.reserve(tot_pch);
  std::vector<EntityHandle> tmp_pch;
  for (rit = all_sets.begin(), i = 0; rit != all_sets.end(); rit++, i++) {
    result = mbImpl->get_parent_meshsets(*rit, tmp_pch);
    RRA("Failed to get parents.");
    std::copy(tmp_pch.begin(), tmp_pch.end(), std::back_inserter(members));
    tmp_pch.clear();
    result = mbImpl->get_child_meshsets(*rit, tmp_pch);
    RRA("Failed to get children.");
    std::copy(tmp_pch.begin(), tmp_pch.end(), std::back_inserter(members));
    tmp_pch.clear();
  }
  assert(members.size() == tot_pch);
  if (!members.empty()) {
    result = get_remote_handles(store_remote_handles,
                                &members[0], &members[0], 
                                members.size(), to_proc,
                                entities);
    RRA("Trouble getting remote handles for set parent/child sets.");
#ifndef NDEBUG
      // check that all handles are either sets or maxtype
    for (unsigned int __j = 0; __j < members.size(); __j++)
      assert((TYPE_FROM_HANDLE(members[__j]) == MBMAXTYPE &&
              ID_FROM_HANDLE(members[__j]) < (int)entities.size()) ||
             TYPE_FROM_HANDLE(members[__j]) == MBENTITYSET);
#endif        
    buff->check_space(members.size()*sizeof(EntityHandle));
    PACK_EH(buff->buff_ptr, &members[0], members.size());
  }
    
    // pack the handles
  if (store_remote_handles && !all_sets.empty()) {
    buff_size = RANGE_SIZE(all_sets);
    buff->check_space(buff_size);
    PACK_RANGE(buff->buff_ptr, all_sets);
  }
  
#ifdef DEBUG_PACKING
  std::cerr << std::endl << "Done packing sets." << std::endl;
#endif

  buff->set_stored_size();
  
  return MB_SUCCESS;
}

ErrorCode ParallelComm::unpack_sets(Buffer *buff,
                                        Range &entities,
                                        const bool store_remote_handles,
                                        const int from_proc)
{
  
    // now the sets; assume any sets the application wants to pass are in the entities list
  ErrorCode result;

  Range new_sets;
  int num_sets;
  UNPACK_INT(buff->buff_ptr, num_sets);

  if (!num_sets) return MB_SUCCESS;
         
  std::vector<EntityHandle> members;
  int num_ents;
  std::vector<unsigned int> options_vec(num_sets);
      // option value
  if (num_sets)
    UNPACK_VOID(buff->buff_ptr, &options_vec[0], num_sets*sizeof(unsigned int));

    // create sets
  int i;
  Range::const_iterator rit;
  for (i = 0; i < num_sets; i++) {
    
      // create the set
    EntityHandle set_handle;
    result = mbImpl->create_meshset(options_vec[i], set_handle);
    RRA("Failed to create set in unpack.");

    // make sure new sets handles are monotonically increasing
    assert(set_handle > *new_sets.rbegin());

    new_sets.insert(set_handle);
  }

  entities.merge(new_sets);
  
  for (rit = new_sets.begin(), i = 0; rit != new_sets.end(); rit++, i++) {
    if (options_vec[i] & MESHSET_SET) {
        // unpack entities as a range
      Range set_range, tmp_range;
      UNPACK_RANGE(buff->buff_ptr, tmp_range);
      result = get_local_handles(tmp_range, set_range, entities);      
      RRA("Failed to get local handles for unordered set contents.");
      result = mbImpl->add_entities(*rit, set_range);
      RRA("Failed to add ents to unordered set in unpack.");
    }
    else if (options_vec[i] & MESHSET_ORDERED) {
        // unpack entities as vector, with length
      UNPACK_INT(buff->buff_ptr, num_ents);
      members.resize(num_ents);
      if (num_ents) UNPACK_EH(buff->buff_ptr, &members[0], num_ents);
      result = get_local_handles(&members[0], num_ents, entities);
      RRA("Failed to get local handles for ordered set contents.");
      result = mbImpl->add_entities(*rit, &members[0], num_ents);
      RRA("Failed to add ents to ordered set in unpack.");
    }
  }

  std::vector<int> num_pch(2*new_sets.size());
  std::vector<int>::iterator vit;
  int tot_pch = 0;
  for (vit = num_pch.begin(); vit != num_pch.end(); vit++) {
    UNPACK_INT(buff->buff_ptr, *vit);
    tot_pch += *vit;
  }
  
  members.resize(tot_pch);
  UNPACK_EH(buff->buff_ptr, &members[0], tot_pch);
  result = get_local_handles(&members[0], tot_pch, entities);
  RRA("Couldn't get local handle for parent/child sets.");

  int num = 0;
  EntityHandle *mem_ptr = &members[0];
  for (rit = new_sets.begin(); rit != new_sets.end(); rit++) {
      // unpack parents/children
    int num_par = num_pch[num++], num_child = num_pch[num++];
    if (num_par+num_child) {
      for (i = 0; i < num_par; i++) {
        assert(0 != mem_ptr[i]);
        result = mbImpl->add_parent_meshset(*rit, mem_ptr[i]);
        RRA("Failed to add parent to set in unpack.");
      }
      mem_ptr += num_par;
      for (i = 0; i < num_child; i++) {
        assert(0 != mem_ptr[i]);
        result = mbImpl->add_child_meshset(*rit, mem_ptr[i]);
        RRA("Failed to add child to set in unpack.");
      }
      mem_ptr += num_child;
    }
  }

    // unpack source handles
  Range dum_range;
  if (store_remote_handles && !new_sets.empty()) {
    UNPACK_RANGE(buff->buff_ptr, dum_range);
    result = update_remote_data(new_sets, dum_range, from_proc, 0);
    RRA("Couldn't set sharing data for sets");
  }

#ifdef DEBUG_PACKING
  std::cerr << std::endl << "Done unpacking sets." << std::endl;
#endif

  return MB_SUCCESS;
}

ErrorCode ParallelComm::pack_adjacencies(Range &entities,
                                         Range::const_iterator &start_rit,
                                         Range &whole_range,
                                         Buffer *buff,
                                         int &count,
                                         const bool just_count,
                                         const bool store_handles,
                                         const int to_proc)
{
  return MB_FAILURE;
}

ErrorCode ParallelComm::unpack_adjacencies(Buffer *buff,
                                           Range &entities,
                                               const bool store_handles,
                                               const int from_proc)
{
  return MB_FAILURE;
}

ErrorCode ParallelComm::pack_tags(Range &entities,
                                      const std::vector<Tag> &src_tags,
                                      const std::vector<Tag> &dst_tags,
                                      const std::vector<Range> &tag_ranges,
                                      Buffer *buff,
                                      const bool store_remote_handles,
                                      const int to_proc)
{
  

  ErrorCode result;
  std::vector<Tag>::const_iterator tag_it, dst_it;
  std::vector<Range>::const_iterator rit;
  int count = 0;
  
  for (tag_it = src_tags.begin(), rit = tag_ranges.begin(); 
       tag_it != src_tags.end(); tag_it++, rit++) {

    result = packed_tag_size( *tag_it, *rit, count );
    if (MB_SUCCESS != result)
      return result;
  }
    
    // number of tags
  count += sizeof(int);

  buff->check_space(count);
  
  PACK_INT(buff->buff_ptr, src_tags.size());
    
  for (tag_it = src_tags.begin(), dst_it = dst_tags.begin(), rit = tag_ranges.begin(); 
       tag_it != src_tags.end(); tag_it++, dst_it++, rit++) {
    
    result = pack_tag( *tag_it, *dst_it, *rit, entities, buff,
                       store_remote_handles, to_proc );
    if (MB_SUCCESS != result)
      return result;
  }
  
#ifdef DEBUG_PACKING
  std::cerr << std::endl << "Done packing tags." << std::endl;
#endif

  buff->set_stored_size();
  
  return MB_SUCCESS;
}
         

ErrorCode ParallelComm::packed_tag_size( Tag tag,
                                             const Range &tagged_entities,
                                             int &count )
{
    // for dense tags, compute size assuming all entities have that tag
    // for sparse tags, get number of entities w/ that tag to compute size

  std::vector<int> var_len_sizes;
  std::vector<const void*> var_len_values;
    
  const TagInfo *tinfo = tagServer->get_tag_info(tag);
    // default value
  count += sizeof(int);
  if (NULL != tinfo->default_value()) 
    count += tinfo->default_value_size();

    // size, type, data type
  count += 3*sizeof(int);

    // name
  count += sizeof(int);
  count += tinfo->get_name().size();

    // range of tag
  count += sizeof(int) + tagged_entities.size() * sizeof(EntityHandle);

  if (tinfo->get_size() == MB_VARIABLE_LENGTH) {
    const int num_ent = tagged_entities.size();
      // send a tag size for each entity
    count += num_ent * sizeof(int);
      // send tag data for each entity
    var_len_sizes.resize( num_ent );
    var_len_values.resize( num_ent );
    ErrorCode result = tagServer->get_data( tag,
                                              tagged_entities, 
                                              &var_len_values[0], 
                                              &var_len_sizes[0] );
    RRA("Failed to get lenghts of variable-length tag values.");
    count += std::accumulate( var_len_sizes.begin(), var_len_sizes.end(), 0 );
  }
  else {
      // tag data values for range or vector
    count += tagged_entities.size() * tinfo->get_size();
  }
  
  return MB_SUCCESS;
}


ErrorCode ParallelComm::pack_tag( Tag src_tag,
                                      Tag dst_tag,
                                      const Range &tagged_entities,
                                      const Range &whole_range,
                                      Buffer *buff,
                                      const bool store_remote_handles,
                                      const int to_proc )
{
  ErrorCode result;
  std::vector<int> var_len_sizes;
  std::vector<const void*> var_len_values;
  unsigned long extra;

  const TagInfo* tinfo = tagServer->get_tag_info(src_tag);
  if (!tinfo)
    return MB_TAG_NOT_FOUND;
    
  const TagInfo* dst_tinfo;
  if (src_tag == dst_tag) {
    dst_tinfo = tinfo;
  }
  else {
    dst_tinfo = tagServer->get_tag_info(dst_tag);
    if (!dst_tinfo)
      return MB_TAG_NOT_FOUND;
    if (dst_tinfo->get_size() != tinfo->get_size())
      return MB_TYPE_OUT_OF_RANGE;
    if (dst_tinfo->get_data_type() != tinfo->get_data_type() && 
        dst_tinfo->get_data_type() != MB_TYPE_OPAQUE &&
            tinfo->get_data_type() != MB_TYPE_OPAQUE)
      return MB_TYPE_OUT_OF_RANGE;
  }
    
    // size, type, data type
  buff->check_space(3*sizeof(int));
  PACK_INT(buff->buff_ptr, tinfo->get_size());
  TagType this_type;
  result = mbImpl->tag_get_type(dst_tag, this_type);
  PACK_INT(buff->buff_ptr, (int)this_type);
  PACK_INT(buff->buff_ptr, (int)(tinfo->get_data_type()));

    // default value
  if (NULL == tinfo->default_value()) {
    buff->check_space(sizeof(int));
    PACK_INT(buff->buff_ptr, 0);
  }
  else {
    PACK_INT(buff->buff_ptr, tinfo->get_size());
    if (tinfo->get_data_type() == MB_TYPE_DOUBLE ||
        tinfo->get_data_type() == MB_TYPE_HANDLE) {
      ALIGN_BUFFER(buff->buff_ptr, extra);
    }
    buff->check_space(tinfo->default_value_size());
    PACK_BYTES(buff->buff_ptr, tinfo->default_value(), tinfo->default_value_size());
  }

    // name
  buff->check_space(tinfo->get_name().size());
  PACK_BYTES(buff->buff_ptr, dst_tinfo->get_name().c_str(), dst_tinfo->get_name().size());

#ifdef DEBUG_PACKING
  std::cerr << "Packing tag \"" << tinfo->get_name() << "\"";
  if (tinfo != dst_tinfo)
    std::cerr << " (as tag \"" << dst_tinfo->get_name() << "\")";
  std::cerr << std::endl;
#endif    
    // pack entities
  buff->check_space(tagged_entities.size()*sizeof(EntityHandle)+sizeof(int));
  PACK_INT(buff->buff_ptr, tagged_entities.size());
  ALIGN_BUFFER(buff->buff_ptr, extra);
  result = get_remote_handles(store_remote_handles,
                              tagged_entities, (EntityHandle*)buff->buff_ptr, to_proc,
                              whole_range);
#ifdef DEBUG_PACKING
  if (MB_SUCCESS != result) {
    std::cerr << "Trouble getting remote handles for tagged entities:" << std::endl;
    tagged_entities.print("  ");
  }
#else
  RRA("Trouble getting remote handles for tagged entities.");
#endif

  buff->buff_ptr += tagged_entities.size() * sizeof(EntityHandle);

  const size_t num_ent = tagged_entities.size();
  if (tinfo->get_size() == MB_VARIABLE_LENGTH) {
    var_len_sizes.resize( num_ent, 0 );
    var_len_values.resize( num_ent, 0 );
    result = mbImpl->tag_get_data(src_tag, tagged_entities, &var_len_values[0], 
                                  &var_len_sizes[0] );
    RRA("Failed to get variable-length tag data in pack_tags.");
    buff->check_space(num_ent*sizeof(int));
    PACK_INTS(buff->buff_ptr, &var_len_sizes[0], num_ent);
    for (unsigned int i = 0; i < num_ent; ++i) {
      buff->check_space(var_len_sizes[i]);
      PACK_VOID(buff->buff_ptr, var_len_values[i], var_len_sizes[i]);
    }
  }
  else {
    if (tinfo->get_data_type() == MB_TYPE_DOUBLE ||
        tinfo->get_data_type() == MB_TYPE_HANDLE) {
      ALIGN_BUFFER(buff->buff_ptr, extra);
    }
    buff->check_space(num_ent * tinfo->get_size());
    result = mbImpl->tag_get_data(src_tag, tagged_entities, buff->buff_ptr);
    RRA("Failed to get tag data in pack_tags.");
    buff->buff_ptr += num_ent * tinfo->get_size();
    PC(num_ent*tinfo->get_size(), " void");
  }

  return MB_SUCCESS;
}

ErrorCode ParallelComm::get_tag_send_list( const Range& whole_range,
                                               std::vector<Tag>& all_tags,
                                               std::vector<Range>& tag_ranges )
{
  std::vector<Tag> tmp_tags;
  ErrorCode result = tagServer->get_tags(tmp_tags);
  RRA("Failed to get tags in pack_tags.");

  std::vector<Tag>::iterator tag_it;
  for (tag_it = tmp_tags.begin(); tag_it != tmp_tags.end(); tag_it++) {
    std::string tag_name;
    result = mbImpl->tag_get_name(*tag_it, tag_name);
    if (tag_name.c_str()[0] == '_' && tag_name.c_str()[1] == '_')
      continue;

    Range tmp_range;
    result = tagServer->get_entities(*tag_it, tmp_range);
    RRA("Failed to get entities for tag in pack_tags.");
    tmp_range = intersect( tmp_range, whole_range);

    if (tmp_range.empty()) continue;
        
      // ok, we'll be sending this tag
    all_tags.push_back( *tag_it );
    tag_ranges.push_back( Range() );
    tag_ranges.back().swap( tmp_range );
  }
  
  return MB_SUCCESS;
}



ErrorCode ParallelComm::unpack_tags(Buffer *buff,
                                    Range &entities,
                                        const bool store_remote_handles,
                                        const int from_proc)
{
    // tags
    // get all the tags
    // for dense tags, compute size assuming all entities have that tag
    // for sparse tags, get number of entities w/ that tag to compute size

  ErrorCode result;
  
  int num_tags = 0;
  UNPACK_INT(buff->buff_ptr, num_tags);
  std::vector<EntityHandle> tag_ents;
  std::vector<const void*> var_len_vals;
  std::vector<int> var_lengths;
  unsigned long extra;

  for (int i = 0; i < num_tags; i++) {
    
        // tag handle
    Tag tag_handle;

      // size, data type
    int tag_size, tag_data_type, tag_type;
    UNPACK_INT(buff->buff_ptr, tag_size);
    UNPACK_INT(buff->buff_ptr, tag_type);
    UNPACK_INT(buff->buff_ptr, tag_data_type);
      
      // default value
    int def_val_size = 0;
    UNPACK_INT(buff->buff_ptr, def_val_size);
    
    void *def_val_ptr = NULL;
    if (def_val_size) {
      if (MB_TYPE_DOUBLE == tag_data_type || MB_TYPE_HANDLE == tag_data_type) {
        ALIGN_BUFFER(buff->buff_ptr, extra);
      }
      def_val_ptr = buff->buff_ptr;
      buff->buff_ptr += def_val_size;
      UPC(tag_size, " void");
    }
    
      // name
    int name_len;
    UNPACK_INT(buff->buff_ptr, name_len);
    std::string tag_name( reinterpret_cast<char*>(buff->buff_ptr), name_len );
    buff->buff_ptr += name_len;
    UPC(64, " chars");
#ifdef DEBUG_PACKING
    std::cerr << "Unpacking tag " << tag_name << std::endl;
#endif    

      // create the tag
    if (tag_size == MB_VARIABLE_LENGTH) 
      result = mbImpl->tag_create_variable_length( tag_name.c_str(), (TagType)tag_type,
                                                   (DataType)tag_data_type, tag_handle,
                                                   def_val_ptr, def_val_size );
    else
      result = mbImpl->tag_create(tag_name.c_str(), tag_size, (TagType) tag_type, 
                                  (DataType) tag_data_type, tag_handle,
                                  def_val_ptr);
    if (MB_ALREADY_ALLOCATED == result) {
        // already allocated tag, check to make sure it's the same size, type, etc.
      const TagInfo *tag_info = tagServer->get_tag_info(tag_name.c_str());
      TagType this_type;
      result = mbImpl->tag_get_type(tag_handle, this_type);
      if (tag_size != tag_info->get_size() ||
          tag_type != this_type ||
          tag_data_type != tag_info->get_data_type() ||
          (def_val_ptr && !tag_info->default_value()) ||
          (!def_val_ptr && tag_info->default_value())) {
        RRA("Didn't get correct tag info when unpacking tag.");
      }
    }
    else if (MB_SUCCESS != result) return result;

      // go through handle vec (in buffer) and convert to local handles in-place
    int num_ents;
    UNPACK_INT(buff->buff_ptr, num_ents);
    ALIGN_BUFFER(buff->buff_ptr, extra);
    EntityHandle *handle_vec = (EntityHandle*)buff->buff_ptr;
    buff->buff_ptr += num_ents * sizeof(EntityHandle);

    if (!store_remote_handles) {
        // in this case handles are indices into new entity range; need to convert
        // to local handles
      result = get_local_handles(handle_vec, num_ents, entities);
      RRA("Unable to convert to local handles.");
    }

      // if it's a handle type, also convert tag vals in-place in buffer
    if (MB_TYPE_HANDLE == tag_type) {
      EntityHandle *val_vec = (EntityHandle*)buff->buff_ptr;
      result = get_local_handles(val_vec, num_ents, entities);
      RRA("Failed to get local handles for tag vals.");
    }

    if (tag_size == MB_VARIABLE_LENGTH) {
        // Be careful of alignment here.  If the integers are aligned
        // in the buffer, we can use them directly.  Otherwise we must
        // copy them.
      const int* size_arr;
      if (((size_t)buff->buff_ptr)%4) {
        var_lengths.resize( num_ents );
        memcpy( &var_lengths[0], buff->buff_ptr, num_ents*sizeof(int) );
        size_arr = &var_lengths[0];
      }
      else {
        size_arr = reinterpret_cast<const int*>(buff->buff_ptr);
      }
      buff->buff_ptr += sizeof(int) * num_ents;
      UPC(sizeof(int) * num_ents, " void");
      
        // get pointers into buffer for each tag value
      var_len_vals.resize(num_ents);
      for (std::vector<EntityHandle>::size_type i = 0; 
           i < (std::vector<EntityHandle>::size_type) num_ents; ++i) {
        var_len_vals[i] = buff->buff_ptr;
        buff->buff_ptr += size_arr[i];
        UPC(size_arr[i], " void");
      }
      result = mbImpl->tag_set_data( tag_handle, handle_vec, num_ents,
                                     &var_len_vals[0], size_arr );
      RRA("Trouble setting tag data when unpacking variable-length tag.");
    }
    else {
      result = mbImpl->tag_set_data(tag_handle, handle_vec,
                                    num_ents, buff->buff_ptr);
      RRA("Trouble setting range-based tag data when unpacking tag.");
      buff->buff_ptr += num_ents * tag_size;
      UPC(num_ents * tag_size, " void");
    }
  }
  
#ifdef DEBUG_PACKING
  std::cerr << std::endl << "Done unpacking tags." << std::endl;
#endif

  return MB_SUCCESS;
}

ErrorCode ParallelComm::resolve_shared_ents(EntityHandle this_set,
                                                int resolve_dim,
                                                int shared_dim,
                                                const Tag* id_tag) 
{
  ErrorCode result;
  Range proc_ents;
      // get the entities in the partition sets
  for (Range::iterator rit = partitionSets.begin(); rit != partitionSets.end(); rit++) {
    Range tmp_ents;
    result = mbImpl->get_entities_by_handle(*rit, tmp_ents, true);
    if (MB_SUCCESS != result) return result;
    proc_ents.merge(tmp_ents);
  }

    // resolve dim is maximal dim of entities in proc_ents
  if (-1 == resolve_dim) {
    if (proc_ents.empty()) 
      return MB_ENTITY_NOT_FOUND;

    resolve_dim = mbImpl->dimension_from_handle(*proc_ents.rbegin()); 
  }

    // proc_ents should all be of same dimension
  if (resolve_dim > shared_dim &&
      mbImpl->dimension_from_handle(*proc_ents.rbegin()) !=
      mbImpl->dimension_from_handle(*proc_ents.begin())) {
    Range::iterator lower = proc_ents.lower_bound(CN::TypeDimensionMap[0].first),
      upper = proc_ents.upper_bound(CN::TypeDimensionMap[resolve_dim-1].second);
    proc_ents.erase(lower, upper);
  }
  
    // must call even if we don't have any entities, to make sure
    // collective comm'n works
  return resolve_shared_ents(this_set, proc_ents, resolve_dim, shared_dim, id_tag);
}
  
ErrorCode ParallelComm::resolve_shared_ents(EntityHandle this_set,
                                                Range &proc_ents,
                                                int resolve_dim,
                                                int shared_dim,
                                                const Tag* id_tag) 
{
#ifdef DEBUG_MPE
  define_mpe();

  MPE_Log_event(RESOLVE_START, procConfig.proc_rank(), "Entering resolve_shared_ents.");
#endif

  ErrorCode result;
  if (debug) std::cerr << "Resolving shared entities." << std::endl;

  if (-1 == shared_dim) {
    if (0 == resolve_dim) {
      result = mbImpl->get_dimension(shared_dim); 
      RRA("Couldn't get dimension.");
    }
    else if (!proc_ents.empty())
      shared_dim = mbImpl->dimension_from_handle(*proc_ents.begin())-1;
    else if (resolve_dim == 3)
      shared_dim = 2;
    else {
      assert(false && "Unable to guess shared_dim.");
      return MB_FAILURE;
    }
  }
  assert(shared_dim >= 0 && resolve_dim >= 0);
  
    // get the skin entities by dimension
  Range skin_ents[4];
  std::vector<int> gid_data;
  std::vector<EntityHandle> handle_vec;
  int skin_dim;

    // get the entities to be skinned
  if (resolve_dim < shared_dim) {
      // for vertex-based partition, it's the elements adj to the vertices
    result = mbImpl->get_adjacencies(proc_ents, shared_dim,
                                     false, skin_ents[resolve_dim],
                                     Interface::UNION);
    RRA("Failed getting skinned entities.");
    skin_dim = shared_dim-1;
  }
  else {
      // for element-based partition, it's just the elements
    skin_ents[resolve_dim] = proc_ents;
    skin_dim = resolve_dim-1;
  }

    // find the skin
  Skinner skinner(mbImpl);
  result = skinner.find_skin(skin_ents[skin_dim+1], false, skin_ents[skin_dim],
                             0, true);
  RRA("Failed to find skin.");
  if (debug) std::cerr << "Found skin, now resolving." << std::endl;

    // get entities adjacent to skin ents from shared_dim down to zero
  for (int this_dim = skin_dim-1; this_dim >= 0; this_dim--) {
    result = mbImpl->get_adjacencies(skin_ents[skin_dim], this_dim,
                                     true, skin_ents[this_dim],
                                     Interface::UNION);
    RRA("Failed getting skin adjacencies.");
  }

    // resolve shared vertices first

    // global id tag
  Tag gid_tag; int def_val = -1;
  if (id_tag)
    gid_tag = *id_tag;
  else {
    result = mbImpl->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int),
                                MB_TAG_DENSE, MB_TYPE_INTEGER, gid_tag,
                                &def_val, true);
    if (MB_FAILURE == result) return result;

    else if (MB_ALREADY_ALLOCATED != result) {
        // just created it, so we need global ids
      result = assign_global_ids(0, skin_dim+1);
      RRA("Failed assigning global ids.");
    }
  }
  
    // store index in temp tag; reuse gid_data 
  gid_data.resize(2*skin_ents[0].size());
  int idx = 0;
  for (Range::iterator rit = skin_ents[0].begin(); 
       rit != skin_ents[0].end(); rit++) 
    gid_data[idx] = idx, idx++;
  Tag idx_tag;
  result = mbImpl->tag_create("__idx_tag", sizeof(int), MB_TAG_DENSE,
                              MB_TYPE_INTEGER, idx_tag, &def_val, true);
  if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) return result;
  result = mbImpl->tag_set_data(idx_tag, skin_ents[0], &gid_data[0]);
  RRA("Couldn't assign index tag.");

    // get gids for skin ents in a vector, to pass to gs
  result = mbImpl->tag_get_data(gid_tag, skin_ents[0], &gid_data[0]);
  RRA("Couldn't get gid tag for skin vertices.");

    // put handles in vector for passing to gs setup
  std::copy(skin_ents[0].begin(), skin_ents[0].end(), 
            std::back_inserter(handle_vec));
  
#ifdef DEBUG_MPE
  MPE_Log_event(SHAREDV_START, procConfig.proc_rank(), "Creating crystal router.");
#endif

    // get a crystal router
  crystal_data *cd = procConfig.crystal_router();

/*  
    // get total number of entities; will overshoot highest global id, but
    // that's ok
  int num_total[2] = {0, 0}, num_local[2] = {0, 0};
  result = mbImpl->get_number_entities_by_dimension(0, 0, num_local);
  if (MB_SUCCESS != result) return result;
  int failure = MPI_Allreduce(num_local, num_total, 1,
                              MPI_INTEGER, MPI_SUM, procConfig.proc_comm());
  if (failure) {
    result = MB_FAILURE;
    RRA("Allreduce for total number of shared ents failed.");
  }
  
*/
    // call gather-scatter to get shared ids & procs
  gs_data *gsd;
  assert(sizeof(ulong_) == sizeof(EntityHandle));
  if (sizeof(int) != sizeof(ulong_)) {
    std::vector<long> lgid_data(gid_data.size());
    std::copy(gid_data.begin(), gid_data.end(), lgid_data.begin());
    gsd = gs_data_setup(skin_ents[0].size(), &lgid_data[0], 
                        (ulong_*)&handle_vec[0], 2, 1, 1, cd);
  }
  else {
    gsd = gs_data_setup(skin_ents[0].size(), (long*)&gid_data[0], 
                        (ulong_*)&handle_vec[0], 2, 1, 1, cd);
  }
  
  if (NULL == gsd) {
    result = MB_FAILURE;
    RRA("Couldn't create gs data.");
  }

    // get shared proc tags
  Tag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
  result = get_shared_proc_tags(sharedp_tag, sharedps_tag, 
                                sharedh_tag, sharedhs_tag, pstatus_tag);
  RRA("Couldn't get shared proc tags.");
  
    // load shared verts into a tuple, then sort by index
  tuple_list shared_verts;
  tuple_list_init_max(&shared_verts, 2, 0, 1, 0, 
                      skin_ents[0].size()*(MAX_SHARING_PROCS+1));
  unsigned int i = 0, j = 0;
  for (unsigned int p = 0; p < gsd->nlinfo->np; p++) 
    for (unsigned int np = 0; np < gsd->nlinfo->nshared[p]; np++) {
      shared_verts.vi[i++] = gsd->nlinfo->sh_ind[j];
      shared_verts.vi[i++] = gsd->nlinfo->target[p];
      shared_verts.vul[j] = gsd->nlinfo->ulabels[j];
      j++;
      shared_verts.n++;
    }
  
  int max_size = skin_ents[0].size()*(MAX_SHARING_PROCS+1);
  buffer sort_buffer;
  buffer_init(&sort_buffer, max_size);
  tuple_list_sort(&shared_verts, 0, &sort_buffer);
  buffer_free(&sort_buffer);

    // set sharing procs and handles tags on skin ents
  int maxp = -1;
  std::vector<int> sharing_procs(MAX_SHARING_PROCS);
  std::fill(sharing_procs.begin(), sharing_procs.end(), maxp);
  j = 0; i = 0;

    // get ents shared by 1 or n procs
  std::map<std::vector<int>, Range> proc_nranges;
  Range proc_verts;
  result = mbImpl->get_adjacencies(proc_ents, 0, false, proc_verts,
                                   Interface::UNION);
  RRA("Couldn't get proc_verts.");
  
  result = tag_shared_verts(shared_verts, skin_ents,
                            proc_nranges, proc_verts);
  RRA("Trouble tagging shared verts.");

#ifdef DEBUG_MPE
  MPE_Log_event(SHAREDV_END, procConfig.proc_rank(), "Finished tag_shared_verts.");
#endif

    // get entities shared by 1 or n procs
  result = tag_shared_ents(resolve_dim, shared_dim, skin_ents,
                           proc_nranges);
  RRA("Trouble tagging shared entities.");

  tuple_list_free(&shared_verts);
  
  if (debug) {
    for (std::map<std::vector<int>, Range>::const_iterator mit = proc_nranges.begin();
         mit != proc_nranges.end(); mit++) {
      std::cout << "Iface: ";
      for (std::vector<int>::const_iterator vit = (mit->first).begin();
           vit != (mit->first).end(); vit++) std::cout << " " << *vit;
      std::cout << std::endl;
    }
  }
  
    // create the sets for each interface; store them as tags on
    // the interface instance
  Range iface_sets;
  result = create_interface_sets(proc_nranges, resolve_dim, shared_dim);
  RRA("Trouble creating iface sets.");

    // establish comm procs and buffers for them
  std::set<unsigned int> procs;
  result = get_interface_procs(procs, true);
  RRA("Trouble getting iface procs.");

#ifndef NDEBUG
  result = check_all_shared_handles(true);
  RRA("Shared handle check failed after iface vertex exchange.");
#endif  

    // resolve shared entity remote handles; implemented in ghost cell exchange
    // code because it's so similar
  result = exchange_ghost_cells(-1, -1, 0, true, true);
  RRA("Trouble resolving shared entity remote handles.");

    // now build parent/child links for interface sets
  result = create_iface_pc_links();
  RRA("Trouble creating interface parent/child links.");

  gs_data_free(gsd);

#ifdef DEBUG_MPE
  MPE_Log_event(RESOLVE_END, procConfig.proc_rank(), "Exiting resolve_shared_ents.");
#endif

    // done
  return result;
}

void ParallelComm::define_mpe() 
{
#ifdef DEBUG_MPE
    // define mpe states used for logging
  int success;
  MPE_Log_get_state_eventIDs( &IFACE_START, &IFACE_END);
  MPE_Log_get_state_eventIDs( &GHOST_START, &GHOST_END);
  MPE_Log_get_state_eventIDs( &SHAREDV_START, &SHAREDV_END);
  MPE_Log_get_state_eventIDs( &RESOLVE_START, &RESOLVE_END);
  MPE_Log_get_state_eventIDs( &ENTITIES_START, &ENTITIES_END);
  MPE_Log_get_state_eventIDs( &RHANDLES_START, &RHANDLES_END);
  success = MPE_Describe_state(IFACE_START, IFACE_END, "Resolve interface ents", "green");
  success = MPE_Describe_state(GHOST_START, GHOST_END, "Exchange ghost ents", "red");
  success = MPE_Describe_state(SHAREDV_START, SHAREDV_END, "Resolve interface vertices", "blue");
  success = MPE_Describe_state(RESOLVE_START, RESOLVE_END, "Resolve shared ents", "purple");
  success = MPE_Describe_state(ENTITIES_START, ENTITIES_END, "Exchange shared ents", "yellow");
  success = MPE_Describe_state(RHANDLES_START, RHANDLES_END, "Remote handles", "cyan");
#endif
}

ErrorCode ParallelComm::resolve_shared_ents(ParallelComm **pc, 
                                                const unsigned int np, 
                                                const int part_dim) 
{
  std::vector<Range> verts(np);
  int tot_verts = 0;
  unsigned int p, i, j, v;
  ErrorCode rval;
  for (p = 0; p < np; p++) {
    Skinner skinner(pc[p]->get_moab());
    Range part_ents, skin_ents;
    rval = pc[p]->get_moab()->get_entities_by_dimension(0, part_dim, part_ents);
    if (MB_SUCCESS != rval) return rval;
    rval = skinner.find_skin(part_ents, false, skin_ents, 0, true);
    if (MB_SUCCESS != rval) return rval;
    rval = pc[p]->get_moab()->get_adjacencies(skin_ents, 0, true, verts[p],
                                              Interface::UNION);
    if (MB_SUCCESS != rval) return rval;
    tot_verts += verts[p].size();
  }
  
  tuple_list shared_ents;
  tuple_list_init_max(&shared_ents, 2, 0, 1, 0, tot_verts);

  i = 0; j = 0;
  std::vector<int> gids;
  Range::iterator rit;
  Tag gid_tag;
  int dum_default = -1;
  for (p = 0; p < np; p++) {
    rval = pc[p]->get_moab()->tag_create(GLOBAL_ID_TAG_NAME, 
                                         sizeof(int), MB_TAG_DENSE,
                                         MB_TYPE_INTEGER, gid_tag, 
                                         &dum_default, true);
    gids.resize(verts[p].size());
    rval = pc[p]->get_moab()->tag_get_data(gid_tag, verts[p], &gids[0]);
    if (MB_SUCCESS != rval) return rval;
    
    for (v = 0, rit = verts[p].begin(); v < gids.size(); v++, rit++) {
      shared_ents.vi[i++] = gids[v];
      shared_ents.vi[i++] = p;
      shared_ents.vul[j] = *rit;
      j++;
      shared_ents.n++;
    }
  }
  
  buffer sort_buffer;
  buffer_init(&sort_buffer, tot_verts);
  tuple_list_sort(&shared_ents, 0, &sort_buffer);
  buffer_free(&sort_buffer);

  j = 0; i = 0;
  std::vector<EntityHandle> handles;
  std::vector<int> procs;
  
  while (i < shared_ents.n) {
    handles.clear();
    procs.clear();
    
      // count & accumulate sharing procs
    int this_gid = shared_ents.vi[j];
    while (i < shared_ents.n && shared_ents.vi[j] == this_gid) {
      j++;
      procs.push_back( shared_ents.vi[j++] );
      handles.push_back( shared_ents.vul[i++] );
    }
    if (1 == procs.size()) continue;
    
    for (v = 0; v < procs.size(); v++) {
      rval = pc[procs[v]]->update_remote_data(handles[v], 
                                              &procs[0], &handles[0], procs.size(),
                                              PSTATUS_INTERFACE);
      if (MB_SUCCESS != rval) return rval;
    }
  }

  std::set<unsigned int> psets;
  for (p = 0; p < np; p++) {
    rval = pc[p]->create_interface_sets(part_dim, part_dim-1);
    if (MB_SUCCESS != rval) return rval;
      // establish comm procs and buffers for them
    psets.clear();
    rval = pc[p]->get_interface_procs(psets, true);
    if (MB_SUCCESS != rval) return rval;
  }

  tuple_list_free(&shared_ents);
  
  return MB_SUCCESS;
}

ErrorCode ParallelComm::tag_iface_entities() 
{
  ErrorCode result = MB_SUCCESS;
  Range iface_ents, tmp_ents, rmv_ents;
  std::vector<unsigned char> pstat;
  unsigned char set_pstat;
  Range::iterator rit2;
  unsigned int i;
  
  for (Range::iterator rit = interfaceSets.begin(); rit != interfaceSets.end(); rit++) {
    iface_ents.clear();
    
    result = mbImpl->get_entities_by_handle(*rit, iface_ents);
    RRA("Couldn't get iface set contents.");
    pstat.resize(iface_ents.size());
    result = mbImpl->tag_get_data(pstatus_tag(), iface_ents, &pstat[0]);
    RRA("Couldn't get pstatus values for set ents.");
    result = mbImpl->tag_get_data(pstatus_tag(), &(*rit), 1, &set_pstat);
    RRA("Couldn't get pstatus values for set.");
    rmv_ents.clear();
    for (rit2 = iface_ents.begin(), i = 0; rit2 != iface_ents.end(); rit2++, i++) {
      if (!(pstat[i] & PSTATUS_INTERFACE)) {
        rmv_ents.insert(*rit2);
        pstat[i] = 0x0;
      }
    }
    result = mbImpl->remove_entities(*rit, rmv_ents);
    RRA("Couldn't remove entities from set.");

    if (!(set_pstat & PSTATUS_NOT_OWNED)) continue;
      // if we're here, we need to set the notowned status on (remaining) set contents

      // remove rmv_ents from the contents list
    iface_ents = subtract(iface_ents, rmv_ents);
      // compress the pstat vector (removing 0x0's)
    std::remove_if(pstat.begin(), pstat.end(), 
                   std::bind2nd(std::equal_to<unsigned char>(), 0x0));
      // fold the not_owned bit into remaining values
    unsigned int sz = iface_ents.size();
    for (i = 0; i < sz; i++)
      pstat[i] |= PSTATUS_NOT_OWNED;

      // set the tag on the entities
    result = mbImpl->tag_set_data(pstatus_tag(), iface_ents, &pstat[0]);
    RRA("Couldn't set pstatus values for set ents.");
  }
  
  return MB_SUCCESS;
}

ErrorCode ParallelComm::set_pstatus_entities(Range &pstatus_ents,
                                                 unsigned char pstatus_val,
                                                 bool lower_dim_ents,
                                                 bool verts_too,
                                                 int operation) 
{
  std::vector<unsigned char> pstatus_vals(pstatus_ents.size());
  Range all_ents, *range_ptr = &pstatus_ents;
  ErrorCode result;
  if (lower_dim_ents || verts_too) {
    all_ents = pstatus_ents;
    range_ptr = &all_ents;
    int start_dim = (lower_dim_ents ? mbImpl->dimension_from_handle(*pstatus_ents.rbegin())-1 : 0);
    for (; start_dim >= 0; start_dim--) {
      result = mbImpl->get_adjacencies(all_ents, start_dim, true, all_ents,
                                       Interface::UNION);
      RRA(" ");
    }
  }
  if (Interface::UNION == operation) {
    result = mbImpl->tag_get_data(pstatus_tag(), *range_ptr, &pstatus_vals[0]);
    RRA("Couldn't get pstatus tag value.");
    for (unsigned int i = 0; i < pstatus_vals.size(); i++)
      pstatus_vals[i] |= pstatus_val;
  }
  else {
    for (unsigned int i = 0; i < pstatus_vals.size(); i++)
      pstatus_vals[i] = pstatus_val;
  }
  result = mbImpl->tag_set_data(pstatus_tag(), *range_ptr, &pstatus_vals[0]);
  RRA("Couldn't set pstatus tag value.");
  
  return MB_SUCCESS;
}
  
ErrorCode ParallelComm::set_pstatus_entities(EntityHandle *pstatus_ents,
                                                 int num_ents,
                                                 unsigned char pstatus_val,
                                                 bool lower_dim_ents,
                                                 bool verts_too,
                                                 int operation) 
{
  std::vector<unsigned char> pstatus_vals(num_ents);
  ErrorCode result;
  if (lower_dim_ents || verts_too) {
      // in this case, call the range-based version
    Range tmp_range;
    std::copy(pstatus_ents, pstatus_ents+num_ents, range_inserter(tmp_range));
    return set_pstatus_entities(tmp_range, pstatus_val, lower_dim_ents, 
                                verts_too, operation);
  }

  if (Interface::UNION == operation) {
    result = mbImpl->tag_get_data(pstatus_tag(), pstatus_ents, num_ents, &pstatus_vals[0]);
    RRA("Couldn't get pstatus tag value.");
    for (unsigned int i = 0; i < (unsigned int) num_ents; i++)
      pstatus_vals[i] |= pstatus_val;
  }
  else {
    for (unsigned int i = 0; i < (unsigned int) num_ents; i++)
      pstatus_vals[i] = pstatus_val;
  }
  result = mbImpl->tag_set_data(pstatus_tag(), pstatus_ents, num_ents, &pstatus_vals[0]);
  RRA("Couldn't set pstatus tag value.");
  
  return MB_SUCCESS;
}
  
ErrorCode ParallelComm::create_interface_sets(int resolve_dim, int shared_dim) 
{
  std::map<std::vector<int>, Range> proc_nranges;
  
    // build up the list of shared entities
  int procs[MAX_SHARING_PROCS];
  EntityHandle handles[MAX_SHARING_PROCS];
  ErrorCode result;
  int nprocs;
  unsigned char pstat;
  for (Range::iterator rit = sharedEnts.begin(); rit != sharedEnts.end(); rit++) {
    if (shared_dim != -1 && mbImpl->dimension_from_handle(*rit) > shared_dim)
      continue;
    result = get_sharing_data(*rit, procs, handles, pstat, nprocs);
    RRA("");
    std::sort(procs, procs+nprocs);
    std::vector<int> tmp_procs(procs, procs + nprocs);
    proc_nranges[tmp_procs].insert(*rit);
  }
                                                  
  Skinner skinner(mbImpl);
  Range skin_ents[4];
  result = mbImpl->get_entities_by_dimension(0, resolve_dim, skin_ents[resolve_dim]);
  RRA("");
  result = skinner.find_skin(skin_ents[resolve_dim], false, 
                             skin_ents[resolve_dim-1], 0, true);
  RRA("Failed to find skin.");
  if (shared_dim > 1) {
    result = mbImpl->get_adjacencies(skin_ents[resolve_dim-1], resolve_dim-2, true,
                                     skin_ents[resolve_dim-2], Interface::UNION);
    RRA("");
  }

  result = tag_shared_ents(resolve_dim, shared_dim, skin_ents,
                           proc_nranges);
    
  return create_interface_sets(proc_nranges, resolve_dim, shared_dim);
}
  
ErrorCode ParallelComm::create_interface_sets(std::map<std::vector<int>, Range> &proc_nranges,
                                                  int resolve_dim, int shared_dim) 
{
  if (proc_nranges.empty()) return MB_SUCCESS;
  
  int proc_ids[MAX_SHARING_PROCS];
  EntityHandle proc_handles[MAX_SHARING_PROCS];
  Tag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
  ErrorCode result = get_shared_proc_tags(sharedp_tag, sharedps_tag, 
                                            sharedh_tag, sharedhs_tag,
                                            pstatus_tag);
  RRA("Trouble getting shared proc tags in create_interface_sets.");
  Range::iterator rit;

    // create interface sets, tag them, and tag their contents with iface set tag
  std::vector<EntityHandle> tag_vals;
  std::vector<unsigned char> pstatus;
  for (std::map<std::vector<int>,Range>::iterator mit = proc_nranges.begin();
       mit != proc_nranges.end(); mit++) {
      // create the set
    EntityHandle new_set;
    result = mbImpl->create_meshset(MESHSET_SET, new_set); 
    RRA("Failed to create interface set.");
    interfaceSets.insert(new_set);

      // add entities
    result = mbImpl->add_entities(new_set, mit->second); 
    RRA("Failed to add entities to interface set.");
      // tag set with the proc rank(s)
    if (mit->first.size() == 1) {
      result = mbImpl->tag_set_data(sharedp_tag, &new_set, 1, 
                                    &(mit->first)[0]); 
      proc_handles[0] = 0;
      result = mbImpl->tag_set_data(sharedh_tag, &new_set, 1, 
                                    proc_handles); 
    }
    else {
      // pad tag data out to MAX_SHARING_PROCS with -1
      assert( mit->first.size() <= MAX_SHARING_PROCS );
      std::copy( mit->first.begin(), mit->first.end(), proc_ids );
      std::fill( proc_ids + mit->first.size(), proc_ids + MAX_SHARING_PROCS, -1 );
      result = mbImpl->tag_set_data(sharedps_tag, &new_set, 1, proc_ids );
      unsigned int ind = std::find(proc_ids, proc_ids+mit->first.size(), procConfig.proc_rank())
          - proc_ids;
      assert(ind < mit->first.size());
      std::fill( proc_handles, proc_handles + MAX_SHARING_PROCS, 0);
      proc_handles[ind] = new_set;
      result = mbImpl->tag_set_data(sharedhs_tag, &new_set, 1, proc_handles); 
    }
    RRA("Failed to tag interface set with procs.");
    
      // get the owning proc, then set the pstatus tag on iface set
    int min_proc = (mit->first)[0];
    unsigned char pval = (PSTATUS_SHARED | PSTATUS_INTERFACE);
    if (min_proc < (int) procConfig.proc_rank()) pval |= PSTATUS_NOT_OWNED;
    if (mit->first.size() > 1) pval |= PSTATUS_MULTISHARED;
    result = mbImpl->tag_set_data(pstatus_tag, &new_set, 1, &pval); 
    RRA("Failed to tag interface set with pstatus.");

      // tag the vertices with the same thing
    pstatus.clear();
    Range verts = (mit->second).subset_by_type(MBVERTEX);
    pstatus.resize(verts.size(), pval);
    result = mbImpl->tag_set_data(pstatus_tag, verts, &pstatus[0]); 
    RRA("Failed to tag interface set vertices with pstatus.");
  }

  return MB_SUCCESS;
}

ErrorCode ParallelComm::create_iface_pc_links() 
{
    // now that we've resolved the entities in the iface sets, 
    // set parent/child links between the iface sets

    // first tag all entities in the iface sets
  Tag tmp_iface_tag;
  EntityHandle tmp_iface_set = 0;
  ErrorCode result = mbImpl->tag_create("__tmp_iface", sizeof(EntityHandle),
                                          MB_TAG_DENSE, MB_TYPE_HANDLE,
                                          tmp_iface_tag, &tmp_iface_set);
  if (MB_ALREADY_ALLOCATED != result && MB_SUCCESS != result) 
    RRA("Failed to create temporary iface set tag.");

  Range iface_ents;
  std::vector<EntityHandle> tag_vals;
  Range::iterator rit;
  
  for (rit = interfaceSets.begin(); rit != interfaceSets.end(); rit++) {
      // tag entities with interface set
    iface_ents.clear();
    result = mbImpl->get_entities_by_handle(*rit, iface_ents);
    RRA("Couldn't get entities in iface set.");
    
    if (iface_ents.empty()) continue;
    
    tag_vals.resize(iface_ents.size());
    std::fill(tag_vals.begin(), tag_vals.end(), *rit);
    result = mbImpl->tag_set_data(tmp_iface_tag, iface_ents, &tag_vals[0]); 
    RRA("Failed to tag iface entities with interface set.");
  }
  
    // now go back through interface sets and add parent/child links
  Range tmp_ents2;
  for (int d = 2; d >= 0; d--) {
    for (rit = interfaceSets.begin(); rit != interfaceSets.end(); rit++) {
        // get entities on this interface
      iface_ents.clear();
      result = mbImpl->get_entities_by_handle(*rit, iface_ents, true);
      RRA("Couldn't get entities by dimension.");
      if (iface_ents.empty() ||
          mbImpl->dimension_from_handle(*iface_ents.rbegin()) != d) continue;

        // get higher-dimensional entities and their interface sets
      result = mbImpl->get_adjacencies(&(*iface_ents.begin()), 1, d+1,
                                       false, tmp_ents2);
      RRA("Couldn't get adjacencies for interface sets.");
      tag_vals.resize(tmp_ents2.size());
      result = mbImpl->tag_get_data(tmp_iface_tag, tmp_ents2, &tag_vals[0]);
      RRA("Couldn't get iface set tag for interface sets.");
      
        // go through and for any on interface make it a parent
      EntityHandle last_set = 0;
      for (unsigned int i = 0; i < tag_vals.size(); i++) {
        if (tag_vals[i] && tag_vals[i] != last_set) {
          result = mbImpl->add_parent_child(tag_vals[i], *rit);
          RRA("Couldn't add parent/child link for interface set.");
          last_set = tag_vals[i];
        }
      }
    }
  }
  
    // delete the temporary tag
  result = mbImpl->tag_delete(tmp_iface_tag);
  RRA("Couldn't delete tmp iface tag.");

  return MB_SUCCESS;
}

ErrorCode ParallelComm::tag_shared_ents(int resolve_dim,
                                            int shared_dim,
                                            Range *skin_ents,
                                            std::map<std::vector<int>, Range> &proc_nranges) 
{
    // set sharing procs tags on other skin ents
  ErrorCode result;
  const EntityHandle *connect; int num_connect;
  std::set<int> sharing_procs;
  std::vector<EntityHandle> dum_connect;
  std::vector<int> sp_vec;

  for (int d = 3; d > 0; d--) {
    if (resolve_dim == d) continue;
    
    for (Range::iterator rit = skin_ents[d].begin();
         rit != skin_ents[d].end(); rit++) {
        // get connectivity
      result = mbImpl->get_connectivity(*rit, connect, num_connect, false,
                                        &dum_connect);
      RRA("Failed to get connectivity on non-vertex skin entities.");
 
      int op = (resolve_dim < shared_dim ? Interface::UNION : Interface::INTERSECT);      
      result = get_sharing_data(connect, num_connect, sharing_procs, op);
      RRA("Failed to get sharing data in tag_shared_ents");
      if (sharing_procs.empty()) continue;

        // intersection is the owning proc(s) for this skin ent
      sp_vec.clear();
      std::copy(sharing_procs.begin(), sharing_procs.end(), std::back_inserter(sp_vec));
      proc_nranges[sp_vec].insert(*rit);
    }
  }

  return MB_SUCCESS;
}

ErrorCode ParallelComm::tag_shared_verts(tuple_list &shared_ents,
                                             Range *skin_ents,
                                             std::map<std::vector<int>, Range> &proc_nranges,
                                             Range &proc_verts) 
{
  Tag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
  ErrorCode result = get_shared_proc_tags(sharedp_tag, sharedps_tag, 
                                            sharedh_tag, sharedhs_tag, pstatus_tag);
  RRA("Trouble getting shared proc tags in tag_shared_verts.");
  
  unsigned int j = 0, i = 0;
  std::vector<int> sharing_procs, sharing_procs2;
  std::vector<EntityHandle> sharing_handles, sharing_handles2;
  
  while (j < 2*shared_ents.n) {
      // count & accumulate sharing procs
    int this_idx = shared_ents.vi[j];
    EntityHandle this_ent = skin_ents[0][this_idx];
    while (j < 2*shared_ents.n && shared_ents.vi[j] == this_idx) {
      j++;
      sharing_procs.push_back( shared_ents.vi[j++] );
      sharing_handles.push_back( shared_ents.vul[i++] );
    }

    if (sharing_procs.size() > 1) {
        // add current proc/handle to list
      sharing_procs.push_back(procConfig.proc_rank());
      sharing_handles.push_back(this_ent);
    }
      
      // sort sharing_procs and sharing_handles such that
      // sharing_procs is in ascending order.  Use temporary
      // lists and binary search to re-order sharing_handles.
    sharing_procs2 = sharing_procs;
    std::sort( sharing_procs2.begin(), sharing_procs2.end() );
    sharing_handles2.resize( sharing_handles.size() );
    for (size_t k = 0; k < sharing_handles.size(); ++k) {
      size_t idx = std::lower_bound( sharing_procs2.begin(), 
                                     sharing_procs2.end(), 
                                     sharing_procs[k] ) - sharing_procs2.begin();
      sharing_handles2[idx] = sharing_handles[k];
    }
    sharing_procs.swap( sharing_procs2 );
    sharing_handles.swap( sharing_handles2 );
    
    
    proc_nranges[sharing_procs].insert(this_ent);

    unsigned char share_flag = PSTATUS_SHARED, 
        ms_flag = (PSTATUS_SHARED | PSTATUS_MULTISHARED);
    if (sharing_procs.size() == 1) {
      result = mbImpl->tag_set_data(sharedp_tag, &this_ent, 1,
                                    &sharing_procs[0]);
      result = mbImpl->tag_set_data(sharedh_tag, &this_ent, 1,
                                    &sharing_handles[0]);
      result = mbImpl->tag_set_data(pstatus_tag, &this_ent, 1, &share_flag);
      RRA("Couldn't set shared tag on shared vertex.");
      sharedEnts.insert(this_ent);
    }
    else {
        // pad lists 
      assert( sharing_procs.size() <= MAX_SHARING_PROCS );
      sharing_procs.resize( MAX_SHARING_PROCS, -1 );
      sharing_handles.resize( MAX_SHARING_PROCS, 0 );
      result = mbImpl->tag_set_data(sharedps_tag, &this_ent, 1,
                                    &sharing_procs[0]);
      result = mbImpl->tag_set_data(sharedhs_tag, &this_ent, 1,
                                    &sharing_handles[0]);
      result = mbImpl->tag_set_data(pstatus_tag, &this_ent, 1, &ms_flag);
      RRA("Couldn't set multi-shared tag on shared vertex.");
      sharedEnts.insert(this_ent);
    }
    RRA("Failed setting shared_procs tag on skin vertices.");

      // reset sharing proc(s) tags
    sharing_procs.clear();
    sharing_handles.clear();
  }

  return MB_SUCCESS;
}
  
  //! get processors with which this processor communicates; sets are sorted by processor
ErrorCode ParallelComm::get_interface_procs(std::set<unsigned int> &procs_set,
                                                bool get_buffs)
{
    // make sure the sharing procs vector is empty
  procs_set.clear();

    // pre-load vector of single-proc tag values
  unsigned int i, j;
  std::vector<int> iface_proc(interfaceSets.size());
  ErrorCode result = mbImpl->tag_get_data(sharedp_tag(), interfaceSets, &iface_proc[0]);
  RRA("Failed to get iface_proc for iface sets.");

    // get sharing procs either from single-proc vector or by getting
    // multi-proc tag value
  int tmp_iface_procs[MAX_SHARING_PROCS];
  std::fill(tmp_iface_procs, tmp_iface_procs+MAX_SHARING_PROCS, -1);
  Range::iterator rit;
  for (rit = interfaceSets.begin(), i = 0; rit != interfaceSets.end(); rit++, i++) {
    if (-1 != iface_proc[i]) procs_set.insert((unsigned int) iface_proc[i]);
    else {
        // get the sharing_procs tag
      result = mbImpl->tag_get_data(sharedps_tag(), &(*rit), 1,
                                    tmp_iface_procs);
      RRA("Failed to get iface_procs for iface set.");
      for (j = 0; j < MAX_SHARING_PROCS; j++) {
        if (-1 != tmp_iface_procs[j] && tmp_iface_procs[j] != (int)procConfig.proc_rank()) 
          procs_set.insert((unsigned int) tmp_iface_procs[j]);
        else if (-1 == tmp_iface_procs[j]) {
          std::fill(tmp_iface_procs, tmp_iface_procs+j, -1);
          break;
        }
      }
    }
  }

  if (get_buffs) {
    for (std::set<unsigned int>::iterator sit = procs_set.begin(); sit != procs_set.end(); sit++)
      get_buffers(*sit);
  }
  
  return MB_SUCCESS;
}
  
ErrorCode ParallelComm::get_pstatus_entities(int dim,
                                                 unsigned char pstatus_val,
                                                 Range &pstatus_ents)
{
  Range ents;
  ErrorCode result;
  
  if (-1 == dim) result = mbImpl->get_entities_by_handle(0, ents);
  else result = mbImpl->get_entities_by_dimension(0, dim, ents);
  RRA(" ");
  
  std::vector<unsigned char> pstatus(ents.size());
  result = mbImpl->tag_get_data(pstatus_tag(), ents, &pstatus[0]);
  RRA("Couldn't get pastatus tag.");
  Range::iterator rit = ents.begin();
  int i = 0;
  if (pstatus_val) {
    for (; rit != ents.end(); i++, rit++)
      if (pstatus[i]&pstatus_val &&
          (-1 == dim || mbImpl->dimension_from_handle(*rit) == dim)) 
        pstatus_ents.insert(*rit);
  }
  else {
    for (; rit != ents.end(); i++, rit++)
      if (!pstatus[i] &&
          (-1 == dim || mbImpl->dimension_from_handle(*rit) == dim)) 
        pstatus_ents.insert(*rit);
  }
  
  return MB_SUCCESS;
}

ErrorCode ParallelComm::check_global_ids(EntityHandle this_set,
                                             const int dimension, 
                                             const int start_id,
                                             const bool largest_dim_only,
                                             const bool parallel)
{
    // global id tag
  Tag gid_tag; int def_val = -1;
  ErrorCode result = mbImpl->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int),
                                          MB_TAG_DENSE, MB_TYPE_INTEGER, gid_tag,
                                          &def_val, true);
  if (MB_ALREADY_ALLOCATED != result &&
      MB_SUCCESS != result) {
    RRA("Failed to create/get gid tag handle.");
  }

  Range dum_range;
  if (MB_ALREADY_ALLOCATED == result) {
    void *tag_ptr = &def_val;
    ErrorCode tmp_result = mbImpl->get_entities_by_type_and_tag(this_set, MBVERTEX, 
                                                                  &gid_tag, &tag_ptr, 1,
                                                                  dum_range);
    if (MB_SUCCESS != tmp_result) {
      result = tmp_result;
      RRA("Failed to get gid tag.");
    }
  }
  
  if (MB_ALREADY_ALLOCATED != result || !dum_range.empty()) {
      // just created it, so we need global ids
    result = assign_global_ids(this_set, dimension, start_id, largest_dim_only,
                               parallel);
    RRA("Failed assigning global ids.");
  }

  return MB_SUCCESS;
}

bool ParallelComm::is_iface_proc(EntityHandle this_set,
                                   int to_proc) 
{
  int sharing_procs[MAX_SHARING_PROCS];
  std::fill(sharing_procs, sharing_procs+MAX_SHARING_PROCS, -1);
  ErrorCode result = mbImpl->tag_get_data(sharedp_tag(), &this_set, 1,
                                            sharing_procs);
  if (to_proc == sharing_procs[0]) return true;
  
  result = mbImpl->tag_get_data(sharedps_tag(), &this_set, 1,
                                sharing_procs);
  for (int i = 0; i < MAX_SHARING_PROCS; i++) {
    if (to_proc == sharing_procs[i]) return true;
    else if (-1 == sharing_procs[i]) return false;
  }
  
  return false;
}

ErrorCode ParallelComm::filter_pstatus( Range &ents,
                                            unsigned char pstat,
                                            unsigned char op,
                                            int to_proc,
                                            Range *returned_ents)
{
  Range tmp_ents;

  //assert(!ents.empty());
  if (ents.empty()) {
    if (returned_ents)
      returned_ents->clear();
    return MB_SUCCESS;
  }

    // Put into tmp_ents any entities which are not owned locally or
    // who are already shared with to_proc
  std::vector<unsigned char> shared_flags(ents.size()), shared_flags2;
  ErrorCode result = mbImpl->tag_get_data(pstatus_tag(), ents,
                                            &shared_flags[0]);
  RRA("Failed to get pstatus flag.");
  Range::const_iterator rit, hint = tmp_ents.begin();;
  int i;
  if (op == PSTATUS_OR) {
    for (rit = ents.begin(), i = 0; rit != ents.end(); rit++, i++) 
      if (((shared_flags[i] & ~pstat)^shared_flags[i]) & pstat) {
        hint = tmp_ents.insert(hint,*rit);
        if (-1 != to_proc) shared_flags2.push_back(shared_flags[i]);
      }
  }
  else if (op == PSTATUS_AND) {
    for (rit = ents.begin(), i = 0; rit != ents.end(); rit++, i++)
      if ((shared_flags[i] & pstat) == pstat) {
        hint = tmp_ents.insert(hint,*rit);
        if (-1 != to_proc) shared_flags2.push_back(shared_flags[i]);
      }
  }
  else if (op == PSTATUS_NOT) {
    for (rit = ents.begin(), i = 0; rit != ents.end(); rit++, i++)
      if (!(shared_flags[i] & pstat)) {
        hint = tmp_ents.insert(hint,*rit);
        if (-1 != to_proc) shared_flags2.push_back(shared_flags[i]);
      }
  }
  else {
    assert(false);
    return MB_FAILURE;
  }

  if (-1 != to_proc) {

    int sharing_procs[MAX_SHARING_PROCS];
    std::fill(sharing_procs, sharing_procs+MAX_SHARING_PROCS, -1);
    Range tmp_ents2;
    hint = tmp_ents2.begin();

    for (rit = tmp_ents.begin(), i = 0; rit != tmp_ents.end(); rit++, i++) {
        // we need to check sharing procs
      if (shared_flags2[i] & PSTATUS_MULTISHARED) {
        result = mbImpl->tag_get_data(sharedps_tag(), &(*rit), 1,
                                      sharing_procs);
        assert(-1 != sharing_procs[0]);
        RRA(" ");
        for (unsigned int j = 0; j < MAX_SHARING_PROCS; j++) {
            // if to_proc shares this entity, add it to list
          if (sharing_procs[j] == to_proc) {
            hint = tmp_ents2.insert(hint, *rit);
          }
          else if (sharing_procs[j] == -1) break;

          sharing_procs[j] = -1;
        }
      }
      else if (shared_flags2[i] & PSTATUS_SHARED) {
        result = mbImpl->tag_get_data(sharedp_tag(), &(*rit), 1,
                                      sharing_procs);
        RRA(" ");
        assert(-1 != sharing_procs[0]);
        if (sharing_procs[0] == to_proc) 
          hint = tmp_ents2.insert(hint,*rit);
        sharing_procs[0] = -1;
      }
      else
        assert("should never get here" && false);
    }

    tmp_ents.swap(tmp_ents2);
  }
  
  if (returned_ents)
    returned_ents->swap(tmp_ents);
  else
    ents.swap(tmp_ents);
  
  return MB_SUCCESS;
}

ErrorCode ParallelComm::exchange_ghost_cells(int ghost_dim, int bridge_dim,
                                                 int num_layers,
                                                 bool store_remote_handles,
                                                 bool wait_all)
{
#ifdef DEBUG_MPE
  if (!num_layers)
    MPE_Log_event(IFACE_START, procConfig.proc_rank(), "Starting interface exchange.");
  else
    MPE_Log_event(GHOST_START, procConfig.proc_rank(), "Starting ghost exchange.");
#endif

#ifdef DEBUG_COMM
//  std::ostringstream pfile("p");
//  pfile << "p" << procConfig.proc_rank() << ".txt";
//  std::cerr.open(pfile.str().c_str(), std::ios_base::trunc);
  std::cerr << "Entering exchange_ghost_cells with num_layers = " 
        << num_layers << std::endl; std::cerr.flush();
#endif
#ifdef DEBUG_MSGS
  msgs.clear();
  msgs.reserve(MAX_SHARING_PROCS);
#endif

    // if we're only finding out about existing ents, we have to be storing
    // remote handles too
  assert(num_layers > 0 || store_remote_handles);
  
  const bool is_iface = !num_layers;

    // get the b-dimensional interface(s) with with_proc, where b = bridge_dim
  
  int success;
  ErrorCode result = MB_SUCCESS;
  int incoming1 = 0, incoming2 = 0;

  reset_all_buffers();
  
    // when this function is called, buffProcs should already have any 
    // communicating procs

    //===========================================
    // post ghost irecv's for ghost entities from all communicating procs
    //===========================================
#ifdef DEBUG_MPE
  MPE_Log_event(ENTITIES_START, procConfig.proc_rank(), "Starting entity exchange.");
#endif
    // index reqs the same as buffer/sharing procs indices
  std::vector<MPI_Request> recv_ent_reqs(2*buffProcs.size(), MPI_REQUEST_NULL),
      recv_remoteh_reqs(2*buffProcs.size(), MPI_REQUEST_NULL);
  std::vector<unsigned int>::iterator proc_it;
  int ind, p;
  sendReqs.resize(2*buffProcs.size(), MPI_REQUEST_NULL);
  for (ind = 0, proc_it = buffProcs.begin(); 
       proc_it != buffProcs.end(); proc_it++, ind++) {
    incoming1++;
    PRINT_DEBUG_IRECV(procConfig.proc_rank(), buffProcs[ind], 
                      remoteOwnedBuffs[ind]->mem_ptr, INITIAL_BUFF_SIZE, 
                      MB_MESG_ENTS_SIZE, incoming1);
    success = MPI_Irecv(remoteOwnedBuffs[ind]->mem_ptr, INITIAL_BUFF_SIZE, 
                        MPI_UNSIGNED_CHAR, buffProcs[ind],
                        MB_MESG_ENTS_SIZE, procConfig.proc_comm(), 
                        &recv_ent_reqs[2*ind]);
    if (success != MPI_SUCCESS) {
      result = MB_FAILURE;
      RRA("Failed to post irecv in ghost exchange.");
    }
  }
  
    //===========================================
    // get entities to be sent to neighbors
    //===========================================

  Range sent_ents[MAX_SHARING_PROCS], allsent, tmp_range;
  std::vector<std::set<unsigned int> > entprocs(allsent.size());
  int dum_ack_buff;
  result = get_sent_ents(is_iface, bridge_dim, ghost_dim, num_layers,
                         sent_ents, allsent, entprocs);
  RRA("get_sent_ents failed.");
  
    //===========================================
    // pack and send ents from this proc to others
    //===========================================
  for (p = 0, proc_it = buffProcs.begin(); 
       proc_it != buffProcs.end(); proc_it++, p++) {

      // reserve space on front for size and for initial buff size
    localOwnedBuffs[p]->reset_buffer(sizeof(int));

      // entities
    result = pack_entities(sent_ents[p], localOwnedBuffs[p], 
                           store_remote_handles, buffProcs[p], is_iface,
                           &entprocs, &allsent); 
    RRA("Packing entities failed.");

#ifdef DEBUG_MSGS
    msgs.resize(msgs.size()+1);
    msgs.back() = new Buffer(*localOwnedBuffs[p]);
      //result = print_buffer(&ownerSBuffs[ind][0], MB_MESG_ENTS_SIZE, *proc_it, true);
#endif

      // send the buffer (size stored in front in send_buffer)
    result = send_buffer(*proc_it, localOwnedBuffs[p], 
                         MB_MESG_ENTS_SIZE, sendReqs[2*p], 
                         recv_ent_reqs[2*p+1], &dum_ack_buff,
                         incoming1,
                         MB_MESG_REMOTEH_SIZE, 
                         (!is_iface && store_remote_handles ? 
                          localOwnedBuffs[p] : NULL),
                         &recv_remoteh_reqs[2*p], &incoming2);
    RRA("Failed to Isend in ghost exchange.");
  }
    

    //===========================================
    // receive/unpack new entities
    //===========================================
    // number of incoming messages for ghosts is the number of procs we 
    // communicate with; for iface, it's the number of those with lower rank
  MPI_Status status;
  std::vector<std::vector<EntityHandle> > recd_ents(buffProcs.size());
  std::vector<std::vector<EntityHandle> > L1hloc(buffProcs.size()), L1hrem(buffProcs.size());
  std::vector<std::vector<int> > L1p(buffProcs.size());
  std::vector<EntityHandle> L2hloc, L2hrem;
  std::vector<unsigned int> L2p;
  Range new_ents;
  
  while (incoming1) {
      // wait for all recvs of ghost ents before proceeding to sending remote handles,
      // b/c some procs may have sent to a 3rd proc ents owned by me;
    PRINT_DEBUG_WAITANY(recv_ent_reqs, MB_MESG_ENTS_SIZE, procConfig.proc_rank());
    
    success = MPI_Waitany(2*buffProcs.size(), &recv_ent_reqs[0], &ind, &status);
    if (MPI_SUCCESS != success) {
      result = MB_FAILURE;
      RRA("Failed in waitany in ghost exchange.");
    }

    PRINT_DEBUG_RECD(status);
    
      // ok, received something; decrement incoming counter
    incoming1--;
    bool done = false;

      // In case ind is for ack, we need index of one before it
    unsigned int base_ind = 2*(ind/2);
    result = recv_buffer(MB_MESG_ENTS_SIZE,
                         status,
                         remoteOwnedBuffs[ind/2],
                         recv_ent_reqs[ind], recv_ent_reqs[ind+1],
                         incoming1,
                         localOwnedBuffs[ind/2], sendReqs[base_ind], sendReqs[base_ind+1],
                         done,
                         (!is_iface && store_remote_handles ? 
                          localOwnedBuffs[ind/2] : NULL),
                         MB_MESG_REMOTEH_SIZE,
                         &recv_remoteh_reqs[base_ind], &incoming2);
    RRA("Failed to receive buffer.");

    if (done) {
#ifdef DEBUG_MSGS
    msgs.resize(msgs.size()+1);
    msgs.back() = new Buffer(*remoteOwnedBuffs[ind/2]);
        //print_buffer(&ghostRBuffs[ind/2][0], MB_MESG_ENTS_SIZE, buffProcs[ind/2], false);
#endif

        // message completely received - process buffer that was sent
    remoteOwnedBuffs[ind/2]->reset_ptr(sizeof(int));
    result = unpack_entities(remoteOwnedBuffs[ind/2],
                             store_remote_handles, ind/2, is_iface,
                             L1hloc, L1hrem, L1p, L2hloc, L2hrem, L2p, new_ents);
      if (MB_SUCCESS != result) {
        std::cout << "Failed to unpack entities.  Buffer contents:" << std::endl;
        print_buffer(remoteOwnedBuffs[ind/2]->mem_ptr, MB_MESG_ENTS_SIZE, buffProcs[ind/2], false);
        return result;
      }

      if (recv_ent_reqs.size() != 2*buffProcs.size()) {
          // post irecv's for remote handles from new proc; shouldn't be iface, 
          // since we know about all procs we share with
        assert(!is_iface);
        recv_remoteh_reqs.resize(2*buffProcs.size(), MPI_REQUEST_NULL);
        for (unsigned int i = recv_ent_reqs.size(); i < 2*buffProcs.size(); i+=2) {
          localOwnedBuffs[i/2]->reset_buffer();
          incoming2++;
          PRINT_DEBUG_IRECV(procConfig.proc_rank(), buffProcs[i/2], 
                            localOwnedBuffs[i/2]->mem_ptr, INITIAL_BUFF_SIZE,
                            MB_MESG_REMOTEH_SIZE, incoming2);
          success = MPI_Irecv(localOwnedBuffs[i/2]->mem_ptr, INITIAL_BUFF_SIZE, 
                              MPI_UNSIGNED_CHAR, buffProcs[i/2],
                              MB_MESG_REMOTEH_SIZE, procConfig.proc_comm(), 
                              &recv_remoteh_reqs[i]);
          if (success != MPI_SUCCESS) {
            result = MB_FAILURE;
            RRA("Failed to post irecv for remote handles in ghost exchange.");
          }
        }
        recv_ent_reqs.resize(2*buffProcs.size(), MPI_REQUEST_NULL);
        sendReqs.resize(2*buffProcs.size(), MPI_REQUEST_NULL);
      }
    }
  }
  
    // add requests for any new addl procs
  if (recv_ent_reqs.size() != 2*buffProcs.size()) {
      // shouldn't get here...
    result = MB_FAILURE;
    RRA("Requests length doesn't match proc count in ghost exchange.");
  }
    
#ifdef DEBUG_MPE
  MPE_Log_event(ENTITIES_END, procConfig.proc_rank(), "Ending entity exchange.");
#endif

  if (is_iface) {
      // need to check over entities I sent and make sure I received 
      // handles for them from all expected procs; if not, need to clean
      // them up
    result = check_clean_iface(allsent);
    if (MB_SUCCESS != result) std::cout << "Failed check." << std::endl;
    
      // now set the shared/interface tag on non-vertex entities on interface
    result = tag_iface_entities();
    RRA("Failed to tag iface entities.");

#ifndef NDEBUG
    result = check_sent_ents(allsent);
    if (MB_SUCCESS != result) std::cout << "Failed check." << std::endl;
    result = check_all_shared_handles(true);
    if (MB_SUCCESS != result) std::cout << "Failed check." << std::endl;
#endif

#ifdef DEBUG_MPE
    MPE_Log_event(IFACE_END, procConfig.proc_rank(), "Ending interface exchange.");
#endif
#ifdef DEBUG_COMM
    std::cerr << "Exiting exchange_ghost_cells" << std::endl; std::cerr.flush();
#endif

      //===========================================
      // wait if requested
      //===========================================
    if (wait_all) {
#ifdef DEBUG_BARRIER
      success = MPI_Barrier(procConfig.proc_comm());
#else
      success = MPI_Waitall(2*buffProcs.size(), &recv_ent_reqs[0], &status);
      success = MPI_Waitall(2*buffProcs.size(), &sendReqs[0], &status);
#endif
      if (MPI_SUCCESS != success) {
        result = MB_FAILURE;
        RRA("Failed in waitall in ghost exchange.");
      }
    }
    return MB_SUCCESS;
  }

    //===========================================
    // send local handles for new ghosts to owner, then add
    // those to ghost list for that owner
    //===========================================
  for (p = 0, proc_it = buffProcs.begin(); 
       proc_it != buffProcs.end(); proc_it++, p++) {

      // reserve space on front for size and for initial buff size
    remoteOwnedBuffs[p]->reset_buffer(sizeof(int));

    result = pack_remote_handles(L1hloc[p], L1hrem[p], L1p[p], *proc_it,
                                 remoteOwnedBuffs[p]);
    RRA("Failed to pack remote handles.");
    remoteOwnedBuffs[p]->set_stored_size();
    
#ifdef DEBUG_MSGS
    msgs.resize(msgs.size()+1);
    msgs.back() = new Buffer(*remoteOwnedBuffs[p]);
      //print_buffer(&ownerSBuffs[ind][0], MB_MESG_REMOTEH_SIZE, buffProcs[ind], true);
#endif  
    result = send_buffer(buffProcs[p], remoteOwnedBuffs[p], 
                         MB_MESG_REMOTEH_SIZE, 
                         sendReqs[2*p], recv_remoteh_reqs[2*p+1], 
                         &dum_ack_buff, incoming2);
    RRA("Failed to send remote handles.");
  }
  
    //===========================================
    // process remote handles of my ghosteds
    //===========================================
  while (incoming2) {
    PRINT_DEBUG_WAITANY(recv_remoteh_reqs, MB_MESG_REMOTEH_SIZE, procConfig.proc_rank());
    success = MPI_Waitany(2*buffProcs.size(), &recv_remoteh_reqs[0], &ind, &status);
    if (MPI_SUCCESS != success) {
      result = MB_FAILURE;
      RRA("Failed in waitany in ghost exchange.");
    }
    
      // ok, received something; decrement incoming counter
    incoming2--;

    PRINT_DEBUG_RECD(status);
    
    bool done = false;
    unsigned int base_ind = 2*(ind/2);
    result = recv_buffer(MB_MESG_REMOTEH_SIZE, status, 
                         localOwnedBuffs[ind/2], 
                         recv_remoteh_reqs[ind], recv_remoteh_reqs[ind+1], incoming2,
                         remoteOwnedBuffs[ind/2], 
                         sendReqs[base_ind], sendReqs[base_ind+1],
                         done);
    RRA("Failed to receive remote handles.");
    if (done) {
        // incoming remote handles
#ifdef DEBUG_MSGS
    msgs.resize(msgs.size()+1);
    msgs.back() = new Buffer(*localOwnedBuffs[ind]);
      //print_buffer(&remotehRBuffs[ind/2][0], MB_MESG_REMOTEH_SIZE, buffProcs[ind/2], false);
#endif  
    localOwnedBuffs[ind/2]->reset_ptr(sizeof(int));
    result = unpack_remote_handles(buffProcs[ind/2], 
                                   localOwnedBuffs[ind/2],
                                   L2hloc, L2hrem, L2p);
      RRA("Failed to unpack remote handles.");
    }
  }
    
#ifdef DEBUG_MPE
  MPE_Log_event(RHANDLES_END, procConfig.proc_rank(), "Ending remote handles.");
#endif
#ifdef DEBUG_MPE
  MPE_Log_event(GHOST_END, procConfig.proc_rank(), 
                "Ending ghost exchange (still doing checks).");
#endif

    //===========================================
    // wait if requested
    //===========================================
  if (wait_all) {
#ifdef DEBUG_BARRIER
    success = MPI_Barrier(procConfig.proc_comm());
#else
    success = MPI_Waitall(2*buffProcs.size(), &recv_remoteh_reqs[0], &status);
    success = MPI_Waitall(2*buffProcs.size(), &sendReqs[0], &status);
#endif
    if (MPI_SUCCESS != success) {
      result = MB_FAILURE;
      RRA("Failed in waitall in ghost exchange.");
    }
  }

#ifndef NDEBUG
  result = check_sent_ents(allsent);
  RRA("Failed check on shared entities.");
  result = check_all_shared_handles(true);
  RRA("Failed check on all shared handles.");
#endif
#ifdef DEBUG_COMM
  std::cerr << "Exiting exchange_ghost_cells" << std::endl; std::cerr.flush();
#endif

  return MB_SUCCESS;
}

ErrorCode ParallelComm::send_buffer(const unsigned int to_proc,
                                        Buffer *send_buff,
                                        int mesg_tag,
                                        MPI_Request &send_req,
                                        MPI_Request &ack_req,
                                        int *ack_buff,
                                        int &this_incoming,
                                        int next_mesg_tag,
                                        Buffer *next_recv_buff,
                                        MPI_Request *next_recv_req,
                                        int *next_incoming) 
{
  ErrorCode result = MB_SUCCESS;
  int success;

    // if small message, post recv for remote handle message
  if (send_buff->get_stored_size() <= (int)INITIAL_BUFF_SIZE && next_recv_buff) {
    (*next_incoming)++;
    PRINT_DEBUG_IRECV(procConfig.proc_rank(), to_proc, next_recv_buff->mem_ptr,
                      INITIAL_BUFF_SIZE, next_mesg_tag, *next_incoming);
    success = MPI_Irecv(next_recv_buff->mem_ptr, INITIAL_BUFF_SIZE, 
                        MPI_UNSIGNED_CHAR, to_proc,
                        next_mesg_tag, procConfig.proc_comm(), 
                        next_recv_req);
    if (success != MPI_SUCCESS) {
      result = MB_FAILURE;
      RRA("Failed to post irecv for next message in ghost exchange.");
    }
  }
    // if large, we'll need an ack before sending the rest
  else if (send_buff->get_stored_size() > (int)INITIAL_BUFF_SIZE) {
    this_incoming++;
    PRINT_DEBUG_IRECV(procConfig.proc_rank(), to_proc, (unsigned char*)ack_buff,
                      sizeof(int), mesg_tag-1, this_incoming);
    success = MPI_Irecv(ack_buff, sizeof(int), 
                        MPI_UNSIGNED_CHAR, to_proc,
                        mesg_tag-1, procConfig.proc_comm(), 
                        &ack_req);
    if (success != MPI_SUCCESS) {
      result = MB_FAILURE;
      RRA("Failed to post irecv for entity ack in ghost exchange.");
    }
  }

    // send the buffer
  PRINT_DEBUG_ISEND(procConfig.proc_rank(), to_proc, send_buff->mem_ptr, mesg_tag,
                    std::min(send_buff->get_stored_size(), (int)INITIAL_BUFF_SIZE));
  assert(0 <= send_buff->get_stored_size() && 
         send_buff->get_stored_size() <= (int)send_buff->alloc_size);
  success = MPI_Isend(send_buff->mem_ptr, 
                      std::min(send_buff->get_stored_size(), 
                               (int)INITIAL_BUFF_SIZE),
                      MPI_UNSIGNED_CHAR, to_proc, 
                      mesg_tag, procConfig.proc_comm(), &send_req);
  if (success != MPI_SUCCESS) return MB_FAILURE;

  return result;
}

ErrorCode ParallelComm::recv_buffer(int mesg_tag_expected,
                                        const MPI_Status &mpi_status,
                                        Buffer *recv_buff,
                                        MPI_Request &recv_req,
                                        MPI_Request &ack_recvd_req,
                                        int &this_incoming,
                                        Buffer *send_buff,
                                        MPI_Request &send_req,
                                        MPI_Request &sent_ack_req,
                                        bool &done,
                                        Buffer *next_buff,
                                        int next_tag,
                                        MPI_Request *next_req,
                                        int *next_incoming) 
{
    // process a received message; if there will be more coming, 
    // post a receive for 2nd part then send an ack message
    //
  int from_proc = mpi_status.MPI_SOURCE;
  int success;
  ErrorCode result = MB_SUCCESS;

    // set the buff_ptr on the recv_buffer; needs to point beyond any
    // valid data already in the buffer
  recv_buff->reset_ptr(std::min(recv_buff->get_stored_size(), 
                                (int)recv_buff->alloc_size));
  
  if (mpi_status.MPI_TAG == mesg_tag_expected &&
      recv_buff->get_stored_size() > (int)INITIAL_BUFF_SIZE) {
      // 1st message & large - allocate buffer, post irecv for 2nd message,
      // then send ack
    recv_buff->reserve(recv_buff->get_stored_size());
    assert(recv_buff->alloc_size > INITIAL_BUFF_SIZE);

      // will expect a 2nd message
    this_incoming++;

    PRINT_DEBUG_IRECV(procConfig.proc_rank(), from_proc, 
                      recv_buff->mem_ptr+INITIAL_BUFF_SIZE,
                      recv_buff->get_stored_size() - INITIAL_BUFF_SIZE,
                      mesg_tag_expected+1, this_incoming);
    success = MPI_Irecv(recv_buff->mem_ptr+INITIAL_BUFF_SIZE, 
                        recv_buff->get_stored_size() - INITIAL_BUFF_SIZE, 
                        MPI_UNSIGNED_CHAR, from_proc,
                        mesg_tag_expected+1, procConfig.proc_comm(), 
                        &recv_req);
    if (success != MPI_SUCCESS) {
      result = MB_FAILURE;
      RRA("Failed to post 2nd iRecv in ghost exchange.");
    }

      // send ack, doesn't matter what data actually is
    PRINT_DEBUG_ISEND(procConfig.proc_rank(), from_proc, recv_buff->mem_ptr, 
                      mesg_tag_expected-1, sizeof(int));
    success = MPI_Isend(recv_buff->mem_ptr, sizeof(int),
                        MPI_UNSIGNED_CHAR, from_proc, 
                        mesg_tag_expected-1, procConfig.proc_comm(), &sent_ack_req);
    if (success != MPI_SUCCESS) {
      result = MB_FAILURE;
      RRA("Failed to send ack in ghost exchange.");
    }
  }

  else if (mpi_status.MPI_TAG == mesg_tag_expected-1) {
      // got an ack back, send the 2nd half of message

      // should be a large message if we got this
    assert(*((size_t*)send_buff->mem_ptr) > INITIAL_BUFF_SIZE);

      // post irecv for next message, then send 2nd message
    if (next_buff) {
        // we'll expect a return message
      (*next_incoming)++;
      PRINT_DEBUG_IRECV(procConfig.proc_rank(), from_proc, next_buff->mem_ptr,
                        INITIAL_BUFF_SIZE, next_tag, *next_incoming);

      success = MPI_Irecv(next_buff->mem_ptr, 
                          INITIAL_BUFF_SIZE, 
                          MPI_UNSIGNED_CHAR, from_proc,
                          next_tag, procConfig.proc_comm(), 
                          next_req);
      if (success != MPI_SUCCESS) {
        result = MB_FAILURE;
        RRA("Failed to post next irecv in ghost exchange.");
      }

    }
    
      // send 2nd message
    PRINT_DEBUG_ISEND(procConfig.proc_rank(), from_proc, 
                      send_buff->mem_ptr+INITIAL_BUFF_SIZE,
                      mesg_tag_expected+1,
                      send_buff->get_stored_size() - INITIAL_BUFF_SIZE);
    
    assert(send_buff->get_stored_size()-INITIAL_BUFF_SIZE < send_buff->alloc_size &&
           0 <= send_buff->get_stored_size());
    success = MPI_Isend(send_buff->mem_ptr+INITIAL_BUFF_SIZE, 
                        send_buff->get_stored_size() - INITIAL_BUFF_SIZE,
                        MPI_UNSIGNED_CHAR, from_proc, mesg_tag_expected+1, 
                        procConfig.proc_comm(), &send_req);
    if (success != MPI_SUCCESS) {
      result = MB_FAILURE;
      RRA("Failed to send 2nd message in ghost exchange.");
    }
  }
  else if ((mpi_status.MPI_TAG == mesg_tag_expected && 
            recv_buff->get_stored_size() <= (int)INITIAL_BUFF_SIZE) ||
           mpi_status.MPI_TAG == mesg_tag_expected+1) {
      // message completely received - signal that we're done
    done = true;
  }
  
  return MB_SUCCESS;
}

ErrorCode ParallelComm::check_clean_iface(Range &allsent) 
{
    // allsent is all entities I think are on interface; go over them, looking
    // for zero-valued handles, and fix any I find

  ErrorCode result = MB_SUCCESS;
  Range::iterator rit;
  unsigned char pstatus;
  int sharedp[MAX_SHARING_PROCS], nump;
  EntityHandle sharedh[MAX_SHARING_PROCS];
  for (rit = allsent.begin(); rit != allsent.end(); rit++) {
    result = get_sharing_data(*rit, sharedp, sharedh, pstatus, nump);
    RRA("");
    assert("Should be shared with at least one other proc" && 
           (nump > 1 || sharedp[0] != (int)procConfig.proc_rank()));
    int numz = 0;
    for (int i = 0; i < nump; i++) {
      if (!sharedh[i]) numz++;
      else if (numz) {
          // shift downward
        sharedh[i-numz] = sharedh[i];
        sharedp[i-numz] = sharedp[i];
      }
    }
    if (numz) {
      for (int i = numz; i > 0; i--) {
        sharedp[nump-i] = -1;
        sharedh[nump-i] = 0;
      }
      result = set_sharing_data(*rit, pstatus, nump, nump-numz, sharedp, sharedh);
      RRA("");
    }
  }
  
  return result;
}

ErrorCode ParallelComm::set_sharing_data(EntityHandle ent, unsigned char pstatus,
                                             int old_nump, int new_nump,
                                             int *ps, EntityHandle *hs) 
{
  assert("Should call this function only when reducing sharing procs." &&
         new_nump < old_nump);

    // set sharing data to what's passed in; may have to clean up existing sharing tags
    // if things changed too much
  
  ErrorCode result;
  if (pstatus & PSTATUS_MULTISHARED && new_nump < 3) {
      // need to remove multishared tags
    result = mbImpl->tag_delete_data(sharedps_tag(), &ent, 1);
    RRA("");
    result = mbImpl->tag_delete_data(sharedhs_tag(), &ent, 1);
    RRA("");
    pstatus ^= PSTATUS_MULTISHARED;
    if (new_nump < 2) pstatus = 0x0;
  }
  else if (pstatus & PSTATUS_SHARED && new_nump < 2) {
    hs[0] = 0;
    ps[0] = -1;
    pstatus = 0x0;
  }

  if (new_nump > 2) {
    result = mbImpl->tag_set_data(sharedps_tag(), &ent, 1, ps);
    RRA("");
    result = mbImpl->tag_set_data(sharedhs_tag(), &ent, 1, hs);
    RRA("");
#ifndef NDEBUG
    {
        // check for duplicates in proc list
      std::set<unsigned int> dumprocs;
      int dp = 0;
      for (; dp < old_nump && -1 != ps[dp]; dp++)
        dumprocs.insert(ps[dp]);
      assert(dp == (int)dumprocs.size());
    }
#endif      
  }
  else {
    unsigned int j = (ps[0] == (int)procConfig.proc_rank() ? 1 : 0);
    assert(-1 != ps[j]);
    result = mbImpl->tag_set_data(sharedp_tag(), &ent, 1, ps+j);
    RRA("");
    result = mbImpl->tag_set_data(sharedh_tag(), &ent, 1, hs+j);
    RRA("");
  }
  
  result = mbImpl->tag_set_data(pstatus_tag(), &ent, 1, &pstatus);
  RRA("");

  if (old_nump > 1 && new_nump < 2) sharedEnts.erase(ent);

  return result;
}

ErrorCode ParallelComm::get_sent_ents(const bool is_iface, 
                                          const int bridge_dim, const int ghost_dim,
                                          const int num_layers,
                                          Range *sent_ents, Range &allsent,
                                          std::vector<std::set<unsigned int> > &entprocs) 
{
  ErrorCode result;
  int ind;
  std::vector<unsigned int>::iterator proc_it;
  Range tmp_range;
  
    // done in a separate loop over procs because sometimes later procs 
    // need to add info to earlier procs' messages
  for (ind = 0, proc_it = buffProcs.begin(); 
       proc_it != buffProcs.end(); proc_it++, ind++) {
    if (!is_iface) {
      result = get_ghosted_entities(bridge_dim, ghost_dim, buffProcs[ind],
                                    num_layers, sent_ents[ind]);
      RRA("Failed to get ghost layers.");
    }
    else {
      result = get_iface_entities(buffProcs[ind], -1, sent_ents[ind]);
      RRA("Failed to get interface layers.");
    }

      // filter out entities already shared with destination
    tmp_range.clear();
    result = filter_pstatus(sent_ents[ind], PSTATUS_SHARED, PSTATUS_AND,
                            buffProcs[ind], &tmp_range);
    RRA("Couldn't filter on owner.");
    if (!tmp_range.empty()) 
      sent_ents[ind] = subtract( sent_ents[ind], tmp_range);

    allsent.merge(sent_ents[ind]);
  }

    //===========================================
    // need to get procs each entity is sent to
    //===========================================
  Range::iterator rit;
  entprocs.resize(allsent.size());
  for (ind = 0, proc_it = buffProcs.begin(); 
       proc_it != buffProcs.end(); proc_it++, ind++) {
    for (rit = sent_ents[ind].begin(); rit != sent_ents[ind].end(); rit++) {
      int rind = allsent.index(*rit);
      assert(rind < (int) allsent.size() && rind >= 0);
      entprocs[rind].insert(*proc_it);
    }
  }

  return MB_SUCCESS;
}

ErrorCode ParallelComm::exchange_ghost_cells(ParallelComm **pcs,
                                                 unsigned int num_procs,
                                                 int ghost_dim, int bridge_dim,
                                                 int num_layers,
                                                 bool store_remote_handles)
{
    // static version of function, exchanging info through buffers rather 
    // than through messages

    // if we're only finding out about existing ents, we have to be storing
    // remote handles too
  assert(num_layers > 0 || store_remote_handles);
  
  const bool is_iface = !num_layers;
  
  unsigned int ind;
  ParallelComm *pc;
  ErrorCode result = MB_SUCCESS;

    // when this function is called, buffProcs should already have any 
    // communicating procs

    //===========================================
    // get entities to be sent to neighbors
    //===========================================

    // done in a separate loop over procs because sometimes later procs 
    // need to add info to earlier procs' messages
  Range sent_ents[MAX_SHARING_PROCS][MAX_SHARING_PROCS], 
      allsent[MAX_SHARING_PROCS];

    //===========================================
    // get entities to be sent to neighbors
    //===========================================

  std::vector<std::set<unsigned int> > entprocs[MAX_SHARING_PROCS];
  for (unsigned int p = 0; p < num_procs; p++) {
    pc = pcs[p];
    result = pc->get_sent_ents(is_iface, bridge_dim, ghost_dim, num_layers,
                               sent_ents[p], allsent[p], entprocs[p]);
    RRAI(pc->get_moab(), "get_sent_ents failed.");
  
    //===========================================
    // pack entities into buffers
    //===========================================

    for (ind = 0; ind < pc->buffProcs.size(); ind++) {
        // entities
      pc->localOwnedBuffs[ind]->reset_ptr(sizeof(int));
      result = pc->pack_entities(sent_ents[p][ind], pc->localOwnedBuffs[ind],
                                 store_remote_handles, pc->buffProcs[ind], is_iface,
                                 &entprocs[p], &allsent[p]); 
      RRAI(pc->get_moab(), "Packing entities failed.");
    }
  }

    //===========================================
    // receive/unpack new entities
    //===========================================
    // number of incoming messages for ghosts is the number of procs we 
    // communicate with; for iface, it's the number of those with lower rank
  std::vector<std::vector<EntityHandle> > L1hloc[MAX_SHARING_PROCS], L1hrem[MAX_SHARING_PROCS];
  std::vector<std::vector<int> > L1p[MAX_SHARING_PROCS];
  std::vector<EntityHandle> L2hloc[MAX_SHARING_PROCS], L2hrem[MAX_SHARING_PROCS];
  std::vector<unsigned int> L2p[MAX_SHARING_PROCS];
  Range new_ents[MAX_SHARING_PROCS];
  
  for (unsigned int p = 0; p < num_procs; p++) {
    L1hloc[p].resize(pcs[p]->buffProcs.size());
    L1hrem[p].resize(pcs[p]->buffProcs.size());
    L1p[p].resize(pcs[p]->buffProcs.size());
  }
  
  for (unsigned int p = 0; p < num_procs; p++) {
  
    ParallelComm *pc = pcs[p];
    
    for (ind = 0; ind < pc->buffProcs.size(); ind++) {
        // incoming ghost entities; unpack; returns entities received
        // both from sending proc and from owning proc (which may be different)

        // buffer could be empty, which means there isn't any message to
        // unpack (due to this comm proc getting added as a result of indirect
        // communication); just skip this unpack
      if (pc->localOwnedBuffs[ind]->get_stored_size() == 0) continue;

      unsigned int to_p = pc->buffProcs[ind];
      pc->localOwnedBuffs[ind]->reset_ptr(sizeof(int));
      result = pcs[to_p]->unpack_entities(pc->localOwnedBuffs[ind],
                                          store_remote_handles, ind, is_iface,
                                          L1hloc[to_p], L1hrem[to_p], L1p[to_p], L2hloc[to_p], 
                                          L2hrem[to_p], L2p[to_p], new_ents[to_p]);
      RRAI(pc->get_moab(), "Failed to unpack entities.");
    }
  }

  if (is_iface) {
      // need to check over entities I sent and make sure I received 
      // handles for them from all expected procs; if not, need to clean
      // them up
    for (unsigned int p = 0; p < num_procs; p++) {
      result = pcs[p]->check_clean_iface(allsent[p]);
      RRAI(pcs[p]->get_moab(), "Failed check on shared entities.");
    }

#ifndef NDEBUG
    for (unsigned int p = 0; p < num_procs; p++) {
      result = pcs[p]->check_sent_ents(allsent[p]);
      RRAI(pcs[p]->get_moab(), "Failed check on shared entities.");
    }
    result = check_all_shared_handles(pcs, num_procs);
    RRAI(pcs[0]->get_moab(), "Failed check on all shared handles.");
#endif
    return MB_SUCCESS;
  }
  
      //===========================================
      // send local handles for new ghosts to owner, then add
      // those to ghost list for that owner
      //===========================================
  std::vector<unsigned int>::iterator proc_it;
  for (unsigned int p = 0; p < num_procs; p++) {
    pc = pcs[p];
  
    for (ind = 0, proc_it = pc->buffProcs.begin(); 
         proc_it != pc->buffProcs.end(); proc_it++, ind++) {
        // skip if iface layer and higher-rank proc
      pc->localOwnedBuffs[ind]->reset_ptr(sizeof(int));
      result = pc->pack_remote_handles(L1hloc[p][ind], L1hrem[p][ind], L1p[p][ind], *proc_it,
                                       pc->localOwnedBuffs[ind]);
      RRAI(pc->get_moab(), "Failed to pack remote handles.");
    }
  }
  
    //===========================================
    // process remote handles of my ghosteds
    //===========================================
  for (unsigned int p = 0; p < num_procs; p++) {
    pc = pcs[p];
  
    for (ind = 0, proc_it = pc->buffProcs.begin(); 
         proc_it != pc->buffProcs.end(); proc_it++, ind++) {
        // incoming remote handles
      unsigned int to_p = pc->buffProcs[ind];
      pc->localOwnedBuffs[ind]->reset_ptr(sizeof(int));
      result = pcs[to_p]->unpack_remote_handles(p, 
                                                pc->localOwnedBuffs[ind],
                                                L2hloc[to_p], L2hrem[to_p], L2p[to_p]);
      RRAI(pc->get_moab(), "Failed to unpack remote handles.");
    }
  }
    
#ifndef NDEBUG
  for (unsigned int p = 0; p < num_procs; p++) {
    result = pcs[p]->check_sent_ents(allsent[p]);
    RRAI(pcs[p]->get_moab(), "Failed check on shared entities.");
  }
  
  result = ParallelComm::check_all_shared_handles(pcs, num_procs);
  RRAI(pcs[0]->get_moab(), "Failed check on all shared handles.");
#endif

  return MB_SUCCESS;
}

ErrorCode ParallelComm::get_iface_entities(int other_proc,
                                               int dim,
                                               Range &iface_ents) 
{
  Range iface_sets;
  ErrorCode result = MB_SUCCESS;
  
  for (Range::iterator rit = interfaceSets.begin(); rit != interfaceSets.end(); rit++) {
    if (-1 != other_proc && !is_iface_proc(*rit, other_proc)) continue;
    
    if (-1 == dim) result = mbImpl->get_entities_by_handle(*rit, iface_ents);
    else result = mbImpl->get_entities_by_dimension(*rit, dim, iface_ents);
    RRA(" Failed to get entities in iface set.");
  }
  
  return MB_SUCCESS;
}

ErrorCode ParallelComm::check_sent_ents(Range &allsent) 
{
    // check entities to make sure there are no zero-valued remote handles
    // where they shouldn't be
  std::vector<unsigned char> pstat(allsent.size());
  ErrorCode result = mbImpl->tag_get_data(pstatus_tag(), allsent, &pstat[0]);
  RRA("Trouble getting pstatus.");
  std::vector<EntityHandle> handles(allsent.size());
  result = mbImpl->tag_get_data(sharedh_tag(), allsent, &handles[0]);
  RRA("Trouble getting shared handles.");
  std::vector<int> procs(allsent.size());
  result = mbImpl->tag_get_data(sharedp_tag(), allsent, &procs[0]);
  RRA("Trouble getting shared procs.");

  Range bad_entities;
  
  Range::iterator rit;
  unsigned int i;
  EntityHandle dum_hs[MAX_SHARING_PROCS];
  int dum_ps[MAX_SHARING_PROCS];
  
  for (rit = allsent.begin(), i = 0; rit != allsent.end(); rit++, i++) {
    if (-1 != procs[i] && 0 == handles[i]) bad_entities.insert(*rit);
    else {
        // might be multi-shared...
      result = mbImpl->tag_get_data(sharedps_tag(), &(*rit), 1, dum_ps);
      if (MB_TAG_NOT_FOUND == result) continue;
      RRA("Trouble getting sharedps.");
      result = mbImpl->tag_get_data(sharedhs_tag(), &(*rit), 1, dum_hs);
      RRA("Trouble getting sharedhs.");

        // find first non-set proc
      int *ns_proc = std::find(dum_ps, dum_ps+MAX_SHARING_PROCS, -1);
      int num_procs = ns_proc-dum_ps;
      assert(num_procs <= MAX_SHARING_PROCS);
        // now look for zero handles in active part of dum_hs
      EntityHandle *ns_handle = std::find(dum_hs, dum_hs+num_procs, 0);
      int num_handles = ns_handle-dum_hs;
      assert(num_handles <= num_procs);
      if (num_handles != num_procs) bad_entities.insert(*rit);
    }
  }
  
  return MB_SUCCESS;
}

ErrorCode ParallelComm::pack_remote_handles(std::vector<EntityHandle> &L1hloc,
                                                std::vector<EntityHandle> &L1hrem,
                                                std::vector<int> &L1p,
                                                unsigned int to_proc,
                                                Buffer *buff) 
{
    // 2 vectors of handles plus ints
  buff->check_space(((L1p.size()+1)*sizeof(int) + 
                     (L1hloc.size()+1)*sizeof(EntityHandle) + 
                     (L1hrem.size()+1)*sizeof(EntityHandle)));
  
    // should be in pairs of handles
  PACK_INT(buff->buff_ptr, L1hloc.size());
  PACK_INTS(buff->buff_ptr, &L1p[0], L1p.size());
  PACK_EH(buff->buff_ptr, &L1hrem[0], L1hrem.size());
  PACK_EH(buff->buff_ptr, &L1hloc[0], L1hloc.size());
  
  buff->set_stored_size();
  
  return MB_SUCCESS;
}

ErrorCode ParallelComm::unpack_remote_handles(unsigned int from_proc,
                                                  Buffer *buff,
                                                  std::vector<EntityHandle> &L2hloc,
                                                  std::vector<EntityHandle> &L2hrem,
                                                  std::vector<unsigned int> &L2p)
{
    // incoming remote handles; use to set remote handles
  int num_eh;
  UNPACK_INT(buff->buff_ptr, num_eh);

  unsigned long extra;

  unsigned char *buff_proc = buff->buff_ptr;

    // buff->buff_ptr is already int-aligned, due to last unpack_int
  buff->buff_ptr += num_eh * sizeof(int);
  ALIGN_BUFFER(buff->buff_ptr, extra);
    // buff_rem is already int-aligned, due to last align_buffer
  unsigned char *buff_rem = buff->buff_ptr + num_eh * sizeof(EntityHandle);
  ErrorCode result;
  EntityHandle hpair[2], dum_h;
  int proc;
  for (int i = 0; i < num_eh; i++) {
    UNPACK_INT(buff_proc, proc);
    UNPACK_EH(buff->buff_ptr, hpair, 1);
    UNPACK_EH(buff_rem, hpair+1, 1);

    if (-1 != proc) {
      result = find_existing_entity(false, proc, hpair[0], 3, NULL, 0,
                                    mbImpl->type_from_handle(hpair[1]),
                                    L2hloc, L2hrem, L2p, dum_h);
      RRA("Didn't get existing entity.");
      if (dum_h) hpair[0] = dum_h;
      else hpair[0] = 0;
    }
    if (!(hpair[0] && hpair[1])) return MB_FAILURE;
    int this_proc = from_proc;
    result = update_remote_data(hpair[0], &this_proc, hpair+1, 1, 0);
    RRA("Trouble setting remote data range on sent entities in ghost exchange.");
  }
  
  return MB_SUCCESS;
}

ErrorCode ParallelComm::get_ghosted_entities(int bridge_dim,
                                                 int ghost_dim,
                                                 int to_proc, 
                                                 int num_layers,
                                                 Range &ghosted_ents) 
{
    // get bridge ents on interface(s)
  Range from_ents;
  ErrorCode result = MB_SUCCESS;
  assert(0 < num_layers);
  for (Range::iterator rit = interfaceSets.begin(); rit != interfaceSets.end();
       rit++) {
    if (!is_iface_proc(*rit, to_proc)) continue;
      
      // get starting "from" entities
    if (bridge_dim == -1)
      result = mbImpl->get_entities_by_handle(*rit, from_ents);
    else
      result = mbImpl->get_entities_by_dimension(*rit, bridge_dim, from_ents);
    RRA("Couldn't get bridge ents in the set.");

      // need to get layers of bridge-adj entities
    if (from_ents.empty()) continue;
    result = MeshTopoUtil(mbImpl).get_bridge_adjacencies(from_ents, bridge_dim,
                                                         ghost_dim, ghosted_ents, 
                                                         num_layers);
    RRA("Couldn't get bridge adjacencies.");
  }
  
  result = add_verts(ghosted_ents);
  RRA("Couldn't add verts.");

  return result;
}

ErrorCode ParallelComm::add_verts(Range &sent_ents) 
{
      // get the verts adj to these entities, since we'll have to send those too

    // first check sets
  std::pair<Range::const_iterator, Range::const_iterator>
      set_range = sent_ents.equal_range(MBENTITYSET);
  ErrorCode result = MB_SUCCESS, tmp_result;
  for (Range::const_iterator rit = set_range.first; rit != set_range.second; rit++) {
    tmp_result = mbImpl->get_entities_by_type(*rit, MBVERTEX, sent_ents);
    if (MB_SUCCESS != tmp_result) result = tmp_result;
  }
  RRA("Failed to get contained verts.");
  
    // now non-sets
  Range tmp_ents;
  std::copy(sent_ents.begin(), set_range.first, range_inserter(tmp_ents));
  result = mbImpl->get_adjacencies(tmp_ents, 0, false, sent_ents,
                                   Interface::UNION);
  RRA("Couldn't get vertices adj to ghosted ents.");

  return result;
}


ErrorCode ParallelComm::exchange_tags( const std::vector<Tag> &src_tags,
                                       const std::vector<Tag> &dst_tags,
                                       const Range &entities_in)
{
  ErrorCode result;
  int success;

#ifdef DEBUG_COMM
  std::cerr << "Entering exchange_tags" << std::endl; std::cerr.flush();
#endif

    // get all procs interfacing to this proc
  std::set<unsigned int> exch_procs;
  result = get_comm_procs(exch_procs);  

    // post ghost irecv's for all interface procs
    // index greqs the same as buffer/sharing procs indices
  std::vector<MPI_Request> recv_tag_reqs(2*buffProcs.size(), MPI_REQUEST_NULL),
      sent_ack_reqs(buffProcs.size(), MPI_REQUEST_NULL);
  std::vector<unsigned int>::iterator sit;
  int ind;

  reset_all_buffers();
  int incoming = 0;

  for (ind = 0, sit = buffProcs.begin(); sit != buffProcs.end(); sit++, ind++) {
    incoming++;
    PRINT_DEBUG_IRECV(*sit, procConfig.proc_rank(), remoteOwnedBuffs[ind]->mem_ptr,
                      INITIAL_BUFF_SIZE, MB_MESG_TAGS_SIZE, incoming);

    success = MPI_Irecv(remoteOwnedBuffs[ind]->mem_ptr, INITIAL_BUFF_SIZE,
                        MPI_UNSIGNED_CHAR, *sit,
                        MB_MESG_TAGS_SIZE, procConfig.proc_comm(), 
                        &recv_tag_reqs[2*ind]);
    if (success != MPI_SUCCESS) {
      result = MB_FAILURE;
      RRA("Failed to post irecv in ghost exchange.");
    }

  }
  
    // pack and send tags from this proc to others
    // make sendReqs vector to simplify initialization
  sendReqs.resize(2*buffProcs.size(), MPI_REQUEST_NULL);
  
    // take all shared entities if incoming list is empty
  const Range& entities = entities_in.empty() ? sharedEnts : entities_in;
  
  int dum_ack_buff;

  for (ind = 0, sit = buffProcs.begin(); sit != buffProcs.end(); sit++, ind++) {
    
    Range tag_ents = entities;
    
      // get ents shared by proc *sit
    result = filter_pstatus(tag_ents, PSTATUS_SHARED, PSTATUS_AND, *sit);
    RRA("Failed pstatus AND check.");
    
      // remote nonowned entities
    if (!tag_ents.empty()) {
      result = filter_pstatus(tag_ents, PSTATUS_NOT_OWNED, PSTATUS_NOT);
      RRA("Failed pstatus NOT check.");
    }
    
      // pack-send; this also posts receives if store_remote_handles is true
    std::vector<Range> tag_ranges;
    for (std::vector<Tag>::const_iterator vit = src_tags.begin(); vit != src_tags.end(); vit++) {
      const void* ptr;
      int size;
      if (tagServer->get_default_data_ref( *vit, ptr, size ) != MB_SUCCESS) {
        Range tagged_ents;
        tagServer->get_entities( *vit, tagged_ents );
        tag_ranges.push_back( intersect( tag_ents, tagged_ents ) );
      } 
      else {
        tag_ranges.push_back(tag_ents);
      }
    }
    
      // pack the data
      // reserve space on front for size and for initial buff size
    localOwnedBuffs[ind]->reset_ptr(sizeof(int));
    
    result = pack_tags(tag_ents,
                       src_tags, dst_tags, tag_ranges, 
                       localOwnedBuffs[ind], true, *sit);
    RRA("Failed to count buffer in pack_send_tag.");

      // now send it
    result = send_buffer(*sit, localOwnedBuffs[ind], MB_MESG_TAGS_SIZE, sendReqs[2*ind],
                         recv_tag_reqs[2*ind+1], &dum_ack_buff, incoming);
    RRA("Failed to send buffer.");
                         
  }
  
    // receive/unpack tags
  while (incoming) {
    int ind;
    MPI_Status status;
    PRINT_DEBUG_WAITANY(recv_tag_reqs, MB_MESG_TAGS_SIZE, procConfig.proc_rank());
    success = MPI_Waitany(2*buffProcs.size(), &recv_tag_reqs[0], &ind, &status);
    if (MPI_SUCCESS != success) {
      result = MB_FAILURE;
      RRA("Failed in waitany in ghost exchange.");
    }
    
    PRINT_DEBUG_RECD(status);

      // ok, received something; decrement incoming counter
    incoming--;
    
    bool done = false;
    Range dum_range;
    result = recv_buffer(MB_MESG_TAGS_SIZE,
                         status,
                         remoteOwnedBuffs[ind/2],
                         recv_tag_reqs[ind/2 * 2], recv_tag_reqs[ind/2 * 2 + 1],
                         incoming,
                         localOwnedBuffs[ind/2], sendReqs[ind/2*2], sendReqs[ind/2*2+1],
                         done);
    RRA("Failed to resize recv buffer.");
    if (done) {
      remoteOwnedBuffs[ind/2]->reset_ptr(sizeof(int));
      result = unpack_tags(remoteOwnedBuffs[ind/2],
                           dum_range, true, buffProcs[ind/2]);
      RRA("Failed to recv-unpack-tag message.");
    }
  }
  
    // ok, now wait
#ifdef DEBUG_BARRIER  
  success = MPI_Barrier(procConfig.proc_comm());
#else
  MPI_Status status;
  success = MPI_Waitall(2*buffProcs.size(), &sendReqs[0], status);
#endif
  if (MPI_SUCCESS != success) {
    result = MB_FAILURE;
    RRA("Failure in waitall in tag exchange.");
  }
  
    // If source tag is not equal to destination tag, then
    // do local copy for owned entities (communicate w/ self)
  assert(src_tags.size() == dst_tags.size());
  if (src_tags != dst_tags) {
    std::vector<unsigned char> data;
    Range owned_ents(entities_in);
    result = filter_pstatus(owned_ents, PSTATUS_NOT_OWNED, PSTATUS_NOT);
    RRA("Failure to get subset of owned entities");
  
    if (!owned_ents.empty()) { // check this here, otherwise we get 
      // unexpected results from get_entities_by_type_and_tag w/ Interface::INTERSECT
  
      for (size_t i = 0; i < src_tags.size(); ++i) {
        if (src_tags[i] == dst_tags[i])
          continue;

        Range tagged_ents(owned_ents);
        result = mbImpl->get_entities_by_type_and_tag( 0, MBMAXTYPE,
                          &src_tags[0], 0, 1, tagged_ents, Interface::INTERSECT );
        RRA("get_entities_by_type_and_tag(type == MBMAXTYPE) failed.");

        int size, size2;
        result = mbImpl->tag_get_size( src_tags[i], size );
        RRA("tag_get_size failed.");
        result = mbImpl->tag_get_size( dst_tags[i], size2 );
        RRA("tag_get_size failed.");
        if (size != size2) {
          result = MB_FAILURE;
          RRA("tag sizes don't match")
        }

        data.resize( size * tagged_ents.size() );
        result = mbImpl->tag_get_data( src_tags[i], tagged_ents, &data[0] );
        RRA("tag_get_data failed.");
        result = mbImpl->tag_set_data( dst_tags[i], tagged_ents, &data[0] );
        RRA("tag_set_data failed.");
      }
    }
  }
  
#ifdef DEBUG_COMM
  std::cerr << "Exiting exchange_tags" << std::endl; std::cerr.flush();
#endif

  return MB_SUCCESS;
}

/*
ErrorCode ParallelComm::exchange_tags( Tag src_tag, 
                                           Tag dst_tag, 
                                           const Range& entities )
{
  ErrorCode result;
  int success;

    // get all procs interfacing to this proc
  std::set<unsigned int> exch_procs;
  result = get_comm_procs(exch_procs);  

    // post ghost irecv's for all interface procs
    // index greqs the same as buffer/sharing procs indices
  std::vector<MPI_Request> recv_reqs(MAX_SHARING_PROCS, MPI_REQUEST_NULL);
  std::vector<MPI_Status> gstatus(MAX_SHARING_PROCS);
  std::vector<unsigned int>::iterator sit;
  int ind;
  for (ind = 0, sit = buffProcs.begin(); sit != buffProcs.end(); sit++, ind++) {
    success = MPI_Irecv(&ghostRBuffs[ind][0], ghostRBuffs[ind].size(), 
                        MPI_UNSIGNED_CHAR, *sit,
                        MB_MESG_ANY, procConfig.proc_comm(), 
                        &recv_reqs[ind]);
    if (success != MPI_SUCCESS) {
      result = MB_FAILURE;
      RRA("Failed to post irecv in ghost exchange.");
    }
  }
  
    // figure out which entities are shared with which processors
  std::map<int,Range> proc_ents;
  int other_procs[MAX_SHARING_PROCS], num_sharing;
  for (Range::const_iterator i = entities.begin(); i != entities.end(); ++i) {
    int owner;
    result = get_owner( *i, owner );
    RRA("Failed to get entity owner.");

      // only send entities that this proc owns
    if ((unsigned)owner != proc_config().proc_rank()) 
      continue;
    
    result = get_sharing_parts( *i, other_procs, num_sharing );
    RRA("Failed to get procs sharing entity.");
    if (num_sharing == 0) // keep track of non-shared entities for later
      proc_ents[proc_config().proc_rank()].insert( *i );
    for (int j = 0; j < num_sharing; ++j)
      proc_ents[other_procs[j]].insert( *i );
  }
  
    // pack and send tags from this proc to others
    // make sendReqs vector to simplify initialization
  std::fill(sendReqs, sendReqs+MAX_SHARING_PROCS, MPI_REQUEST_NULL);
  std::map<unsigned int,Range>::const_iterator mit;
  
  for (ind = 0, sit = buffProcs.begin(); sit != buffProcs.end(); sit++, ind++) {
    
      // count first
      // buffer needs to begin with the number of tags (one)
    int buff_size = sizeof(int);
    result = packed_tag_size( src_tag, proc_ents[*sit], buff_size );
    RRA("Failed to count buffer in pack_send_tag.");

    unsigned char *buff_ptr = &ownerSBuffs[ind][0];
    buff->check_space(ownerSBuffs[ind], buff_ptr, buff_size);
    PACK_INT( buff_ptr, 1 ); // number of tags
    result = pack_tag( src_tag, dst_tag, proc_ents[*sit], proc_ents[*sit],
                       ownerSBuffs[ind], buff_ptr, true, *sit );
    RRA("Failed to pack buffer in pack_send_tag.");

      // if the message is large, send a first message to tell how large
    if (INITIAL_BUFF_SIZE < buff_size) {
      int tmp_buff_size = -buff_size;
      int success = MPI_Send(&tmp_buff_size, sizeof(int), MPI_UNSIGNED_CHAR, 
                             *sit, MB_MESG_SIZE, procConfig.proc_comm());
      if (success != MPI_SUCCESS) return MB_FAILURE;
    }
    
      // send the buffer
    success = MPI_Isend(&ownerSBuffs[ind][0], buff_size, MPI_UNSIGNED_CHAR, *sit, 
                        MB_MESG_TAGS, procConfig.proc_comm(), &sendReqs[ind]);
    if (success != MPI_SUCCESS) return MB_FAILURE;
  }
  
    // receive/unpack tags
  int num_incoming = exch_procs.size();
  
  while (num_incoming) {
    int ind;
    MPI_Status status;
    success = MPI_Waitany(MAX_SHARING_PROCS, &recv_reqs[0], &ind, &status);
    if (MPI_SUCCESS != success) {
      result = MB_FAILURE;
      RRA("Failed in waitany in ghost exchange.");
    }
    
      // ok, received something; decrement incoming counter
    num_incoming--;
    
    int new_size;
    unsigned char *buff_ptr;
    Range dum_range;
    
      // branch on message type
    switch (status.MPI_TAG) {
      case MB_MESG_SIZE:
          // incoming message just has size; resize buffer and re-call recv,
          // then re-increment incoming count
        assert(ind < MAX_SHARING_PROCS);
        new_size = *((int*)&ghostRBuffs[ind][0]);
        assert(0 > new_size);
        result = recv_size_buff(buffProcs[ind], ghostRBuffs[ind], recv_reqs[ind],
                                MB_MESG_TAGS);
        RRA("Failed to resize recv buffer.");
        num_incoming++;
        break;
      case MB_MESG_TAGS:
          // incoming ghost entities; process
          buff_ptr = &ghostRBuffs[ind][0];
          result = unpack_tags(buff_ptr, dum_range, true,
                               buffProcs[ind]);
        RRA("Failed to recv-unpack-tag message.");
        break;
      default:
        result = MB_FAILURE;
        RRA("Failed to get message of correct type in exch_tags.");
        break;
    }
  }
  
    // ok, now wait
  MPI_Status status[MAX_SHARING_PROCS];
  success = MPI_Waitall(MAX_SHARING_PROCS, &sendReqs[0], status);
  if (MPI_SUCCESS != success) {
    result = MB_FAILURE;
    RRA("Failure in waitall in tag exchange.");
  }
  
    // if src and destination tags aren't the same, need to copy 
    // values for local entities
  if (src_tag != dst_tag) {
    const Range& myents = proc_ents[proc_config().proc_rank()];
    std::vector<const void*> data_ptrs(myents.size());
    std::vector<int> data_sizes(myents.size());
    result = get_moab()->tag_get_data( src_tag, myents, &data_ptrs[0], &data_sizes[0] );
    RRA("Failure to get pointers to local data.");
    result = get_moab()->tag_set_data( dst_tag, myents, &data_ptrs[0], &data_sizes[0] );
    RRA("Failure to get pointers to local data.");
  }  
  
  return MB_SUCCESS;
}
*/

ErrorCode ParallelComm::update_shared_mesh()
{
//  ErrorCode result;
//  int success;

    // ,,,
    /*

    // get all procs interfacing to this proc
  std::set<unsigned int> iface_procs;
  result = get_interface_procs(iface_procs);
  RRA("Failed to get iface sets, procs");

    // post ghost irecv's for all interface procs
    // index greqs the same as buffer/sharing procs indices
  std::vector<MPI_Request> recv_reqs(2*MAX_SHARING_PROCS, MPI_REQUEST_NULL);
  std::vector<MPI_Status> gstatus(MAX_SHARING_PROCS);
  std::set<unsigned int>::iterator sit;
  for (sit = iface_procs.begin(); sit != iface_procs.end(); sit++) {
    int ind = get_buffers(*sit);
    success = MPI_Irecv(&ghostRBuffs[ind][0], ghostRBuffs[ind].size(), 
                        MPI_UNSIGNED_CHAR, *sit,
                        MB_MESG_ANY, procConfig.proc_comm(), 
                        &recv_reqs[ind]);
    if (success != MPI_SUCCESS) {
      result = MB_FAILURE;
      RRA("Failed to post irecv in ghost exchange.");
    }
  }
  
    // pack and send vertex coordinates from this proc to others
    // make sendReqs vector to simplify initialization
  std::fill(sendReqs, sendReqs+2*MAX_SHARING_PROCS, MPI_REQUEST_NULL);
  Range recd_ents[MAX_SHARING_PROCS];
  
  for (sit = iface_procs.begin(); sit != iface_procs.end(); sit++) {
    int ind = get_buffers(*sit);
    
    Range vertices;
    for (Range::iterator rit = interfaceSets.begin(); rit != interfaceSets.end();
         rit++) {
      if (!is_iface_proc(*rit, *sit)) 
        continue;
      
      result = mbImpl->get_entities_by_type( *rit, MBVERTEX, vertices );
      RRA("Bad interface set.");
    }
    std::map<unsigned int,Range>::iterator ghosted = ghostedEnts.find(*sit);
    if (ghosted != ghostedEnts.end()) {
      Range::iterator e = ghosted->second.upper_bound(MBVERTEX);
      vertices.merge( ghosted->second.begin(), e );
    }

      // pack-send; this also posts receives if store_remote_handles is true
    Range sent;
    result = pack_send_entities(*sit, vertices, false, false, 
                                false, true,
                                ownerSBuffs[ind], ownerRBuffs[MAX_SHARING_PROCS+ind], 
                                sendReqs[ind], recv_reqs[MAX_SHARING_PROCS+ind], 
                                sent);
    RRA("Failed to pack-send in mesh update exchange.");
  }
  
    // receive/unpack entities
    // number of incoming messages depends on whether we're getting back
    // remote handles
  int num_incoming = iface_procs.size();
  
  while (num_incoming) {
    int ind;
    MPI_Status status;
    success = MPI_Waitany(2*MAX_SHARING_PROCS, &recv_reqs[0], &ind, &status);
    if (MPI_SUCCESS != success) {
      result = MB_FAILURE;
      RRA("Failed in waitany in ghost exchange.");
    }
    
      // ok, received something; decrement incoming counter
    num_incoming--;
    
    std::vector<EntityHandle> remote_handles_v, sent_ents_tmp;
    Range remote_handles_r;
    int new_size;
    
      // branch on message type
    switch (status.MPI_TAG) {
      case MB_MESG_SIZE:
          // incoming message just has size; resize buffer and re-call recv,
          // then re-increment incoming count
        assert(ind < MAX_SHARING_PROCS);
        new_size = *((int*)&ghostRBuffs[ind][0]);
        assert(0 > new_size);
        result = recv_size_buff(buffProcs[ind], ghostRBuffs[ind], recv_reqs[ind],
                                MB_MESG_ENTS);
        RRA("Failed to resize recv buffer.");
        num_incoming++;
        break;
      case MB_MESG_ENTS:
          // incoming ghost entities; process
        result = recv_unpack_entities(buffProcs[ind], true,
                                      false, 
                                      ghostRBuffs[ind], ghostSBuffs[ind], 
                                      sendReqs[ind], recd_ents[ind]);
        RRA("Failed to recv-unpack message.");
        break;
    }
  }
  
    // ok, now wait if requested
  MPI_Status status[2*MAX_SHARING_PROCS];
  success = MPI_Waitall(2*MAX_SHARING_PROCS, &sendReqs[0], status);
  if (MPI_SUCCESS != success) {
    result = MB_FAILURE;
    RRA("Failure in waitall in ghost exchange.");
  }
  
  return MB_SUCCESS;
}
ErrorCode ParallelComm::update_iface_sets(Range &sent_ents,
                                              std::vector<EntityHandle> &remote_handles, 
                                              int from_proc) 
{
  std::vector<EntityHandle>::iterator remote_it = remote_handles.begin();
  Range::iterator sent_it = sent_ents.begin();
  Range ents_to_remove;
  for (; sent_it != sent_ents.end(); sent_it++, remote_it++) {
    if (!*remote_it) ents_to_remove.insert(*sent_it);
  }
  
  for (Range::iterator set_it = interfaceSets.begin(); set_it != interfaceSets.end(); set_it++) {
    if (!is_iface_proc(*set_it, from_proc)) continue;
    ErrorCode result = mbImpl->remove_entities(*set_it, ents_to_remove);
    RRA("Couldn't remove entities from iface set in update_iface_sets.");
  }

*/
  
  return MB_SUCCESS;
}

  //! return sharedp tag
Tag ParallelComm::sharedp_tag()
{
  if (!sharedpTag) {
    int def_val = -1;
    ErrorCode result = mbImpl->tag_create(PARALLEL_SHARED_PROC_TAG_NAME, 
                                            sizeof(int), 
                                            MB_TAG_DENSE,
                                            MB_TYPE_INTEGER, sharedpTag, 
                                            &def_val, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) 
      return 0;
  }
  
  return sharedpTag;
}

  //! return sharedps tag
Tag ParallelComm::sharedps_tag()
{
  if (!sharedpsTag) {
    ErrorCode result = mbImpl->tag_create(PARALLEL_SHARED_PROCS_TAG_NAME, 
                                            MAX_SHARING_PROCS*sizeof(int), 
                                            MB_TAG_SPARSE,
                                            MB_TYPE_INTEGER, sharedpsTag, NULL, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) 
      return 0;
  }
  
  return sharedpsTag;
}
  
  //! return sharedh tag
Tag ParallelComm::sharedh_tag()
{
  if (!sharedhTag) {
    EntityHandle def_val = 0;
    ErrorCode result = mbImpl->tag_create(PARALLEL_SHARED_HANDLE_TAG_NAME, 
                                            sizeof(EntityHandle), 
                                            MB_TAG_DENSE,
                                            MB_TYPE_HANDLE, sharedhTag, 
                                            &def_val, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result)
      return 0;
  }
  
  return sharedhTag;
}
  
  //! return sharedhs tag
Tag ParallelComm::sharedhs_tag()
{  
  if (!sharedhsTag) {
    ErrorCode result = mbImpl->tag_create(PARALLEL_SHARED_HANDLES_TAG_NAME, 
                                            MAX_SHARING_PROCS*sizeof(EntityHandle), 
                                            MB_TAG_SPARSE,
                                            MB_TYPE_INTEGER, sharedhsTag, NULL, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) 
      return 0;
  }

  return sharedhsTag;
}
  
  //! return pstatus tag
Tag ParallelComm::pstatus_tag()
{  
  if (!pstatusTag) {
    unsigned char tmp_pstatus = 0;
    ErrorCode result = mbImpl->tag_create(PARALLEL_STATUS_TAG_NAME, 
                                            sizeof(unsigned char),
                                            MB_TAG_DENSE,
                                            MB_TYPE_OPAQUE, pstatusTag, 
                                            &tmp_pstatus, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result)
      return 0;
  }
  
  return pstatusTag;
}
  
  //! return partition set tag
Tag ParallelComm::partition_tag()
{  
  if (!partitionTag) {
    ErrorCode result = mbImpl->tag_create(PARALLEL_PARTITION_TAG_NAME, 
                                            sizeof(int),
                                            MB_TAG_SPARSE,
                                            MB_TYPE_INTEGER, 
                                            partitionTag, 
                                            NULL, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result)
      return 0;
  }
  
  return partitionTag;
}
  
  //! return pcomm tag; passes in impl 'cuz this is a static function
Tag ParallelComm::pcomm_tag(Interface *impl,
                                bool create_if_missing)
{
  Tag this_tag = 0;
  ErrorCode result;
  result = impl->tag_get_handle(PARALLEL_COMM_TAG_NAME, this_tag);
  if ((MB_TAG_NOT_FOUND == result || 0 == this_tag) &&
      create_if_missing) {
    result = impl->tag_create(PARALLEL_COMM_TAG_NAME, 
                              MAX_SHARING_PROCS*sizeof(ParallelComm*),
                              MB_TAG_SPARSE,
                              MB_TYPE_OPAQUE, this_tag,
                              NULL, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result)
      return 0;
  }
  
  return this_tag;
}

    //! get the indexed pcomm object from the interface
ParallelComm *ParallelComm::get_pcomm(Interface *impl, const int index) 
{
  Tag pc_tag = pcomm_tag(impl, false);
  if (0 == pc_tag) return NULL;
  
  ParallelComm *pc_array[MAX_SHARING_PROCS];
  ErrorCode result = impl->tag_get_data(pc_tag, 0, 0, (void*)pc_array);
  if (MB_SUCCESS != result) return NULL;
  
  return pc_array[index];
}

ErrorCode ParallelComm::get_all_pcomm( Interface* impl, std::vector<ParallelComm*>& list )
{
  Tag pc_tag = pcomm_tag(impl, false);
  if (0 == pc_tag)
    return MB_TAG_NOT_FOUND;
  
  ParallelComm *pc_array[MAX_SHARING_PROCS];
  ErrorCode rval = impl->tag_get_data( pc_tag, 0, 0, pc_array );
  if (MB_SUCCESS != rval)
    return rval;
  
  for (int i = 0; i < MAX_SHARING_PROCS; ++i)
    if (pc_array[i])
      list.push_back( pc_array[i] );
  
  return MB_SUCCESS;
}
  

    //! get the indexed pcomm object from the interface
ParallelComm *ParallelComm::get_pcomm( Interface *impl, 
                                           EntityHandle prtn,
                                           const MPI_Comm* comm ) 
{
  ErrorCode rval;
  ParallelComm* result = 0;
  
  Tag prtn_tag;
  rval = impl->tag_create( PARTITIONING_PCOMM_TAG_NAME, 
                           sizeof(int),
                           MB_TAG_SPARSE,
                           MB_TYPE_INTEGER,
                           prtn_tag,
                           0, true );
  if (MB_SUCCESS != rval)
    return 0;
  
  int pcomm_id;
  rval = impl->tag_get_data( prtn_tag, &prtn, 1, &pcomm_id );
  if (MB_SUCCESS == rval) {
    result= get_pcomm( impl, pcomm_id );
  }
  else if (MB_TAG_NOT_FOUND == rval && comm) {
    result = new ParallelComm( impl, *comm, &pcomm_id );
    if (!result)
      return 0;
    result->set_partitioning( prtn );
    
    rval = impl->tag_set_data( prtn_tag, &prtn, 1, &pcomm_id );
    if (MB_SUCCESS != rval) {
      delete result;
      result = 0;
    }
  }
  
  return result;
}

ErrorCode ParallelComm::set_partitioning( EntityHandle set) 
{
  ErrorCode rval;
  Tag prtn_tag;
  rval = mbImpl->tag_create( PARTITIONING_PCOMM_TAG_NAME, 
                           sizeof(int),
                           MB_TAG_SPARSE,
                           MB_TYPE_INTEGER,
                           prtn_tag,
                           0, true );
  if (MB_SUCCESS != rval)
    return rval;

    // get my id
  ParallelComm* pcomm_arr[MAX_SHARING_PROCS];
  Tag pc_tag = pcomm_tag(mbImpl, false);
  if (0 == pc_tag) 
    return MB_FAILURE;
  ErrorCode result = mbImpl->tag_get_data(pc_tag, 0, 0, pcomm_arr);
  if (MB_SUCCESS != result) 
    return MB_FAILURE;  
  int id = std::find(pcomm_arr,pcomm_arr+MAX_SHARING_PROCS,this) - pcomm_arr;
  if (id == MAX_SHARING_PROCS)
    return MB_FAILURE;

  EntityHandle old = partitioningSet;
  if (old) {
    rval = mbImpl->tag_delete_data( prtn_tag, &old, 1 );
    if (MB_SUCCESS != rval)
      return rval;
    partitioningSet = 0;
  }
  
  if (!set) 
    return MB_SUCCESS;
  
  Range contents;
  if (old) {
    rval = mbImpl->get_entities_by_handle( old, contents );
    if (MB_SUCCESS != rval)
      return rval;
  }
  else {
    contents = partition_sets();
  }

  rval = mbImpl->add_entities( set, contents );
  if (MB_SUCCESS != rval)
    return rval;
  
    // store pcomm id on new partition set
  rval = mbImpl->tag_set_data( prtn_tag, &set, 1, &id );
  if (MB_SUCCESS != rval)
    return rval;
  
  partitioningSet = set;
  return MB_SUCCESS;
}
  

  //! return all the entities in parts owned locally
ErrorCode ParallelComm::get_part_entities(Range &ents, int dim) 
{
  ErrorCode result;
  
  for (Range::iterator rit = partitionSets.begin(); 
       rit != partitionSets.end(); rit++) {
    Range tmp_ents;
    if (-1 == dim) 
      result = mbImpl->get_entities_by_handle(*rit, tmp_ents, true);
    else
      result = mbImpl->get_entities_by_dimension(*rit, dim, tmp_ents, true);

    if (MB_SUCCESS != result) return result;
    ents.merge(tmp_ents);
  }
  
  return MB_SUCCESS;
}

  /** \brief Return the rank of the entity owner
   */
ErrorCode ParallelComm::get_owner_handle(EntityHandle entity,
                                             int &owner,
                                             EntityHandle &handle) 
{
  unsigned char pstat;
  int sharing_procs[MAX_SHARING_PROCS];
  EntityHandle sharing_handles[MAX_SHARING_PROCS];

  ErrorCode result = mbImpl->tag_get_data(pstatus_tag(), &entity, 1,
                                            &pstat);
  if (!(pstat & PSTATUS_NOT_OWNED)) {
    owner = proc_config().proc_rank();
    handle = entity;
  }
  
  else if (pstat & PSTATUS_MULTISHARED) {
    result = mbImpl->tag_get_data(sharedps_tag(), &entity, 1,
                                  sharing_procs);
    owner = sharing_procs[0];
    result = mbImpl->tag_get_data(sharedhs_tag(), &entity, 1,
                                  sharing_handles);
    handle = sharing_handles[0];
  }
  else if (pstat & PSTATUS_SHARED) {
    result = mbImpl->tag_get_data(sharedp_tag(), &entity, 1,
                                  sharing_procs);
    RRA(" ");
    owner = sharing_procs[0];
    result = mbImpl->tag_get_data(sharedh_tag(), &entity, 1,
                                  sharing_handles);
    handle = sharing_handles[0];
  }
  else {
    owner = -1;
    handle = 0;
  }
  
  return MB_SUCCESS;
}

ErrorCode ParallelComm::get_global_part_count( int& count_out ) const
{
  count_out = globalPartCount;
  return count_out < 0 ? MB_FAILURE : MB_SUCCESS;
}

ErrorCode ParallelComm::get_part_owner( int part_id, int& owner ) const
{
  // FIXME: assumes only 1 local part
  owner = part_id;
  return MB_SUCCESS;
}

ErrorCode ParallelComm::get_part_id( EntityHandle /*part*/, int& id_out ) const
{
  // FIXME: assumes only 1 local part
  id_out = proc_config().proc_rank();
  return MB_SUCCESS;
}

ErrorCode ParallelComm::get_part_handle( int id, EntityHandle& handle_out ) const
{
  // FIXME: assumes only 1 local part
  if ((unsigned)id != proc_config().proc_rank())
    return MB_ENTITY_NOT_FOUND;
  handle_out = partition_sets().front();
  return MB_SUCCESS;
}

ErrorCode ParallelComm::create_part( EntityHandle& set_out )
{
    // mark as invalid so we know that it needs to be updated
  globalPartCount = -1;
  
    // create set representing part
  ErrorCode rval = mbImpl->create_meshset( MESHSET_SET, set_out );
  if (MB_SUCCESS != rval)
    return rval;
  
    // set tag on set
    // FIXME: need to assign valid global id
  int val = 0;
  rval = mbImpl->tag_set_data( part_tag(), &set_out, 1, &val );
  if (MB_SUCCESS != rval) {
    mbImpl->delete_entities( &set_out, 1 );
    return rval;
  }
  
  if (get_partitioning()) {
    rval = mbImpl->add_entities( get_partitioning(), &set_out, 1 );
    if (MB_SUCCESS != rval) {
      mbImpl->delete_entities( &set_out, 1 );
      return rval;
    }
  }
  
  return MB_SUCCESS;
}

ErrorCode ParallelComm::destroy_part( EntityHandle part_id )
{
    // mark as invalid so we know that it needs to be updated
  globalPartCount = -1;
  
  ErrorCode rval;
  if (get_partitioning()) {
    rval = mbImpl->remove_entities( get_partitioning(), &part_id, 1 );
    if (MB_SUCCESS != rval)
      return rval;
  }
  return mbImpl->delete_entities( &part_id, 1 );
}

ErrorCode ParallelComm::collective_sync_partition()
{
  int count = partition_sets().size();
  globalPartCount = 0;
  int err = MPI_Allreduce( &count, &globalPartCount, 1, MPI_INT, MPI_SUM, 
                           proc_config().proc_comm() );
  return err ? MB_FAILURE : MB_SUCCESS;
}

ErrorCode ParallelComm::get_part_neighbor_ids( EntityHandle part,
                                                   int neighbors_out[MAX_SHARING_PROCS],
                                                   int& num_neighbors_out )
{
  ErrorCode rval;
  Range iface;
  rval = get_interface_sets( part, iface );
  if (MB_SUCCESS != rval)
    return rval;
  
  num_neighbors_out = 0;
  int n, j = 0;
  int tmp[MAX_SHARING_PROCS], curr[MAX_SHARING_PROCS];
  int *parts[2] = { neighbors_out, tmp };
  for (Range::iterator i = iface.begin(); i != iface.end(); ++i) {
    unsigned char pstat;
    rval = get_sharing_data( *i, curr, NULL, pstat, n);
    if (MB_SUCCESS != rval)
      return rval;
    std::sort( curr, curr+n );
    assert( num_neighbors_out < MAX_SHARING_PROCS );
    int* k = std::set_union( parts[j], parts[j]+num_neighbors_out,
                             curr, curr + n, parts[1-j] );
    j = 1-j;
    num_neighbors_out = k - parts[j];
  }
  if (parts[j] != neighbors_out)
    std::copy( parts[j], parts[j]+num_neighbors_out, neighbors_out );
    
    
    // remove input part from list
  int id;
  rval = get_part_id( part, id );
  if (MB_SUCCESS == rval) 
    num_neighbors_out = std::remove( neighbors_out, neighbors_out+num_neighbors_out, id ) - neighbors_out;
  return rval;
}

ErrorCode ParallelComm::get_interface_sets( EntityHandle ,
                                                Range& iface_sets_out,
                                                int* adj_part_id )
{
    // FIXME : assumes one part per processor.
    // Need to store part iface sets as children to implement
    // this correctly.
  iface_sets_out = interface_sets();

  if (adj_part_id) {
    int part_ids[MAX_SHARING_PROCS], num_parts;
    Range::iterator i = iface_sets_out.begin();
    while (i != iface_sets_out.end()) {
      unsigned char pstat;
      ErrorCode rval = get_sharing_data( *i, part_ids, NULL, pstat, num_parts );
      if (MB_SUCCESS != rval)
        return rval;
      
      if (std::find(part_ids, part_ids+num_parts, *adj_part_id) - part_ids != num_parts)
        ++i;
      else
        i = iface_sets_out.erase( i );
    }
  }
    
  return MB_SUCCESS;
}

ErrorCode ParallelComm::get_owning_part( EntityHandle handle,
                                             int& owning_part_id,
                                             EntityHandle* remote_handle )
{

  // FIXME : assumes one part per proc, and therefore part_id == rank
  
    // If entity is not shared, then we're the owner.
  unsigned char pstat;
  ErrorCode result = mbImpl->tag_get_data(pstatus_tag(), &handle, 1,
                                            &pstat);
  if (!(pstat & PSTATUS_NOT_OWNED)) {
    owning_part_id = proc_config().proc_rank();
    if (remote_handle)
      *remote_handle = handle;
    return MB_SUCCESS;
  }
  
    // If entity is shared with one other proc, then
    // sharedp_tag will contain a positive value.
  result = mbImpl->tag_get_data( sharedp_tag(), &handle, 1, &owning_part_id );
  if (MB_SUCCESS != result)
    return result;
  if (owning_part_id != -1) {
      // done?
    if (!remote_handle)
      return MB_SUCCESS;
      
      // get handles on remote processors (and this one)
    return mbImpl->tag_get_data( sharedh_tag(), &handle, 1, remote_handle );
  }
  
    // If here, then the entity is shared with at least two other processors.
    // Get the list from the sharedps_tag
  const void* part_id_list = 0;
  result = mbImpl->tag_get_data( sharedps_tag(), &handle, 1, &part_id_list );
  if (MB_SUCCESS != result)
    return result;
  owning_part_id = ((const int*)part_id_list)[0];
 
    // done?
  if (!remote_handle)
    return MB_SUCCESS;
  
    // get remote handles
  const void* handle_list = 0;
  result = mbImpl->tag_get_data( sharedhs_tag(), &handle, 1, &handle_list );
  if (MB_SUCCESS != result)
    return result;
  
  *remote_handle = ((const EntityHandle*)handle_list)[0];
  return MB_SUCCESS;
}    

ErrorCode ParallelComm::get_sharing_parts( EntityHandle entity,
                                               int part_ids_out[MAX_SHARING_PROCS],
                                               int& num_part_ids_out,
                                               EntityHandle remote_handles[MAX_SHARING_PROCS] )
{

  // FIXME : assumes one part per proc, and therefore part_id == rank
  
    // If entity is not shared, then we're the owner.
  unsigned char pstat;
  ErrorCode result = mbImpl->tag_get_data(pstatus_tag(), &entity, 1,
                                            &pstat);
  if (!(pstat & PSTATUS_SHARED)) {
    part_ids_out[0] = proc_config().proc_rank();
    if (remote_handles)
      remote_handles[0] = entity;
    num_part_ids_out = 1;
    return MB_SUCCESS;
  }
  
    // If entity is shared with one other proc, then
    // sharedp_tag will contain a positive value.
  result = mbImpl->tag_get_data( sharedp_tag(), &entity, 1, part_ids_out );
  if (MB_SUCCESS != result)
    return result;
  if (part_ids_out[0] != -1) {
    
    num_part_ids_out = 2;
    part_ids_out[1] = proc_config().proc_rank();

      // done?
    if (!remote_handles)
      return MB_SUCCESS;
      
      // get handles on remote processors (and this one)
    remote_handles[1] = entity;
    return mbImpl->tag_get_data( sharedh_tag(), &entity, 1, remote_handles );
  }
  
    // If here, then the entity is shared with at least two other processors.
    // Get the list from the sharedps_tag
  result = mbImpl->tag_get_data( sharedps_tag(), &entity, 1, part_ids_out );
  if (MB_SUCCESS != result)
    return result;
    // Count number of valid (positive) entries in sharedps_tag
  for (num_part_ids_out = 0; num_part_ids_out < MAX_SHARING_PROCS &&
       part_ids_out[num_part_ids_out] >= 0; ++num_part_ids_out);
  //part_ids_out[num_part_ids_out++] = proc_config().proc_rank();
#ifndef NDEBUG
  int my_idx = std::find(part_ids_out, part_ids_out+num_part_ids_out, proc_config().proc_rank()) - part_ids_out;
  assert(my_idx < num_part_ids_out);
#endif
  
    // done?
  if (!remote_handles)
    return MB_SUCCESS;
  
    // get remote handles
  result = mbImpl->tag_get_data( sharedhs_tag(), &entity, 1, remote_handles );
  //remote_handles[num_part_ids_out-1] = entity;
  assert(remote_handles[my_idx] == entity);

  return result;
}

ErrorCode ParallelComm::pack_shared_handles(
    std::vector<std::vector<SharedEntityData> > &send_data) 
{
    // get all shared entities
  Range all_shared, dum_range;
  ErrorCode rval = mbImpl->get_entities_by_handle(0, all_shared);
  if (MB_SUCCESS != rval)
    return rval;
  rval = get_pstatus_entities(-1, 0, dum_range);
  if (MB_SUCCESS != rval)
    return rval;
  all_shared = subtract( all_shared, dum_range);
  all_shared.erase(all_shared.upper_bound(MBPOLYHEDRON), all_shared.end());
  assert(sharedEnts == all_shared);

    // build up send buffers
  int ent_procs[MAX_SHARING_PROCS];
  EntityHandle handles[MAX_SHARING_PROCS];
  int num_sharing, tmp_int;
  SharedEntityData tmp;
  send_data.resize(buffProcs.size());
  for (Range::iterator i = all_shared.begin(); i != all_shared.end(); ++i) {
    tmp.remote = *i; // swap local/remote so they're correct on the remote proc.
    rval = get_owner( *i, tmp_int );
    tmp.owner = tmp_int;
    if (MB_SUCCESS != rval)
      return rval;

    unsigned char pstat;
    rval = get_sharing_data( *i, ent_procs, handles, pstat, num_sharing );
    if (MB_SUCCESS != rval)
      return rval;
    for (int j = 0; j < num_sharing; ++j) {
      if (ent_procs[j] == (int)proc_config().proc_rank())
        continue;
      tmp.local = handles[j];
      int ind = get_buffers(ent_procs[j]);
      assert(-1 != ind);
      if ((int)send_data.size() < ind+1) send_data.resize(ind+1);
      send_data[ind].push_back( tmp );
    }
  }

  return MB_SUCCESS;
}

ErrorCode ParallelComm::exchange_all_shared_handles(  
    std::vector<std::vector<SharedEntityData> > &send_data, 
    std::vector<std::vector<SharedEntityData> > &result)
{
  int ierr;
  const int tag = 0x4A41534E;
  const MPI_Comm comm = procConfig.proc_comm();
  const int num_proc = buffProcs.size();
  const std::vector<int> procs( buffProcs.begin(), buffProcs.end() );
  std::vector<MPI_Request> recv_req(buffProcs.size(), MPI_REQUEST_NULL);
  
    // set up to receive sizes
  std::vector<int> sizes_send(num_proc), sizes_recv(num_proc);
  for (int i = 0; i < num_proc; ++i) {
    ierr = MPI_Irecv( &sizes_recv[i], 1, MPI_INT, procs[i], tag, comm, &recv_req[i] );
    if (ierr) 
      return MB_FILE_WRITE_ERROR;
  }
  
    // send sizes
  assert(num_proc == (int)send_data.size());
  
  sendReqs.resize(buffProcs.size(), MPI_REQUEST_NULL);
  result.resize(num_proc);
  for (int i = 0; i < num_proc; ++i) {
    sizes_send[i] = send_data[i].size();
    ierr = MPI_Isend( &sizes_send[i], 1, MPI_INT, buffProcs[i], tag, comm, &sendReqs[i] );
    if (ierr) 
      return MB_FILE_WRITE_ERROR;
  }
  
    // receive sizes
  std::vector<MPI_Status> stat(num_proc);
  ierr = MPI_Waitall( num_proc, &recv_req[0], &stat[0] );
  if (ierr)
    return MB_FILE_WRITE_ERROR;
  
    // wait until all sizes are sent (clean up pending req's)
  ierr = MPI_Waitall( num_proc, &sendReqs[0], &stat[0] );
  if (ierr)
    return MB_FILE_WRITE_ERROR;
  
    // set up to receive data
  for (int i = 0; i < num_proc; ++i) {
    result[i].resize( sizes_recv[i] );
    ierr = MPI_Irecv( &result[i][0], 
                      sizeof(SharedEntityData)*sizes_recv[i], 
                      MPI_UNSIGNED_CHAR, 
                      buffProcs[i], tag, comm, &recv_req[i] );
    if (ierr) 
      return MB_FILE_WRITE_ERROR;
  }
  
    // send data
  for (int i = 0; i < num_proc; ++i) {
    ierr = MPI_Isend( &send_data[i][0], 
                      sizeof(SharedEntityData)*sizes_send[i], 
                      MPI_UNSIGNED_CHAR, 
                      buffProcs[i], tag, comm, &sendReqs[i] );
    if (ierr) 
      return MB_FILE_WRITE_ERROR;
  }
  
    // receive data
  ierr = MPI_Waitall( num_proc, &recv_req[0], &stat[0] );
  if (ierr)
    return MB_FILE_WRITE_ERROR;
  
    // wait until everything is sent to release send buffers
  ierr = MPI_Waitall( num_proc, &sendReqs[0], &stat[0] );
  if (ierr)
    return MB_FILE_WRITE_ERROR;
  
  return MB_SUCCESS;
}

ErrorCode ParallelComm::check_all_shared_handles(bool print_em) 
{
    // get all shared ent data from other procs
  std::vector<std::vector<SharedEntityData> > shents(buffProcs.size()),
      send_data(buffProcs.size());

  ErrorCode result;
  bool done = false;
  
  while (!done) {
    result = check_local_shared();
    if (MB_SUCCESS != result) {
      done = true;
      continue;
    }

    result = pack_shared_handles(send_data);
    if (MB_SUCCESS != result) {
      done = true;
      continue;
    }
  
    result = exchange_all_shared_handles(send_data, shents);
    if (MB_SUCCESS != result) {
      done = true;
      continue;
    }

    if (!shents.empty()) check_my_shared_handles(shents);
    done = true;
  }
  
  if (MB_SUCCESS != result && print_em) {
    std::ostringstream ent_str;
    ent_str << "mesh." << procConfig.proc_rank() << ".h5m";
    mbImpl->write_mesh(ent_str.str().c_str());
  }
  
  return result;
}

ErrorCode ParallelComm::check_local_shared() 
{
    // do some local checks on shared entities to make sure things look
    // consistent

    // check that non-vertex shared entities are shared by same procs as all
    // their vertices
  std::pair<Range::const_iterator,Range::const_iterator> vert_it =
      sharedEnts.equal_range(MBVERTEX);
  std::vector<EntityHandle> dum_connect;
  const EntityHandle *connect;
  int num_connect;
  int tmp_procs[MAX_SHARING_PROCS];
  EntityHandle tmp_hs[MAX_SHARING_PROCS];
  std::set<int> tmp_set, vset;
  int num_ps;
  ErrorCode result;
  unsigned char pstat;
  Range bad_ents;
  std::vector<std::string> errors;
  std::string dum_err;
  
  Range::const_iterator rit;
  for (rit = sharedEnts.begin(); rit != sharedEnts.end(); rit++) {

      // get sharing procs for this ent
    result = get_sharing_data(*rit, tmp_procs, tmp_hs, pstat, num_ps);
    if (MB_SUCCESS != result) {
      bad_ents.insert(*rit);
      errors.push_back(std::string("Failure getting sharing data."));
      continue;
    }

    bool bad = false;
      // entity must be shared
    if (!(pstat & PSTATUS_SHARED))
      errors.push_back(std::string("Entity should be shared but isn't.")), bad = true;

      // if entity is not owned this must not be first proc
    if (pstat & PSTATUS_NOT_OWNED && tmp_procs[0] == (int)procConfig.proc_rank())
      errors.push_back(std::string("Entity not owned but is first proc.")), bad = true;

      // if entity is owned and multishared, this must be first proc
    if (!(pstat & PSTATUS_NOT_OWNED) && pstat & PSTATUS_MULTISHARED && 
        (tmp_procs[0] != (int)procConfig.proc_rank() || tmp_hs[0] != *rit))
      errors.push_back(std::string("Entity owned and multishared but not first proc or not first handle.")), bad = true;

    if (bad) {
      bad_ents.insert(*rit);
      continue;
    }
    
    if (mbImpl->type_from_handle(*rit) == MBVERTEX) continue;

      // copy element's procs to vset and save size
    int orig_ps = num_ps; vset.clear(); 
    std::copy(tmp_procs, tmp_procs+num_ps, std::inserter(vset, vset.begin()));
    
      // get vertices for this ent and intersection of sharing procs
    result = mbImpl->get_connectivity(*rit, connect, num_connect, false, &dum_connect);
    if (MB_SUCCESS != result) {
      bad_ents.insert(*rit); 
      errors.push_back(std::string("Couldn't get connectivity."));
      continue;
    }
    
    for (int i = 0; i < num_connect; i++) {
      result = get_sharing_data(connect[i], tmp_procs, NULL, pstat, num_ps);
      if (MB_SUCCESS != result) {bad_ents.insert(*rit); continue;}
      if (!num_ps) {vset.clear(); break;}
      std::sort(tmp_procs, tmp_procs+num_ps);
      tmp_set.clear();
      std::set_intersection(tmp_procs, tmp_procs+num_ps,
                            vset.begin(), vset.end(), std::inserter(tmp_set, tmp_set.end()));
      vset.swap(tmp_set);
      if (vset.empty()) break;
    }
    
      // intersect them; should be the same size as orig_ps
    tmp_set.clear();
    std::set_intersection(tmp_procs, tmp_procs+num_ps,
                          vset.begin(), vset.end(), std::inserter(tmp_set, tmp_set.end()));
    if (orig_ps != (int)tmp_set.size()) {
      errors.push_back(std::string("Vertex proc set not same size as entity proc set."));
      bad_ents.insert(*rit);
    }
  }
  
  if (!bad_ents.empty()) {
    std::cout << "Found bad entities in check_local_shared, proc rank "
              << procConfig.proc_rank() << "," << std::endl;
    std::vector<std::string>::iterator vit;
    for (rit = bad_ents.begin(), vit = errors.begin(); rit != bad_ents.end(); rit++, vit++) {
      list_entities(&(*rit), 1);
      std::cout << "Reason: " << *vit << std::endl;
    }
    return MB_FAILURE;
  }

    // to do: check interface sets

  return MB_SUCCESS;
}

ErrorCode ParallelComm::check_all_shared_handles(ParallelComm **pcs,
                                                     int num_pcs) 
{
  std::vector<std::vector<std::vector<SharedEntityData> > > shents, send_data;
  ErrorCode result = MB_SUCCESS, tmp_result;

    // get all shared ent data from each proc to all other procs
  send_data.resize(num_pcs);
  for (int p = 0; p < num_pcs; p++) {
    tmp_result = pcs[p]->pack_shared_handles(send_data[p]);
    if (MB_SUCCESS != tmp_result) result = tmp_result;
  }
  if (MB_SUCCESS != result) return result;

    // move the data sorted by sending proc to data sorted by receiving proc
  shents.resize(num_pcs);
  for (int p = 0; p < num_pcs; p++)
    shents[p].resize(pcs[p]->buffProcs.size());
    
  for (int p = 0; p < num_pcs; p++) {
    for (unsigned int idx_p = 0; idx_p < pcs[p]->buffProcs.size(); idx_p++) {
        // move send_data[p][to_p] to shents[to_p][idx_p]
      int to_p = pcs[p]->buffProcs[idx_p];
      int top_idx_p = pcs[to_p]->get_buffers(p);
      assert(-1 != top_idx_p);
      shents[to_p][top_idx_p] = send_data[p][idx_p];
    }
  }
  
  for (int p = 0; p < num_pcs; p++) {
    std::ostringstream ostr;
    ostr << "Processor " << p << " bad entities:";
    tmp_result = pcs[p]->check_my_shared_handles(shents[p], ostr.str().c_str());
    if (MB_SUCCESS != tmp_result) result = tmp_result;
  }
  
  return result;
}

ErrorCode ParallelComm::check_my_shared_handles(
    std::vector<std::vector<SharedEntityData> > &shents,
                                                    const char *prefix) 
{
    // now check against what I think data should be
    // get all shared entities
  ErrorCode result;
  Range dum_range, all_shared = sharedEnts;
  all_shared.erase(all_shared.upper_bound(MBPOLYHEDRON), all_shared.end());

  Range bad_ents, local_shared;
  std::vector<SharedEntityData>::iterator vit;
  for (unsigned int i = 0; i < shents.size(); i++) {
    int other_proc = buffProcs[i];
    result = get_shared_entities(other_proc, local_shared);
    if (MB_SUCCESS != result) return result;
    for (vit = shents[i].begin(); vit != shents[i].end(); vit++) {
      EntityHandle localh = vit->local, remoteh = vit->remote, dumh;
      local_shared.erase(localh);
      result = get_remote_handles(true, &localh, &dumh, 1, other_proc, dum_range);
      if (MB_SUCCESS != result || dumh != remoteh) 
        bad_ents.insert(localh);
    }

    if (!local_shared.empty()) 
      bad_ents.merge(local_shared);
  }
  
  if (!bad_ents.empty()) {
    if (prefix) {
      std::cout << prefix << std::endl;
      list_entities(bad_ents);
    }
    return MB_FAILURE;
  }

  else return MB_SUCCESS;
}

ErrorCode ParallelComm::get_shared_entities(int other_proc,
                                                Range &shared_ents,
                                                int dim,
                                                const bool iface,
                                                const bool owned_filter) 
{
  shared_ents.clear();
  ErrorCode result = MB_SUCCESS;
  
    // dimension
  if (-1 != dim) {
    DimensionPair dp = CN::TypeDimensionMap[dim];
    Range dum_range;
    shared_ents.merge(sharedEnts.lower_bound(dp.first), 
                      sharedEnts.upper_bound(dp.second));
  }
  else shared_ents = sharedEnts;

    // filter by iface
  if (iface) {
    result = filter_pstatus(shared_ents, PSTATUS_INTERFACE, PSTATUS_AND);
    RRA("");
  }
  
    // filter by owned
  if (owned_filter) {
    result = filter_pstatus(shared_ents, PSTATUS_NOT_OWNED, PSTATUS_NOT);
    RRA("");
  }

    // filter by proc
  if (-1 != other_proc) {
    result = filter_pstatus(shared_ents, PSTATUS_SHARED, PSTATUS_AND, other_proc);
    RRA("");
  }
  
  return result;
}


} // namespace moab

#ifdef TEST_PARALLELCOMM

#include <iostream>

#include "moab/Core.hpp"
#include "moab/ParallelComm.hpp"
#include "moab/Range.hpp"

#define PM {std::cerr << "Test failed; error message:" << std::endl;\
          std::string errmsg; \
          dynamic_cast<Core*>(my_impl)->get_last_error(errmsg); \
          std::cerr << errmsg << std::endl;\
          return 1;}

using namespace moab;

int main(int argc, char* argv[])
{

    // Check command line arg
  if (argc < 2)
  {
    std::cout << "Usage: " << argv[0] << " <mesh_file_name>" << std::endl;
    exit(1);
  }

  const char* file = argv[1];
  Core *my_impl = new Core(0, 2);
  Interface* mbImpl = my_impl;

    // create a communicator class, which will start mpi too
  ParallelComm pcomm(mbImpl, my_impl->tag_server(), my_impl->sequence_manager());
  ErrorCode result;

    // load the mesh
  result = mbImpl->load_mesh(file, 0, 0);
  if (MB_SUCCESS != result) return result;

    // get the mesh
  Range all_mesh, whole_range;
  result = mbImpl->get_entities_by_dimension(0, 3, all_mesh);
  if (MB_SUCCESS != result) return result;
    
  int buff_size;
  result = pcomm.pack_buffer(all_mesh, false, true, true, false, whole_range, buff_size);
  PM;


    // allocate space in the buffer
  pcomm.buffer_size(buff_size);

    // pack the actual buffer
  int actual_buff_size;
  result = pcomm.pack_buffer(whole_range, false, true, false, false, all_mesh, 
                             actual_buff_size);
  PM;

    // list the entities that got packed
  std::cout << "ENTITIES PACKED:" << std::endl;
  mbImpl->list_entities(all_mesh);

    // get the buffer
  std::vector<unsigned char> tmp_buffer;
  pcomm.take_buffer(tmp_buffer);
    
    // stop and restart MOAB
  delete mbImpl;
  my_impl = new Core(1, 2);
  mbImpl = my_impl;
    
    // create a new communicator class, using our old buffer
  ParallelComm pcomm2(mbImpl, my_impl->tag_server(), my_impl->sequence_manager(),
                        tmp_buffer);

    // unpack the results
  all_mesh.clear();
  result = pcomm2.unpack_buffer(all_mesh, store_remote_handles, from_proc);
  PM;
  
  std::cout << "ENTITIES UNPACKED:" << std::endl;
  mbImpl->list_entities(all_mesh);
  
  std::cout << "Success, processor " << mbImpl->proc_rank() << "." << std::endl;
  
  return 1;
}

#endif
