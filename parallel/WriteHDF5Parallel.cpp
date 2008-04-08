
#undef DEBUG

#ifdef DEBUG
#  include <stdio.h>
#  include <stdarg.h>
#endif

#ifndef HDF5_FILE
#  error Attempt to compile WriteHDF5Parallel with HDF5 support disabled
#endif

#include <stdlib.h>
#include <string.h>

#include <vector>
#include <set>
#include <map>
#include <utility>

#include <mpi.h>

#include <H5Tpublic.h>
#include <H5Ppublic.h>
#include <H5FDmpi.h>
#include <H5FDmpio.h>

#include "mhdf.h"

#include "MBInterface.hpp"
#include "MBInternals.hpp"
#include "MBTagConventions.hpp"
#include "MBParallelConventions.h"
#include "MBCN.hpp"
#include "MBWriteUtilIface.hpp"
#include "MBRange.hpp"

#include "WriteHDF5Parallel.hpp"


#ifdef DEBUG
#  define START_SERIAL                     \
     for (int _x = 0; _x < numProc; ++_x) {\
       MPI_Barrier( MPI_COMM_WORLD );      \
       if (_x != myRank) continue     
#  define END_SERIAL                       \
     }                                     \
     MPI_Barrier( MPI_COMM_WORLD )
#else
#  define START_SERIAL
#  define END_SERIAL
#endif


#define DEBUG_OUT_STREAM stdout

#ifndef DEBUG
static void printdebug( const char*, ... ) {}
#else
static void printdebug( const char* fmt, ... )
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  fprintf( DEBUG_OUT_STREAM, "[%d] ", rank );
  va_list args;
  va_start( args, fmt );
  vfprintf( DEBUG_OUT_STREAM, fmt, args );
  va_end( args );
  fflush( DEBUG_OUT_STREAM );
}
#endif


#ifdef NDEBUG
#  define assert(A)
#else
#  define assert(A) if (!(A)) do_assert(__FILE__, __LINE__, #A)
   static void do_assert( const char* file, int line, const char* condstr )
   {
     int rank;
     MPI_Comm_rank( MPI_COMM_WORLD, &rank );
     fprintf( DEBUG_OUT_STREAM, "[%d] Assert(%s) failed at %s:%d\n", rank, condstr, file, line );
     fflush( DEBUG_OUT_STREAM );
     abort();
   }
#endif


#ifndef DEBUG
void WriteHDF5Parallel::printrange( MBRange& ) {}
#else
void WriteHDF5Parallel::printrange( MBRange& r )
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MBEntityType type = MBMAXTYPE;
  for (MBRange::const_pair_iterator i = r.const_pair_begin(); i != r.const_pair_end(); ++i)
  {
    MBEntityHandle a, b;
    a = (*i).first;
    b = (*i).second;
    MBEntityType mytype = iFace->type_from_handle(a);
    if (mytype != type)
    {
      type = mytype;
      fprintf(DEBUG_OUT_STREAM, "%s[%d]  %s", type == MBMAXTYPE ? "" : "\n", rank, MBCN::EntityTypeName( type ) );
    }
    unsigned long id1 = iFace->id_from_handle( a );
    unsigned long id2 = iFace->id_from_handle( b );
    if (id1 == id2)
      fprintf(DEBUG_OUT_STREAM, " %lu", id1 );
    else
      fprintf(DEBUG_OUT_STREAM, " %lu-%lu", id1, id2 );
  }
  fprintf(DEBUG_OUT_STREAM, "\n");
  fflush( DEBUG_OUT_STREAM );
}
#endif


#ifndef DEBUG
static void print_type_sets( MBInterface* , int , int , MBRange& ) {}
#else
static void print_type_sets( MBInterface* iFace, int myRank, int numProc, MBRange& sets )
{
  MBTag gid, did, bid, sid, nid, iid;
  iFace->tag_get_handle( GLOBAL_ID_TAG_NAME, gid ); 
  iFace->tag_get_handle( GEOM_DIMENSION_TAG_NAME, did );
  iFace->tag_get_handle( MATERIAL_SET_TAG_NAME, bid );
  iFace->tag_get_handle( DIRICHLET_SET_TAG_NAME, nid );
  iFace->tag_get_handle( NEUMANN_SET_TAG_NAME, sid );
  iFace->tag_get_handle( PARALLEL_PARTITION_TAG_NAME, iid );
  MBRange typesets[10];
  const char* typenames[] = {"Block", "Sideset", "NodeSet", "Vertex", "Curve", "Surface", "Volume", "Body", "Partition", "Other"};
  for (MBRange::iterator riter = sets.begin(); riter != sets.end(); ++riter)
  {
    unsigned dim, id, proc[2], oldsize;
    if (MB_SUCCESS == iFace->tag_get_data(bid, &*riter, 1, &id)) 
      dim = 0;
    else if (MB_SUCCESS == iFace->tag_get_data(sid, &*riter, 1, &id))
      dim = 1;
    else if (MB_SUCCESS == iFace->tag_get_data(nid, &*riter, 1, &id))
      dim = 2;
    else if (MB_SUCCESS == iFace->tag_get_data(did, &*riter, 1, &dim)) {
      id = 0;
      iFace->tag_get_data(gid, &*riter, 1, &id);
      dim += 3;
    }
    else if (MB_SUCCESS == iFace->tag_get_data(iid, &*riter, 1, proc)) {
      assert(proc[0] == (unsigned)myRank || proc[1] == (unsigned)myRank);
      id = proc[proc[0] == (unsigned)myRank];
      dim = 8;
    }
    else {
      id = *riter;
      dim = 9;
    }

    oldsize = typesets[dim].size();
    typesets[dim].insert( id );
    assert( typesets[dim].size() - oldsize == 1 );  
  }
  for (int ii = 0; ii < 10; ++ii)
  {
    char num[16];
    std::string line(typenames[ii]);
    if (typesets[ii].empty())
      continue;
    sprintf(num, "(%lu):",(unsigned long)typesets[ii].size());
    line += num;
    for (MBRange::const_pair_iterator piter = typesets[ii].const_pair_begin();
         piter != typesets[ii].const_pair_end(); ++piter)
    {
      sprintf(num," %lx", (unsigned long)(*piter).first);
      line += num;
      if ((*piter).first != (*piter).second) {
        sprintf(num,"-%lx", (unsigned long)(*piter).second);
        line += num;
      }
    }

    printdebug ("%s\n", line.c_str());
  }
  printdebug("Total: %lu\n", (unsigned long)sets.size());
}
#endif


void range_remove( MBRange& from, const MBRange& removed )
{
  
/* The following should be more efficient, but isn't due
   to the inefficient implementation of MBRange::erase(iter,iter)
  MBRange::const_iterator s, e, n = from.begin();
  for (MBRange::const_pair_iterator p = removed.pair_begin();
       p != removed.pair_end(); ++p)
  {
    e = s = MBRange::lower_bound(n, from.end(), (*p).first);
    e = MBRange::lower_bound(s, from.end(), (*p).second);
    if (e != from.end() && *e == (*p).second)
      ++e;
    n = from.erase( s, e );
  }
*/

  if (removed.size())
  {
    MBRange tmp = from.subtract(removed);
    from.swap( tmp );
  }
}

MBWriterIface* WriteHDF5Parallel::factory( MBInterface* iface )
  { return new WriteHDF5Parallel( iface ); }

void WriteHDF5Parallel::MultiProcSetTags::add( const std::string& name )
  { list.push_back( Data(name) ); }

void WriteHDF5Parallel::MultiProcSetTags::add( const std::string& filter, 
                                               const std::string& data )
  { list.push_back( Data(filter,data) ); }

void WriteHDF5Parallel::MultiProcSetTags::add( const std::string& filter, 
                                               int filterval,
                                               const std::string& data )
  { list.push_back( Data(filter,data,filterval) ); }


WriteHDF5Parallel::WriteHDF5Parallel( MBInterface* iface )
  : WriteHDF5(iface)
{
  multiProcSetTags.add(  MATERIAL_SET_TAG_NAME );
  multiProcSetTags.add( DIRICHLET_SET_TAG_NAME );
  multiProcSetTags.add(   NEUMANN_SET_TAG_NAME );
  multiProcSetTags.add( GEOM_DIMENSION_TAG_NAME, 0, GLOBAL_ID_TAG_NAME );
  multiProcSetTags.add( GEOM_DIMENSION_TAG_NAME, 1, GLOBAL_ID_TAG_NAME );
  multiProcSetTags.add( GEOM_DIMENSION_TAG_NAME, 2, GLOBAL_ID_TAG_NAME );
  multiProcSetTags.add( GEOM_DIMENSION_TAG_NAME, 3, GLOBAL_ID_TAG_NAME );
}

WriteHDF5Parallel::WriteHDF5Parallel( MBInterface* iface,
                                      const std::vector<std::string>& tag_names )
  : WriteHDF5(iface)
{
  for(std::vector<std::string>::const_iterator i = tag_names.begin();
      i != tag_names.end(); ++i)
    multiProcSetTags.add( *i );
}

WriteHDF5Parallel::WriteHDF5Parallel( MBInterface* iface,
                                      const MultiProcSetTags& set_tags )
  : WriteHDF5(iface), multiProcSetTags(set_tags)
{}

// The parent WriteHDF5 class has ExportSet structs that are
// populated with the entities to be written, grouped by type
// (and for elements, connectivity length).  This function:
//  o determines which entities are to be written by a remote processor
//  o removes those entities from the ExportSet structs in WriteMesh
//  o puts them in the 'remoteMesh' array of MBRanges in this class
//  o sets their file Id to '1'
MBErrorCode WriteHDF5Parallel::gather_interface_meshes()
{
  MBRange range;
  MBErrorCode result;
  MBTag iface_tag, geom_tag;
  int i, proc_pair[2];
  
  //START_SERIAL;
  printdebug( "Pre-interface mesh:\n");
  printrange(nodeSet.range);
  for (std::list<ExportSet>::iterator eiter = exportList.begin();
           eiter != exportList.end(); ++eiter )
    printrange(eiter->range);
  printrange(setSet.range);
  
    // Allocate space for remote mesh data
  remoteMesh.resize( numProc );
  
    // Get tag handles
  result = iFace->tag_get_handle( PARALLEL_SHARED_PROC_TAG_NAME, iface_tag );
  if (MB_SUCCESS != result) {
    iface_tag = 0;
    return MB_SUCCESS;
  }
  result = iFace->tag_get_handle( GEOM_DIMENSION_TAG_NAME, geom_tag );
  if (MB_SUCCESS != result) return result;
  
  
    // Get interface mesh sets
  result = iFace->get_entities_by_type_and_tag( 0,
                                                MBENTITYSET,
                                                &iface_tag,
                                                0,
                                                1,
                                                range );
  if (MB_SUCCESS != result) return result;
  
  
    // Populate lists of interface mesh entities
  for (MBRange::iterator iiter = range.begin(); iiter != range.end(); ++iiter)
  {
    result = iFace->tag_get_data( iface_tag, &*iiter, 1, proc_pair );
    if (MB_SUCCESS != result) return result;
    const int remote_proc = proc_pair[0];
    
      // Get list of all entities in interface and 
      // the subset of that list that are meshsets.
    MBRange entities, sets;
    result = iFace->get_entities_by_handle( *iiter, entities );
    if (MB_SUCCESS != result) return result;
    result = iFace->get_entities_by_type( *iiter, MBENTITYSET, sets );
    if (MB_SUCCESS != result) return result;

      // Put any non-meshset entities in the list directly.
    //range_remove( entities, sets ); //not necessary, get_entities_by_handle doesn't return sets
    remoteMesh[remote_proc].merge( entities );
    //remoteMesh[remote_proc].insert( *iiter );
    
    for (MBRange::iterator siter = sets.begin(); siter != sets.end(); ++siter)
    {
        // For current parallel meshing code, root processor owns
        // all curve and geometric vertex meshes.  
      int dimension;
      result = iFace->tag_get_data( geom_tag, &*siter, 1, &dimension );
      if (result == MB_SUCCESS && dimension < 2)
        continue;
        
        // Put entities in list for appropriate processor.
      //remoteMesh[remote_proc].insert( *siter );
      entities.clear();
      result = iFace->get_entities_by_handle( *siter, entities );
      if (MB_SUCCESS != result) return result;
      remoteMesh[remote_proc].merge( entities );
    }
  }
  
    // For current parallel meshing code, root processor owns
    // all curve and geometric vertex meshes.  Find them and
    // allocate them appropriately.
  MBRange curves_and_verts;
  MBTag tags[] = { geom_tag, geom_tag };
  int value_ints[] = { 0, 1 };
  const void* values[] = {value_ints, value_ints + 1};
  result = iFace->get_entities_by_type_and_tag( 0, MBENTITYSET,
                                                tags, values, 2,
                                                curves_and_verts, 
                                                MBInterface::UNION );
                                                assert(MB_SUCCESS == result);
  MBRange edges, nodes;
  for (MBRange::iterator riter = curves_and_verts.begin();
       riter != curves_and_verts.end(); ++riter)
  {
    result = iFace->get_entities_by_type( *riter, MBVERTEX, nodes ); assert(MB_SUCCESS == result);
    result = iFace->get_entities_by_type( *riter, MBEDGE, edges ); assert(MB_SUCCESS == result);
  }
  std::list<ExportSet>::iterator eiter = exportList.begin();
  for ( ; eiter != exportList.end() && eiter->type != MBEDGE; ++eiter );
  
  remoteMesh[0].merge( nodes );
  remoteMesh[0].merge( edges );
  //remoteMesh[0].merge( curves_and_verts );
  if (myRank == 0)
  {
    nodeSet.range.merge( nodes );
    //setSet.range.merge(curves_and_verts);
    eiter->range.merge( edges );
  } 
  edges.merge(nodes);
  //edges.merge(curves_and_verts);
  for (i = 1; i < numProc; i++)
  {
    MBRange diff = edges.intersect( remoteMesh[i] );
    range_remove(remoteMesh[i], diff);
  }
  
  
  
    // For all remote mesh entities, remove them from the
    // lists of local mesh to be exported and the ID map
    // (they will be put back into the ID map with different
    //  IDs later.)
  for (i = 0; i < numProc; i++)
  {
    if (i == myRank) continue;
    
    const MBRange& range = remoteMesh[i];
    
    range_remove( nodeSet.range, range );
    //range_remove( setSet.range, range );
    for (std::list<ExportSet>::iterator eiter = exportList.begin();
         eiter != exportList.end(); ++eiter )
      range_remove( eiter->range, range );
    
    //for (MBRange::const_pair_iterator pi = range.const_pair_begin();
    //     pi != range.const_pair_end(); ++pi) 
    //  idMap.erase( pi->first, pi->second - pi->first + 1 );
  }
  
    // print some debug output summarizing what we've accomplished
  
  printdebug("Remote mesh:\n");
  for (int ii = 0; ii < numProc; ++ii)
  {
    printdebug("  proc %d : %d\n", ii, remoteMesh[ii].size());
    printrange( remoteMesh[ii] );
  }

  printdebug( "Post-interface mesh:\n");
  printrange(nodeSet.range);
  for (std::list<ExportSet>::iterator eiter = exportList.begin();
           eiter != exportList.end(); ++eiter )
    printrange(eiter->range);
  printrange(setSet.range);

  //END_SERIAL;
  
  return MB_SUCCESS;
}



MBErrorCode WriteHDF5Parallel::create_file( const char* filename,
                                            bool overwrite,
                                            std::vector<std::string>& qa_records,
                                            int dimension,
                                            bool parallel )
{
  if (!parallel)
    return WriteHDF5::create_file(filename, overwrite, qa_records, dimension, false );

  MBErrorCode rval;
  int result;
  mhdf_Status status;
    
  result = MPI_Comm_rank( MPI_COMM_WORLD, &myRank );
  assert(MPI_SUCCESS == result);
  result = MPI_Comm_size( MPI_COMM_WORLD, &numProc );
  assert(MPI_SUCCESS == result);
  
  rval = gather_interface_meshes();
  if (MB_SUCCESS != rval) return rval;
  
    /**************** Create actual file and write meta info ***************/

  if (myRank == 0)
  {
      // create the file
    const char* type_names[MBMAXTYPE];
    memset( type_names, 0, MBMAXTYPE * sizeof(char*) );
    for (MBEntityType i = MBEDGE; i < MBENTITYSET; ++i)
      type_names[i] = MBCN::EntityTypeName( i );
   
    filePtr = mhdf_createFile( filename, overwrite, type_names, MBMAXTYPE, &status );
    if (!filePtr)
    {
      writeUtil->report_error( "%s\n", mhdf_message( &status ) );
      return MB_FAILURE;
    }
    
    rval = write_qa( qa_records );
    if (MB_SUCCESS != rval) return rval;
  }
  
  
     /**************** Create node coordinate table ***************/
 
  rval = create_node_table( dimension );
  if (MB_SUCCESS != rval) return rval;
  rval = communicate_remote_ids( MBVERTEX );
  if (MB_SUCCESS != rval) return rval;
  
    /**************** Create element tables ***************/

  rval = negotiate_type_list();
  if (MB_SUCCESS != rval) return rval;
  rval = create_element_tables();
  if (MB_SUCCESS != rval) return rval;
  for (std::list<ExportSet>::iterator ex_itor = exportList.begin(); 
       ex_itor != exportList.end(); ++ex_itor)
  {
    rval = communicate_remote_ids( ex_itor->type );
    assert(MB_SUCCESS == rval);
  }
  
  
    /**************** Create adjacency tables *********************/
  
  rval = create_adjacency_tables();
  if (MB_SUCCESS != rval) return rval;
  
    /**************** Create meshset tables *********************/
  
  rval = create_meshset_tables();
  if (MB_SUCCESS != rval) return rval;
  
  
    /* Need to write tags for shared sets this proc is responsible for */
  
  MBRange parallel_sets;
  for (std::list<ParallelSet>::const_iterator psiter = parallelSets.begin();
       psiter != parallelSets.end(); ++psiter)
    if (psiter->description)
      parallel_sets.insert( psiter->handle );
  
  setSet.range.merge( parallel_sets );
  rval = gather_tags();
  if (MB_SUCCESS != rval)
    return rval;
  range_remove( setSet.range, parallel_sets );   
  

    /**************** Create tag data *********************/
  
  std::list<SparseTag>::iterator tag_iter;
  sort_tags_by_name();
  const int num_tags = tagList.size();
  
    // Construct vector (tag_counts) containing a pair of values for each
    // tag, where the first value in the pair is the number of entities on
    // this processor for which the tag has been set.  The second value is
    // zero for normal tags.  For variable-length tags it is the total number
    // of tag values set for all entities on this processor.
  std::vector<unsigned long> tag_offsets(2*num_tags), tag_counts(2*num_tags);
  std::vector<unsigned long>::iterator tag_off_iter = tag_counts.begin();
  for (tag_iter = tagList.begin(); tag_iter != tagList.end(); ++tag_iter) {
    int s;
    *tag_off_iter = tag_iter->range.size();
    ++tag_off_iter;
    if (MB_VARIABLE_DATA_LENGTH == iFace->tag_get_size( tag_iter->tag_id, s )) {
      unsigned long total_len;
      rval = get_tag_data_length( *tag_iter, total_len );
      if (MB_SUCCESS != rval)
        return rval;
      
      *tag_off_iter = total_len;
      assert(total_len == *tag_off_iter);
    }
    else {
      *tag_off_iter = 0;
    }
    ++tag_off_iter;
  }
  
    // Populate proc_tag_offsets on root processor with the values from
    // tag_counts on each processor.
  printdebug("Exchanging tag data for %d tags.\n", num_tags);
  std::vector<unsigned long> proc_tag_offsets(2*num_tags*numProc);
  result = MPI_Gather( &tag_counts[0], 2*num_tags, MPI_UNSIGNED_LONG,
                 &proc_tag_offsets[0], 2*num_tags, MPI_UNSIGNED_LONG,
                       0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
    // Calculate the total counts over all processors (tag_counts)
    // and the offset at which each processor should begin writing
    // its data (proc_tag_offsets).  Both lists contain a pair of
    // values, where the first is the number of entities and the second
    // is either unused for fixed-length tags or the total data table
    // size for variable-length tags.
  tag_iter = tagList.begin();
  for (int i = 0; i < num_tags; ++i, ++tag_iter)
  {
    tag_counts[2*i] = tag_counts[2*i+1] = 0;
    unsigned long next_offset = 0;
    unsigned long next_var_len_offset = 0;
    for (int j = 0; j < numProc; j++)
    {
      unsigned long count = proc_tag_offsets[2*i + j*2*num_tags];
      proc_tag_offsets[2*i + j*2*num_tags] = next_offset;
      next_offset += count;
      tag_counts[2*i] += count;
      
      count = proc_tag_offsets[2*i + 1 + j*2*num_tags];
      proc_tag_offsets[2*i + 1 + j*2*num_tags] = next_var_len_offset;
      next_var_len_offset += count;
      tag_counts[2*i + 1] += count;
    }

    if (0 == myRank)
    {
      rval = create_tag(tag_iter->tag_id, next_offset, next_var_len_offset);
      assert(MB_SUCCESS == rval);
      printdebug( "Creating table of size %lu for tag 0x%lx\n", 
                  next_var_len_offset ? next_var_len_offset : next_offset, 
                  (unsigned long)tag_iter->tag_id );
    }
  }
  
    // Send total counts to all processors.  This is necessary because all 
    // processors need to know if we are not writing anything for the tag (count == 0).  
  result = MPI_Bcast( &tag_counts[0], 2*num_tags, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
    // Send to each processor its per-tag offset values.
  result = MPI_Scatter( &proc_tag_offsets[0], 2*num_tags, MPI_UNSIGNED_LONG,
                             &tag_offsets[0], 2*num_tags, MPI_UNSIGNED_LONG,
                             0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);


  tag_iter = tagList.begin();
  for (int i = 0; i < num_tags; ++i, ++tag_iter)
  {
    tag_iter->offset = tag_offsets[2*i];
    tag_iter->write = tag_counts[2*i] > 0;
    tag_iter->varDataOffset = tag_offsets[2*i + 1];
  }

  #ifdef DEBUG
  START_SERIAL;  
  printdebug("Tags: %12s %8s %8s %8s %8s %8s\n", "Name", "Count", "Offset", "Var Off", "Var Len", "Handle");

  tag_iter = tagList.begin();
  for (int i = 0; i < num_tags; ++i, ++tag_iter)
  {
    std::string name;
    iFace->tag_get_name( tag_iter->tag_id, name );
    printdebug("%18s %8lu %8lu %8lu %8lu 0x%7lx\n", name.c_str(), tag_counts[2*i], tag_offsets[2*i], tag_offsets[2*i+1], tag_counts[2*i+1], (unsigned long)tag_iter->tag_id );
  }
  END_SERIAL;  
  #endif
  
  /************** Close serial file and reopen parallel *****************/
  
  if (0 == myRank)
  {
    mhdf_closeFile( filePtr, &status );
  }
  
  unsigned long junk;
  hid_t hdf_opt = H5Pcreate( H5P_FILE_ACCESS );
  H5Pset_fapl_mpio( hdf_opt, MPI_COMM_WORLD, MPI_INFO_NULL );
  filePtr = mhdf_openFileWithOpt( filename, 1, &junk, hdf_opt, &status );
  if (!filePtr)
  {
    writeUtil->report_error( "%s\n", mhdf_message( &status ) );
    return MB_FAILURE;
  }
  
  
  return MB_SUCCESS;
}


MBErrorCode WriteHDF5Parallel::create_node_table( int dimension )
{
  int result;
  mhdf_Status status;
 
    // gather node counts for each processor
  std::vector<long> node_counts(numProc);
  long num_nodes = nodeSet.range.size();
  result = MPI_Gather( &num_nodes, 1, MPI_LONG, &node_counts[0], 1, MPI_LONG, 0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
    // create node data in file
  long first_id;
  if (myRank == 0)
  {
    int total = 0;
    for (int i = 0; i < numProc; i++)
      total += node_counts[i];
      
    hid_t handle = mhdf_createNodeCoords( filePtr, dimension, total, &first_id, &status );
    if (mhdf_isError( &status ))
    {
      writeUtil->report_error( "%s\n", mhdf_message( &status ) );
      return MB_FAILURE;
    }
    mhdf_closeData( filePtr, handle, &status );
 }
    
    // send id offset to every proc
  result = MPI_Bcast( &first_id, 1, MPI_LONG, 0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  nodeSet.first_id = (id_t)first_id;
   
      // calculate per-processor offsets
  if (myRank == 0)
  {
    int prev_size = node_counts[0];
    node_counts[0] = 0;
    for (int i = 1; i < numProc; ++i)
    {
      int mysize = node_counts[i];
      node_counts[i] = node_counts[i-1] + prev_size;
      prev_size = mysize;
    }
  }
  
    // send each proc it's offset in the node table
  int offset;
  result = MPI_Scatter( &node_counts[0], 1, MPI_INT, 
                        &offset, 1, MPI_INT,
                        0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  nodeSet.offset = offset;
  
  return assign_ids( nodeSet.range, nodeSet.first_id + nodeSet.offset );
}



struct elemtype {
  int mbtype;
  int numnode;
  
  elemtype( int vals[2] ) : mbtype(vals[0]), numnode(vals[1]) {}
  elemtype( int t, int n ) : mbtype(t), numnode(n) {}
  
  bool operator==( const elemtype& other ) const
  {
    return mbtype == other.mbtype &&
            (mbtype == MBENTITYSET ||
             numnode == other.numnode);
  }
  bool operator<( const elemtype& other ) const
  {
    if (mbtype > other.mbtype)
      return false;
   
    return mbtype < other.mbtype ||
           (mbtype != MBENTITYSET &&
            numnode < other.numnode);
  }
  bool operator!=( const elemtype& other ) const
    { return !this->operator==(other); }
};


MBErrorCode WriteHDF5Parallel::negotiate_type_list()
{
  int result;
  
  exportList.sort();
  
    // Get number of types each processor has
  int num_types = 2*exportList.size();
  std::vector<int> counts(numProc);
  result = MPI_Gather( &num_types, 1, MPI_INT, &counts[0], 1, MPI_INT, 0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
    // Get list of types on this processor
  std::vector<int> my_types(num_types);
  std::vector<int>::iterator viter = my_types.begin();
  for (std::list<ExportSet>::iterator eiter = exportList.begin();
       eiter != exportList.end(); ++eiter)
  {
    *viter = eiter->type;      ++viter;
    *viter = eiter->num_nodes; ++viter;
  }

  #ifdef DEBUG
  START_SERIAL;
  printdebug( "Local Element Types:\n");
  viter = my_types.begin();
  while (viter != my_types.end())
  {
    int type = *viter; ++viter;
    int count = *viter; ++viter;
    printdebug("  %s : %d\n", MBCN::EntityTypeName((MBEntityType)type), count);
  }
  END_SERIAL;
  #endif

    // Get list of types from each processor
  std::vector<int> displs(numProc + 1);
  displs[0] = 0;
  for (long i = 1; i <= numProc; ++i)
    displs[i] = displs[i-1] + counts[i-1];
  int total = displs[numProc];
  std::vector<int> alltypes(total);
  result = MPI_Gatherv( &my_types[0], my_types.size(), MPI_INT,
                        &alltypes[0], &counts[0], &displs[0], MPI_INT,
                        0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
    // Merge type lists
  std::list<elemtype> type_list;
  std::list<elemtype>::iterator liter;
  for (int i = 0; i < numProc; ++i)
  {
    int* proc_type_list = &alltypes[displs[i]];
    liter = type_list.begin();
    for (int j = 0; j < counts[i]; j += 2)
    {
      elemtype type( &proc_type_list[j] );
        // skip until insertion spot
      for (; liter != type_list.end() && *liter < type; ++liter);
      
      if (liter == type_list.end() || *liter != type)
        liter = type_list.insert( liter, type );
    }
  }
  
    // Send total number of types to each processor
  total = type_list.size();
  result = MPI_Bcast( &total, 1, MPI_INT, 0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
    // Send list of types to each processor
  std::vector<int> intlist(total * 2);
  viter = intlist.begin();
  for (liter = type_list.begin(); liter != type_list.end(); ++liter)
  {
    *viter = liter->mbtype;  ++viter;
    *viter = liter->numnode; ++viter;
  }
  result = MPI_Bcast( &intlist[0], 2*total, MPI_INT, 0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);

  #ifdef DEBUG
  START_SERIAL;
  printdebug( "Global Element Types:\n");
  viter = intlist.begin();
  while (viter != intlist.end())
  {
    int type = *viter; ++viter;
    int count = *viter; ++viter;
    printdebug("  %s : %d\n", MBCN::EntityTypeName((MBEntityType)type), count);
  }
  END_SERIAL;
  #endif
  
    // Insert missing types into exportList, with an empty
    // range of entities to export.
  std::list<ExportSet>::iterator ex_iter = exportList.begin();
  viter = intlist.begin();
  for (long i = 0; i < total; ++i)
  {
    long mbtype = *viter; ++viter;
    long numnode = *viter; ++viter;
    while (ex_iter != exportList.end() && ex_iter->type < mbtype)
      ++ex_iter;
    
    bool equal = ex_iter != exportList.end() && ex_iter->type == mbtype;
    if (equal)
    {
      while (ex_iter != exportList.end() && ex_iter->num_nodes < numnode)
        ++ex_iter;
        
      equal = ex_iter != exportList.end() && ex_iter->num_nodes == numnode;
    }
    
    if (!equal)
    {
      ExportSet insert;
      insert.type = (MBEntityType)mbtype;
      insert.num_nodes = numnode;
      insert.first_id = 0;
      insert.offset = 0;
      insert.adj_offset = 0;
      ex_iter = exportList.insert( ex_iter, insert );
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode WriteHDF5Parallel::create_element_tables()
{
  int result;
  MBErrorCode rval;
  std::list<ExportSet>::iterator ex_iter;
  std::vector<long>::iterator viter;
  
    // Get number of each element type from each processor
  const int numtypes = exportList.size();
  std::vector<long> my_counts(numtypes);
  std::vector<long> counts(numtypes * numProc + numtypes);
  viter = my_counts.begin();
  for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter)
    { *viter = ex_iter->range.size(); ++viter; }
  
  result = MPI_Gather( &my_counts[0], numtypes, MPI_LONG,
                       &counts[0],    numtypes, MPI_LONG, 0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
    // Convert counts to offsets
  for (int i = 0; i < numtypes; i++) 
  {
    long prev = 0;
    for (int j = 0; j <= numProc; j++)
    {
      long tmp = counts[j*numtypes + i];
      counts[j*numtypes+i] = prev;
      prev += tmp;
    }
  }
  
    // Send offsets to each processor
  result = MPI_Scatter( &counts[0],    numtypes, MPI_LONG,
                        &my_counts[0], numtypes, MPI_LONG,
                        0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
    // Update store offsets in ExportSets
  viter = my_counts.begin();
  for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter)
    ex_iter->offset = (id_t)*(viter++);
  
    // Create element tables
  std::vector<long> start_ids(numtypes);
  if (myRank == 0)
  {
    viter = start_ids.begin();
    long* citer = &counts[numtypes * numProc];
    for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter)
    {
      rval = create_elem_tables( ex_iter->type,
                                 ex_iter->num_nodes,
                                 *citer,
                                 *viter );
      assert(MB_SUCCESS == rval);
      ++citer;
      ++viter;
    }
  }
  
    // send start IDs to each processor
  result = MPI_Bcast( &start_ids[0], numtypes, MPI_LONG, 0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
    // Assign IDs to local elements
  viter = start_ids.begin();
  for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter)
  {
    ex_iter->first_id = *(viter++);
    id_t myfirst = (id_t)(ex_iter->first_id + ex_iter->offset);
    rval = assign_ids( ex_iter->range, myfirst );
    assert(MB_SUCCESS == rval);
  }
  
  return MB_SUCCESS;
}
  
MBErrorCode WriteHDF5Parallel::create_adjacency_tables()
{
  MBErrorCode rval;
  mhdf_Status status;
  int i, j, result;
#ifdef WRITE_NODE_ADJACENCIES  
  const int numtypes = exportList.size()+1;
#else
  const int numtypes = exportList.size();
#endif
  std::vector<long>::iterator viter;
  std::list<ExportSet>::iterator ex_iter;
  std::vector<long> local(numtypes), all(numProc * numtypes + numtypes);
  
    // Get adjacency counts for local processor
  viter = local.begin();
  id_t num_adj;
#ifdef WRITE_NODE_ADJACENCIES  
  rval = count_adjacencies( nodeSet.range, num_adj );
  assert (MB_SUCCESS == rval);
  *viter = num_adj; ++viter;
#endif

  for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter)
  {
    rval = count_adjacencies( ex_iter->range, num_adj );
    assert (MB_SUCCESS == rval);
    *viter = num_adj; ++viter;
  }
  
    // Send local adjacency counts to root processor
  result = MPI_Gather( &local[0], numtypes, MPI_LONG,
                       &all[0],   numtypes, MPI_LONG, 
                       0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
    // Convert counts to offsets
  for (i = 0; i < numtypes; i++) 
  {
    long prev = 0;
    for (j = 0; j <= numProc; j++)
    {
      long tmp = all[j*numtypes + i];
      all[j*numtypes+i] = prev;
      prev += tmp;
    }
  }
  
    // For each element type for which there is no adjacency data,
    // send -1 to all processors as the offset
  for (i = 0; i < numtypes; ++i)
    if (all[numtypes*numProc+i] == 0)
      for (j = 0; j < numProc; ++j)
        all[j*numtypes+i] = -1;
  
    // Send offsets back to each processor
  result = MPI_Scatter( &all[0],   numtypes, MPI_LONG,
                        &local[0], numtypes, MPI_LONG,
                        0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
    // Record the adjacency offset in each ExportSet
  viter = local.begin();
#ifdef WRITE_NODE_ADJACENCIES  
  nodeSet.adj_offset = *viter; ++viter;
#endif
  for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter)
    { ex_iter->adj_offset = *viter; ++viter; }
  
    // Create data tables in file
  if (myRank == 0)
  {
    viter = all.begin() + (numtypes * numProc);
#ifdef WRITE_NODE_ADJACENCIES  
    if (*viter) {
      hid_t handle = mhdf_createAdjacency( filePtr, 
                                           mhdf_node_type_handle(),
                                           *viter,
                                           &status );
      if (mhdf_isError( &status ))
      {
        writeUtil->report_error( "%s\n", mhdf_message( &status ) );
        return MB_FAILURE;
      }
      mhdf_closeData( filePtr, handle, &status );
    }
    ++viter;
#endif
    for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter, ++viter)
    {
      if (!*viter) 
        continue;
      
      hid_t handle = mhdf_createAdjacency( filePtr,
                                           ex_iter->name(),
                                           *viter,
                                           &status );
      if (mhdf_isError( &status ))
      {
        writeUtil->report_error( "%s\n", mhdf_message( &status ) );
        return MB_FAILURE;
      }
      mhdf_closeData( filePtr, handle, &status );
    }
  }

  return MB_SUCCESS;
}

/*
MBErrorCode WriteHDF5Parallel::get_interface_set_data( RemoteSetData& data,
                                                       long& offset )
{
  const char* PROC_ID_TAG = "HDF5Writer_Rank";
  MBTag iface_tag, proc_tag;
  MBErrorCode rval;
  
  rval = iFace->tag_get_handle( PARALLEL_INTERFACE_TAG_NAME, iface_tag );
  if (MB_SUCCESS != rval) return rval;
  
  rval = iFace->tag_get_handle( PROC_ID_TAG, proc_tag );
  if (MB_SUCCESS == rval) 
    iFace->tag_delete( proc_tag );
  rval = iFace->tag_create( PROC_ID_TAG, sizeof(int), MB_TAG_DENSE, MB_TYPE_INTEGER, proc_tag, 0 );
  if (MB_SUCCESS != rval) return rval;
    
  MBRange interface_sets, sets;
  rval = iFace->get_entities_by_type_and_tag( 0, MBENTITYSET, &iface_tag, 0, 1, interface_sets );
  if (MB_SUCCESS != rval) return rval;
  
  std::vector<int> list;
  for (MBRange::iterator i = interface_sets.begin(); i != interface_sets.end(); ++i)
  {
    int proc_ids[2];
    rval = iFace->tag_get_data( iface_tag, &*i, 1, proc_ids );
    if (MB_SUCCESS != rval) return rval;
    
    sets.clear();
    rval = iFace->get_entities_by_type( *i, MBENTITYSET, sets );
    if (MB_SUCCESS != rval) return rval;
  
    list.clear();
    list.resize( sets.size(), proc_ids[0] );
    rval = iFace->tag_set_data( proc_tag, sets, &list[0] );
    if (MB_SUCCESS != rval) return rval;
  }
  
  return get_remote_set_data( PROC_ID_TAG, PARALLEL_GLOBAL_ID_TAG_NAME, data, offset );
}
*/
  
/** Working data for group of global sets identified by an ID tag */
struct RemoteSetData {
  MBTag data_tag;    //!< The ID tag for matching sets across processors
  MBTag filter_tag;  //!< Optional tag to filter on (e.g. geometric dimension)
  int filter_value;  //!< Value of filter_tag for this group
  MBRange range;     //!< Set handles with data_tag set (and optionally filter_tag == filter_value)
  std::vector<int> counts;       //!< Number of sets with tag on each proc, indexed by MPI rank
  std::vector<int> displs;       //!< Offset in all_values at which the data_tag values for
                                 //!< each processor begin. displs[n] = sum(i from 0 to n-1)(counts[i])
  std::vector<int> all_values;   //!< data_tag values for sets on all processors, 
                                 //!< counts[0] values for proc 0, then counts[1] values for proc 
                                 //!< 1, etc.
  std::vector<int> local_values; //!< data_tag values for sets that exist on this processor
};

MBErrorCode WriteHDF5Parallel::get_remote_set_data( 
                        const WriteHDF5Parallel::MultiProcSetTags::Data& tags,
                        RemoteSetData& data, long& offset )
{
  MBErrorCode rval;
  int i, result;
  MBRange::iterator riter;

  rval = iFace->tag_get_handle( tags.filterTag.c_str(), data.filter_tag );
  if (rval != MB_SUCCESS) return rval;
  if (tags.useFilterValue) 
  {
    i = 0;
    iFace->tag_get_size( data.filter_tag, i );
    if (i != sizeof(int)) {
      fprintf(stderr, "Cannot use non-int tag data for filtering remote sets.\n" );
      assert(0);
      return MB_FAILURE;
    }  
    data.filter_value = tags.filterValue;
  }
  else
  {
    data.filter_value = 0;
  }
  
  rval = iFace->tag_get_handle( tags.dataTag.c_str(), data.data_tag );
  if (rval != MB_SUCCESS) return rval;
  i = 0;
  iFace->tag_get_size( data.data_tag, i );
  if (i != sizeof(int)) {
    fprintf(stderr, "Cannot use non-int tag data for matching remote sets.\n" );
    assert(0);
    return MB_FAILURE;
  }  
    

  printdebug("Negotiating multi-proc meshsets for tag: \"%s\"\n", tags.filterTag.c_str());

    // Get sets with tag, or leave range empty if the tag
    // isn't defined on this processor.
  if (rval != MB_TAG_NOT_FOUND)
  {
    MBTag handles[] = { data.filter_tag, data.data_tag };
    const void* values[] = { tags.useFilterValue ? &tags.filterValue : 0, 0 };
    rval = iFace->get_entities_by_type_and_tag( 0, 
                                                MBENTITYSET, 
                                                handles,
                                                values,
                                                2,
                                                data.range );
    if (rval != MB_SUCCESS) return rval;
    MBRange tmp = data.range.intersect( setSet.range );
    data.range.swap( tmp );
    range_remove( setSet.range, data.range );
  }
  
  printdebug("Found %d meshsets with \"%s\" tag.\n", data.range.size(), tags.filterTag.c_str() );

    // Exchange number of sets with tag between all processors
  data.counts.resize(numProc);
  int count = data.range.size();
  result = MPI_Allgather( &count,          1, MPI_INT, 
                          &data.counts[0], 1, MPI_INT,
                          MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);

    // Exchange tag values for sets between all processors
  data.displs.resize(numProc+1);
  data.displs[0] = 0;
  for (i = 1; i <= numProc; i++)
    data.displs[i] = data.displs[i-1] + data.counts[i-1];
  int total = data.displs[numProc];
  data.all_values.resize(total);
  data.local_values.resize(count);
  rval = iFace->tag_get_data( data.data_tag, data.range, &data.local_values[0] );
  assert( MB_SUCCESS == rval );
  result = MPI_Allgatherv( &data.local_values[0], count, MPI_INT,
                           &data.all_values[0], &data.counts[0], &data.displs[0], MPI_INT,
                           MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);


    // Remove from the list any sets that are unique to one processor
  std::vector<int> sorted( data.all_values );
  std::sort( sorted.begin(), sorted.end() );
  int r = 0, w = 0;
  for (i = 0; i < numProc; ++i)
  {
    const int start = w;
    for (int j = 0; j < data.counts[i]; ++j)
    {
      std::vector<int>::iterator p 
        = std::lower_bound( sorted.begin(), sorted.end(), data.all_values[r] );
      ++p;
      if (p != sorted.end() && *p == data.all_values[r])
      {
        data.all_values[w] = data.all_values[r];
        ++w;
      }
      ++r;
    }
    data.counts[i] = w - start;
  }
  total = w;
  data.all_values.resize( total );
  r = w = 0;
  for (i = 0; i < count; ++i)
  {
    std::vector<int>::iterator p 
      = std::lower_bound( sorted.begin(), sorted.end(), data.local_values[r] );
    ++p;
    if (p != sorted.end() && *p == data.local_values[r])
    {
      data.local_values[w] = data.local_values[r];
      ++w;
    }
    else
    {
      riter = data.range.begin();
      riter += w;
      setSet.range.insert( *riter );
      data.range.erase( riter );
    }
    ++r;
  }
  count = data.range.size();
  assert( count == data.counts[myRank] );
  assert( count == w );
  data.local_values.resize( count );
  sorted.clear(); // release storage
  
    // recalculate displacements for updated counts
  data.displs[0] = 0;
  for (i = 1; i <= numProc; i++)
    data.displs[i] = data.displs[i-1] + data.counts[i-1];
  
    // Find sets that span multple processors and update appropriately.
    // The first processor (sorted by MPI rank) that contains a given set
    // will be responsible for writing the set description.  All multi-
    // processor sets will be written at the beginning of the set tables.
    // Processors will write set contents/children for a given set in
    // the order of their MPI rank.
    //
    // Identify which meshsets will be managed by this processor and
    // the corresponding offset in the set description table. 
  std::map<int,int> val_id_map; // Map from tag value to file ID for set
  int cpu = 0;
  for (i = 0; i < total; ++i)
  {
    if (data.displs[cpu+1] == i)
      ++cpu;

    int id = 0;
    std::map<int,int>::iterator p = val_id_map.find( data.all_values[i] );
    if (p == val_id_map.end())
    {
      id = (int)++offset;
      val_id_map[data.all_values[i]] = id;
      //const unsigned int values_offset = (unsigned)i - (unsigned)data.displs[myRank];
      //if (values_offset < (unsigned)count)
      //{
      //  riter = data.range.begin();
      //  riter += values_offset;
      //  myParallelSets.insert( *riter );
      //}
    }
    std::vector<int>::iterator loc 
      = std::find( data.local_values.begin(), data.local_values.end(), data.all_values[i] );
    if (loc != data.local_values.end()) 
    {
      riter = data.range.begin();
      riter += loc - data.local_values.begin();
      cpuParallelSets[cpu].insert( *riter );
    }
  }
  riter = data.range.begin();
  for (i = 0; i < count; ++i, ++riter)
  {
    std::map<int,int>::iterator p = val_id_map.find( data.local_values[i] );
    assert( p != val_id_map.end() );
    int id = p->second;
    if (idMap.end() == idMap.insert( *riter, id, 1 )) {
      assert(false);
      return MB_FAILURE;
    }
  }
  
  return MB_SUCCESS;
}


MBErrorCode WriteHDF5Parallel::create_meshset_tables()
{
  MBErrorCode rval;
  int result, i;
  long total_offset = 0;
  MBRange::const_iterator riter;

  START_SERIAL;
  print_type_sets( iFace, myRank, numProc, setSet.range );
  END_SERIAL;

    // Gather data about multi-processor meshsets - removes sets from setSet.range
  cpuParallelSets.resize( numProc );
  std::vector<RemoteSetData> remote_set_data( multiProcSetTags.list.size() );
  for (i = 0; i< (int)multiProcSetTags.list.size(); i++)
  {
    rval = get_remote_set_data( multiProcSetTags.list[i],
                                remote_set_data[i],
                                total_offset ); assert(MB_SUCCESS == rval);
  }
  //rval = get_interface_set_data( remote_set_data[i], total_offset );
  if (MB_SUCCESS != rval) return rval;

  START_SERIAL;
  printdebug("myLocalSets\n");
  print_type_sets( iFace, myRank, numProc, setSet.range );
  END_SERIAL;

    // Gather counts of non-shared sets from each proc
    // to determine total table size.
  std::vector<long> set_offsets(numProc + 1);
  long local_count = setSet.range.size();
  result = MPI_Gather( &local_count,    1, MPI_LONG,
                       &set_offsets[0], 1, MPI_LONG,
                       0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  for (i = 0; i <= numProc; i++)
  {
    long tmp = set_offsets[i];
    set_offsets[i] = total_offset;
    total_offset += tmp;
  }
  
    // Send each proc its offsets in the set description table.
  long sets_offset;
  result = MPI_Scatter( &set_offsets[0], 1, MPI_LONG,
                        &sets_offset,    1, MPI_LONG, 0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  setSet.offset = (id_t)(sets_offset);

    // Create the set description table
  long total_count_and_start_id[2] = { set_offsets[numProc], 0 };
  if (myRank == 0 && total_count_and_start_id[0] > 0)
  {
    rval = create_set_meta( (id_t)total_count_and_start_id[0], total_count_and_start_id[1] );
    assert (MB_SUCCESS == rval);
  }
  
    // Send totals to all procs.
  result = MPI_Bcast( total_count_and_start_id, 2, MPI_LONG, 0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  setSet.first_id = total_count_and_start_id[1];
  writeSets = total_count_and_start_id[0] > 0;

  START_SERIAL;  
  printdebug("Non-shared sets: %ld local, %ld global, offset = %ld, first_id = %ld\n",
    local_count, total_count_and_start_id[0], sets_offset, total_count_and_start_id[1] );
  printdebug("my Parallel Sets:\n");
  print_type_sets(iFace, myRank, numProc, cpuParallelSets[myRank] );
  END_SERIAL;
  
    // Not writing any sets??
  if (!writeSets)
    return MB_SUCCESS;
  
    // Assign set IDs
  assign_ids( setSet.range, setSet.first_id + setSet.offset );
  for (i = 0; i < (int)remote_set_data.size(); ++i)
    fix_remote_set_ids( remote_set_data[i], setSet.first_id );
  
    // Communicate sizes for remote sets
  long data_offsets[3] = { 0, 0, 0 };
  for (i = 0; i < (int)remote_set_data.size(); ++i)
  {
    rval = negotiate_remote_set_contents( remote_set_data[i], data_offsets ); 
    assert(MB_SUCCESS == rval);
  }
  remote_set_data.clear();
  
    // Exchange IDs for remote/adjacent sets not shared between procs
  //rval = communicate_remote_ids( MBENTITYSET ); assert(MB_SUCCESS == rval);
  
    // Communicate counts for local sets
  long data_counts[3];
  rval = count_set_size( setSet.range, rangeSets, data_counts[0], data_counts[1], data_counts[2] );
  if (MB_SUCCESS != rval) return rval;
  std::vector<long> set_counts(3*numProc);
  result = MPI_Gather( data_counts,    3, MPI_LONG,
                       &set_counts[0], 3, MPI_LONG,
                       0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  for (i = 0; i < 3*numProc; ++i)
  {
    long tmp = set_counts[i];
    set_counts[i] = data_offsets[i%3];
    data_offsets[i%3] += tmp;
  }
  long all_counts[] = {data_offsets[0], data_offsets[1], data_offsets[2]};
  result = MPI_Scatter( &set_counts[0], 3, MPI_LONG,
                        data_offsets,   3, MPI_LONG,
                        0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  setContentsOffset = data_offsets[0];
  setChildrenOffset = data_offsets[1];
  setParentsOffset = data_offsets[2];
  
    // Create set contents and set children tables
  if (myRank == 0)
  {
    rval = create_set_tables( all_counts[0], all_counts[1], all_counts[2] );
    if (MB_SUCCESS != rval) return rval;
  }
  
    // Send totals to all processors
  result = MPI_Bcast( all_counts, 3, MPI_LONG, 0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  writeSetContents = all_counts[0] > 0;
  writeSetChildren = all_counts[1] > 0;
  writeSetParents  = all_counts[2] > 0;

  START_SERIAL;  
  printdebug("Non-shared set contents: %ld local, %ld global, offset = %ld\n",
    data_counts[0], all_counts[0], data_offsets[0] );
  printdebug("Non-shared set children: %ld local, %ld global, offset = %ld\n",
    data_counts[1], all_counts[1], data_offsets[1] );
  printdebug("Non-shared set parents: %ld local, %ld global, offset = %ld\n",
    data_counts[2], all_counts[2], data_offsets[2] );
  END_SERIAL;
  
  return MB_SUCCESS;
}

void WriteHDF5Parallel::remove_remote_entities( MBEntityHandle relative,
                                                MBRange& range )
{
  MBRange result;
  result.merge( range.intersect( nodeSet.range ) );
  result.merge( range.intersect( setSet.range ) );  
  for (std::list<ExportSet>::iterator eiter = exportList.begin();
           eiter != exportList.end(); ++eiter )
  {
    result.merge( range.intersect( eiter->range ) );
  }
  //result.merge( range.intersect( myParallelSets ) );
  MBRange sets;
  int junk;
  sets.merge( MBRange::lower_bound( range.begin(), range.end(), CREATE_HANDLE(MBENTITYSET, 0, junk )), range.end() );
  remove_remote_sets( relative, sets );
  result.merge( sets );
  range.swap(result);
}

void WriteHDF5Parallel::remove_remote_sets( MBEntityHandle relative, 
                                            MBRange& range )
{
  MBRange result( range.intersect( setSet.range ) );
  //result.merge( range.intersect( myParallelSets ) );
  MBRange remaining( range.subtract( result ) );
  
  for(MBRange::iterator i = remaining.begin(); i != remaining.end(); ++i)
  {
      // Look for the first CPU which knows about both sets.
    int cpu;
    for (cpu = 0; cpu < numProc; ++cpu)
      if (cpuParallelSets[cpu].find(relative) != cpuParallelSets[cpu].end() &&
          cpuParallelSets[cpu].find(*i) != cpuParallelSets[cpu].end())
        break;
      // If we didn't find one, it may indicate a bug.  However,
      // it could also indicate that it is a link to some set that
      // exists on this processor but is not being written, because
      // the caller requested that some subset of the mesh be written.
    //assert(cpu < numProc);
      // If I'm the first set that knows about both, I'll handle it.
    if (cpu == myRank)
      result.insert( *i );
  }
  
  range.swap( result );
}
  
  

void WriteHDF5Parallel::remove_remote_entities( MBEntityHandle relative,
                                                std::vector<MBEntityHandle>& vect )
{
  MBRange intrsct;
  for (std::vector<MBEntityHandle>::const_iterator iter = vect.begin();
       iter != vect.end(); ++iter)
    intrsct.insert(*iter);
  remove_remote_entities( relative, intrsct );
  
  unsigned int read, write;
  for (read = write = 0; read < vect.size(); ++read)
  {
    if (intrsct.find(vect[read]) != intrsct.end())
    {
      if (read != write)
        vect[write] = vect[read];
      ++write;
    }
  }
  if (write != vect.size())
    vect.resize(write);
}

  

void WriteHDF5Parallel::remove_remote_sets( MBEntityHandle relative,
                                            std::vector<MBEntityHandle>& vect )
{
  MBRange intrsct;
  for (std::vector<MBEntityHandle>::const_iterator iter = vect.begin();
       iter != vect.end(); ++iter)
    intrsct.insert(*iter);
  remove_remote_sets( relative, intrsct );
  
  unsigned int read, write;
  for (read = write = 0; read < vect.size(); ++read)
  {
    if (intrsct.find(vect[read]) != intrsct.end())
    {
      if (read != write)
        vect[write] = vect[read];
      ++write;
    }
  }
  if (write != vect.size())
    vect.resize(write);
}

// Given a RemoteSetData object describing the set information for a 
// single tag (or tag pair), populate the list of parallel sets
// (this->parallelSets) with the per-entityset data.
MBErrorCode WriteHDF5Parallel::negotiate_remote_set_contents( RemoteSetData& data,
                                                              long* offsets /* long[3] */ )
{
  unsigned i;
  MBErrorCode rval;
  MBRange::const_iterator riter;
  int result;
  const unsigned count = data.range.size();
  const unsigned total = data.all_values.size();
  std::vector<int>::iterator viter, viter2;

    // Calculate counts for each meshset
  std::vector<long> local_sizes(3*count);
  std::vector<long>::iterator sizes_iter = local_sizes.begin();
  MBRange tmp_range;
  std::vector<MBEntityHandle> child_list;
  for (riter = data.range.begin(); riter != data.range.end(); ++riter)
  {
      // Count contents
    *sizes_iter = 0;
    tmp_range.clear();
    rval = iFace->get_entities_by_handle( *riter, tmp_range );
    remove_remote_entities( *riter, tmp_range );
    assert (MB_SUCCESS == rval);
    for (MBRange::iterator iter = tmp_range.begin(); iter != tmp_range.end(); ++iter)
      if (0 != idMap.find( *iter ))
        ++*sizes_iter;
    ++sizes_iter;
    
      // Count children
    *sizes_iter = 0;
    child_list.clear();
    rval = iFace->get_child_meshsets( *riter, child_list );
    remove_remote_sets( *riter, child_list );
    assert (MB_SUCCESS == rval);
    for (std::vector<MBEntityHandle>::iterator iter = child_list.begin();
         iter != child_list.end(); ++iter)
      if (0 != idMap.find( *iter ))
        ++*sizes_iter;
    
      // Count parents
    *sizes_iter = 0;
    child_list.clear();
    rval = iFace->get_parent_meshsets( *riter, child_list );
    remove_remote_sets( *riter, child_list );
    assert (MB_SUCCESS == rval);
    for (std::vector<MBEntityHandle>::iterator iter = child_list.begin();
         iter != child_list.end(); ++iter)
      if (0 != idMap.find( *iter ))
        ++*sizes_iter;
  }
  
    // Exchange sizes for sets between all processors.
  std::vector<long> all_sizes(3*total);
  std::vector<int> counts(numProc), displs(numProc);
  for (i = 0; i < (unsigned)numProc; i++)
    counts[i] = 3 * data.counts[i];
  displs[0] = 0;
  for (i = 1; i < (unsigned)numProc; i++)
    displs[i] = displs[i-1] + counts[i-1];
  result = MPI_Allgatherv( &local_sizes[0], 3*count, MPI_LONG,
                           &all_sizes[0], &counts[0], &displs[0], MPI_LONG,
                           MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);

  
    // Update information in-place in the array from the Allgatherv.
    
    // Change the corresponding sizes for the first instance of a tag
    // value such that it ends up being the total size of the set.
    // Change the size to -1 for the later instances of a tag value.
    //
    // For the sets that this processor has, update the offsets at
    // which the set data is to be written.  Store the offset of the data
    // on this processor for the set *relative* to the start of the
    // data of *the set*.
  std::vector<long> local_offsets(3*count);
  std::map<int,int> tagsort;  // Map of {tag value, index of first set w/ value}
  for (i = 0; i < total; ++i)
  {
    const std::map<int,int>::iterator p = tagsort.find( data.all_values[i] );
    const unsigned r = (unsigned)(i - data.displs[myRank]);  // offset in "local" array
    
      // If this is the first instance of this tag value, 
      // then the processor with this instance is responsible
      // for writing the tag description
    if ( p == tagsort.end() )  
    {
      tagsort[data.all_values[i]] = i;
        // If within the range for this processor, save offsets
      if (r < (unsigned)count) 
      {
        local_offsets[3*r] = local_offsets[3*r+1] = local_offsets[3*r+2] = 0;
      }
    }
      // Otherwise update the total size in the table
      // for the processor that is responsible for writing
      // the data and mark the data for the current processor
      // with a -1.
    else 
    {
        // If within the range for this processor, save offsets
      int j = p->second;
      if (r < (unsigned)count) 
      {
          // the offset for this processor, from the start of the data
          // for this group of sets, is the current total count for the
          // group of sets.
        local_offsets[3*r  ] = all_sizes[3*j  ];  // contents
        local_offsets[3*r+1] = all_sizes[3*j+1];  // children
        local_offsets[3*r+2] = all_sizes[3*j+2];  // parents
      }
      
        // update the total count for the set in the first position in
        // all_sizes at which the set occurs (the one corresponding to
        // the processor that owns the set.)
      all_sizes[3*j  ] += all_sizes[3*i  ]; // contents
      all_sizes[3*j+1] += all_sizes[3*i+1]; // children
      all_sizes[3*j+2] += all_sizes[3*i+2]; // parents
        // set the size to -1 in the positions corresponding to non-owning processor
      all_sizes[3*i  ] = all_sizes[3*i+1] = all_sizes[3*i+2] = -1;
    }
  }  
    
  
    // Store the total size of each set (rather than the
    // number of entities local to this processor) in the
    // local_sizes array for each meshset.  Only need this
    // for the sets this processor is writing the description
    // for, but it's easier to get it for all of them.
  sizes_iter = local_sizes.begin();
  viter = data.local_values.begin();
  for (riter = data.range.begin(); riter != data.range.end(); ++riter, ++viter)
  {
    const std::map<int,int>::iterator p = tagsort.find( *viter ); 
    assert( p != tagsort.end() );
    int j = 3 * p->second;
    *sizes_iter = all_sizes[j  ]; ++sizes_iter;  // contents
    *sizes_iter = all_sizes[j+1]; ++sizes_iter;  // children
    *sizes_iter = all_sizes[j+2]; ++sizes_iter;  // parents
  }
  
    // Now calculate the offset of the data for each (entire, parallel) set in
    // the set contents, children and parents tables.  offsets is long[3], and
    // is both input and output of this function.  We increment offsets by the
    // total count (over all processors) for each set such that it contains
    // the next open row in the table.  This will be passed back into this
    // function for the next tag (or tag pair) such that ultimately it will
    // contain the beginning of the non-shared set data in each of the three tables.
    // all_sizes is re-used to store the global offset in each table for each 
    // set with the tag.
  for (i = 0; i < all_sizes.size(); ++i)
  {
    if (all_sizes[i] >= 0) // value is -1 (from above) if not this processor
    {
      int j = i % 3;              // contents, children or parents list ?
      long tmp = offsets[j];      // save current, running offset
      offsets[j] += all_sizes[i]; // next set's offset is current plus the size of this set
      all_sizes[i] = tmp;         // size of this set is running offset.
    }
  }
  
    // Local offsets for this processor are stored as values relative to the
    // start of each set's data.  Convert them to offsets relative to the
    // start of all the set data.  Add the offset *from* the start of the set
    // data (local_offsets) to the offset *of* the start of the set data 
    // (stored in all_sizes in the previous loop) 
  std::vector<long>::iterator offset_iter = local_offsets.begin();
  viter = data.local_values.begin();
  for (riter = data.range.begin(); riter != data.range.end(); ++riter, ++viter)
  {
    const std::map<int,int>::iterator p = tagsort.find( *viter );
    assert( p != tagsort.end() );
    int j = 3 * p->second;
    *offset_iter += all_sizes[j  ]; ++offset_iter; // contents
    *offset_iter += all_sizes[j+1]; ++offset_iter; // children
    *offset_iter += all_sizes[j+2]; ++offset_iter; // parents
  }

#ifdef DEBUG  
START_SERIAL; if (counts[myRank]) {
std::string name1, name2;
iFace->tag_get_name( data.data_tag, name1 );
iFace->tag_get_name( data.filter_tag, name2 );
printdebug("Remote set data\n" );
printdebug("    %13s %13s owner local_offsets total_counts\n", name1.c_str(), name2.c_str());
for (unsigned d = 0; d < (unsigned)counts[myRank]; ++d) {
switch(d%3) {
  case 0: // data/contents
printdebug("   %13d %13d %5s %13d %12d\n", data.all_values[(d+displs[myRank])/3], 
 data.filter_value, 
 all_sizes[d+displs[myRank]] < 0 ? "no" : "yes", 
 local_offsets[d], local_sizes[d] );
  break;
  case 1: // children
printdebug("                          (children) %13d %12d\n", local_offsets[d], local_sizes[d] );
  break;
  case 2: // parents
printdebug("                           (parents) %13d %12d\n", local_offsets[d], local_sizes[d] );
  break;
} 
}
} 
END_SERIAL;
#endif
  
    // Store each parallel meshset in the list
  sizes_iter = local_sizes.begin();
  offset_iter = local_offsets.begin();
  std::vector<long>::iterator all_iter = all_sizes.begin() + displs[myRank];
  for (riter = data.range.begin(); riter != data.range.end(); ++riter)
  {
    ParallelSet info;
    info.handle = *riter;
    info.contentsOffset = *offset_iter; ++offset_iter;
    info.childrenOffset = *offset_iter; ++offset_iter;
    info.parentsOffset = *offset_iter; ++offset_iter;
    info.contentsCount = *sizes_iter; ++sizes_iter;
    info.childrenCount = *sizes_iter; ++sizes_iter;
    info.parentsCount = *sizes_iter; ++sizes_iter;
    info.description = *all_iter >= 0; all_iter += 3;
    parallelSets.push_back( info );
  }
  
  return MB_SUCCESS;
}

MBErrorCode WriteHDF5Parallel::fix_remote_set_ids( RemoteSetData& data, long first_id )
{
  const id_t id_diff = (id_t)(first_id - 1);
  MBRange::const_iterator i;
  std::vector<id_t> ids(data.range.size());
  std::vector<id_t>::iterator j = ids.begin();
  for (i = data.range.begin(); i != data.range.end(); ++i, ++j)
    *j = idMap.find( *i ) + id_diff;
  
  for (MBRange::const_pair_iterator pi = data.range.const_pair_begin();
       pi != data.range.const_pair_end(); ++pi) 
    idMap.erase( pi->first, pi->second - pi->first + 1);
  
  j = ids.begin();
  for (i = data.range.begin(); i != data.range.end(); ++i) 
    idMap.insert( *i, *j, 1 );
  
  return MB_SUCCESS;
}   


MBErrorCode WriteHDF5Parallel::write_shared_set_descriptions( hid_t table )
{
  const id_t start_id = setSet.first_id;
  MBErrorCode rval;
  mhdf_Status status;
  
  for( std::list<ParallelSet>::iterator iter = parallelSets.begin();
        iter != parallelSets.end(); ++iter)
  {
    if (!iter->description)
      continue;  // handled by a different processor
    
      // Get offset in table at which to write data
    id_t file_id = idMap.find( iter->handle );
    assert( file_id >= start_id );
    file_id -= start_id;
    
      // Get flag data
    unsigned int flags;
    rval = iFace->get_meshset_options( iter->handle, flags );
    assert( MB_SUCCESS == rval );
      
      // Write the data
    long data[4] = { iter->contentsOffset + iter->contentsCount - 1, 
                     iter->childrenOffset + iter->childrenCount - 1, 
                     iter->parentsOffset  + iter->parentsCount  - 1,
                     flags };
    mhdf_writeSetMeta( table, file_id, 1, H5T_NATIVE_LONG, data, &status );
    if (mhdf_isError(&status))
      printdebug("Meshset %d : %s\n", ID_FROM_HANDLE(iter->handle), mhdf_message(&status));
    assert( !mhdf_isError( &status ) );
  }

  return MB_SUCCESS;
}
    

MBErrorCode WriteHDF5Parallel::write_shared_set_contents( hid_t table )
{
  MBErrorCode rval;
  mhdf_Status status;
  std::vector<MBEntityHandle> handle_list;
  std::vector<id_t> id_list;
  
  for( std::list<ParallelSet>::iterator iter = parallelSets.begin();
        iter != parallelSets.end(); ++iter)
  {
    handle_list.clear();
    rval = iFace->get_entities_by_handle( iter->handle, handle_list );
    assert( MB_SUCCESS == rval );
    remove_remote_entities( iter->handle, handle_list );
    
    id_list.clear();
    vector_to_id_list( handle_list, id_list, true );
    if (id_list.empty())
      continue;
    
    mhdf_writeSetData( table, 
                       iter->contentsOffset, 
                       id_list.size(),
                       id_type,
                       &id_list[0],
                       &status );
    assert(!mhdf_isError(&status));
  }
  
  return MB_SUCCESS;
}
    

MBErrorCode WriteHDF5Parallel::write_shared_set_children( hid_t table )
{
  MBErrorCode rval;
  mhdf_Status status;
  std::vector<MBEntityHandle> handle_list;
  std::vector<id_t> id_list;
  
  printdebug("Writing %d parallel sets.\n", parallelSets.size());
  for( std::list<ParallelSet>::iterator iter = parallelSets.begin();
        iter != parallelSets.end(); ++iter)
  {
    handle_list.clear();
    rval = iFace->get_child_meshsets( iter->handle, handle_list );
    assert( MB_SUCCESS == rval );
    remove_remote_sets( iter->handle, handle_list );
    
    id_list.clear();
    vector_to_id_list( handle_list, id_list, true );
    if (!id_list.empty())
    {
      mhdf_writeSetParentsChildren( table, 
                                    iter->childrenOffset, 
                                    id_list.size(),
                                    id_type,
                                    &id_list[0],
                                    &status );
      assert(!mhdf_isError(&status));
    }
  }

  return MB_SUCCESS;
}
    

MBErrorCode WriteHDF5Parallel::write_shared_set_parents( hid_t table )
{
  MBErrorCode rval;
  mhdf_Status status;
  std::vector<MBEntityHandle> handle_list;
  std::vector<id_t> id_list;
  
  printdebug("Writing %d parallel sets.\n", parallelSets.size());
  for( std::list<ParallelSet>::iterator iter = parallelSets.begin();
        iter != parallelSets.end(); ++iter)
  {
    handle_list.clear();
    rval = iFace->get_parent_meshsets( iter->handle, handle_list );
    assert( MB_SUCCESS == rval );
    remove_remote_sets( iter->handle, handle_list );
    
    id_list.clear();
    vector_to_id_list( handle_list, id_list, true );
    if (!id_list.empty())
    {
      mhdf_writeSetParentsChildren( table, 
                                    iter->parentsOffset, 
                                    id_list.size(),
                                    id_type,
                                    &id_list[0],
                                    &status );
      assert(!mhdf_isError(&status));
    }
  }

  return MB_SUCCESS;
}

MBErrorCode WriteHDF5Parallel::write_finished()
{
  parallelSets.clear();
  cpuParallelSets.clear();
  //myParallelSets.clear();
  return WriteHDF5::write_finished();
}


class TagNameCompare {
  MBInterface* iFace;
  std::string name1, name2;
public:
  TagNameCompare( MBInterface* iface ) : iFace(iface) {}
  bool operator() (const WriteHDF5::SparseTag& t1, 
                   const WriteHDF5::SparseTag& t2);
};
bool TagNameCompare::operator() (const WriteHDF5::SparseTag& t1, 
                                 const WriteHDF5::SparseTag& t2)
{
  MBErrorCode rval;
  rval = iFace->tag_get_name( t1.tag_id, name1 );
  rval = iFace->tag_get_name( t2.tag_id, name2 );
  return name1 < name2;
}  

void WriteHDF5Parallel::sort_tags_by_name( )
{
  tagList.sort( TagNameCompare( iFace ) );
}


MBErrorCode WriteHDF5Parallel::communicate_remote_ids( MBEntityType type )
{
  int result;
  MBErrorCode rval;

    // Get the export set for the specified type
  ExportSet* export_set = 0;
  if (type == MBVERTEX)
    export_set = &nodeSet;
  else if(type == MBENTITYSET)
    export_set = &setSet;
  else
  {
    for (std::list<ExportSet>::iterator esiter = exportList.begin();
         esiter != exportList.end(); ++esiter)
      if (esiter->type == type)
      {
        export_set = &*esiter;
        break;
      }
  }
  assert(export_set != NULL);
  
    // Get the ranges in the set
  std::vector<unsigned long> myranges;
  MBRange::const_pair_iterator p_iter = export_set->range.const_pair_begin();
  const MBRange::const_pair_iterator p_end = export_set->range.const_pair_end();
  for ( ; p_iter != p_end; ++p_iter)
  {
    myranges.push_back( (*p_iter).first );
    myranges.push_back( (*p_iter).second );
  }

  START_SERIAL;
  printdebug("%s ranges to communicate:\n", MBCN::EntityTypeName(type));
  for (unsigned long xx = 0; xx != myranges.size(); xx+=2)
    printdebug("  %lu - %lu\n", myranges[xx], myranges[xx+1] );
  END_SERIAL;
  
    // Communicate the number of ranges and the start_id for
    // each processor.
  std::vector<int> counts(numProc), offsets(numProc), displs(numProc);
  int mycount = myranges.size();
  int mystart = export_set->first_id + export_set->offset;
  result = MPI_Allgather( &mycount, 1, MPI_INT, &counts[0], 1, MPI_INT, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  result = MPI_Allgather( &mystart, 1, MPI_INT, &offsets[0], 1, MPI_INT, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
    // Communicate the ranges 
  displs[0] = 0;
  for (int i = 1; i < numProc; ++i)
    displs[i] = displs[i-1] + counts[i-1];
  std::vector<unsigned long> allranges( displs[numProc-1] + counts[numProc-1] );
  result = MPI_Allgatherv( &myranges[0], myranges.size(), MPI_UNSIGNED_LONG,
                           &allranges[0], &counts[0], &displs[0],
                           MPI_UNSIGNED_LONG, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
  MBTag global_id_tag;
  rval = iFace->tag_get_handle( PARALLEL_GID_TAG_NAME, global_id_tag );
  assert(MB_SUCCESS == rval);
  
    // Set file IDs for each communicated entity
    
    // For each processor
  for (int proc = 0; proc < numProc; ++proc)
  {
    if (proc == myRank)
      continue;
    
      // Get data for corresponding processor
    const int offset = offsets[proc];
    const int count = counts[proc];
    const unsigned long* const ranges = &allranges[displs[proc]];
    
      // For each geometry meshset in the interface
    MBRange::iterator r_iter = MBRange::lower_bound( remoteMesh[proc].begin(),
                                                     remoteMesh[proc].end(),
                                                     CREATE_HANDLE(type,0,result) );
    MBRange::iterator r_stop = MBRange::lower_bound( r_iter,
                                                     remoteMesh[proc].end(),
                                                     CREATE_HANDLE(type+1,0,result) );
    for ( ; r_iter != r_stop; ++r_iter)
    {
      MBEntityHandle entity = *r_iter;

        // Get handle on other processor
      MBEntityHandle global;
      rval = iFace->tag_get_data( global_id_tag, &entity, 1, &global );
      assert(MB_SUCCESS == rval);

        // Find corresponding fileid on other processor.
        // This could potentially be n**2, but we will assume that
        // the range list from each processor is short (typically 1).
      int j, steps = 0;
      unsigned long low, high;
      for (j = 0; j < count; j += 2)
      {
        low = ranges[j];
        high = ranges[j+1];
        if (low <= global && high >= global)
          break;
        steps += (high - low) + 1;
      }
      if (j >= count) {
      printdebug("*** handle = %u, type = %d, id = %d, proc = %d\n",
      (unsigned)global, (int)(iFace->type_from_handle(global)), (int)(iFace->id_from_handle(global)), proc);
      for (int ii = 0; ii < count; ii+=2) 
      printdebug("***  %u to %u\n", (unsigned)ranges[ii], (unsigned)ranges[ii+1] );
      MBRange junk; junk.insert( global );
      print_type_sets( iFace, myRank, numProc, junk );
      }
      assert(j < count);
      int fileid = offset + steps + (global - low);
      RangeMap<MBEntityHandle,id_t>::iterator ri = idMap.insert( entity, fileid, 1 );
      assert( ri != idMap.end() );
    } // for(r_iter->range)
  } // for(each processor)
  
  return MB_SUCCESS;
}
