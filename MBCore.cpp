/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

#ifdef WIN32
// turn off warnings that say they debugging identifier has been truncated
// this warning comes up when using some STL containers
#pragma warning(disable : 4786)
#endif

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include "MBVersion.h"
#include "MBCore.hpp"
#include "TagServer.hpp"
#include "MeshSetSequence.hpp"
#include "assert.h"
#include "AEntityFactory.hpp"
#include "MBReadUtil.hpp"
#include "MBWriteUtil.hpp"
#include "MBCN.hpp"
#include "HigherOrderFactory.hpp"
#include "EntitySequenceManager.hpp"
#include "MBError.hpp"
#include "MBReaderWriterSet.hpp"
#include "MBReaderIface.hpp"
#include "MBWriterIface.hpp"
#include "ReadParallel.hpp"

#ifdef HDF5_FILE
#  include "WriteHDF5.hpp"
   typedef WriteHDF5 DefaultWriter;
#elif defined(NETCDF_FILE)
#  include "WriteNCDF.hpp"
   typedef WriteNCDF DefaultWriter;
#else
#  include "WriteVtk.hpp"
   typedef WriteVtk DefaultWriter;
#endif
#ifdef USE_MPI
#include "mpi.h"
#endif
#include "MBTagConventions.hpp"
#include "ExoIIUtil.hpp"
#include "EntitySequence.hpp"
#include "FileOptions.hpp"
#include "MBParallelComm.hpp"
#ifdef LINUX
# include <dlfcn.h>
# include <dirent.h>
#endif


#ifdef XPCOM_MB
#include "nsMemory.h"
#endif

using namespace std;

const char *MBCore::errorStrings[] = {
  "MB_SUCCESS",
  "MB_INDEX_OUT_OF_RANGE",
  "MB_TYPE_OUT_OF_RANGE",
  "MB_MEMORY_ALLOCATION_FAILED",
  "MB_ENTITY_NOT_FOUND",
  "MB_MULTIPLE_ENTITIES_FOUND",
  "MB_TAG_NOT_FOUND",
  "MB_FILE_DOES_NOT_EXIST",
  "MB_FILE_WRITE_ERROR",
  "MB_NOT_IMPLEMENTED",
  "MB_ALREADY_ALLOCATED",
  "MB_FAILURE",
};

static inline MBMeshSet* get_mesh_set( const EntitySequenceManager* sm,
                                       MBEntityHandle h )
{
  MBEntitySequence* seq;
  if (MBENTITYSET != TYPE_FROM_HANDLE(h) ||
      MB_SUCCESS != sm->find( h, seq ) || 
      !seq->is_valid_entity(h))
    return 0;
  return reinterpret_cast<MeshSetSequence*>(seq)->get_set(h);
}

//! Constructor
MBCore::MBCore( int rank, int num_procs ) 
  : procInfo( rank, num_procs )
{
#ifdef XPCOM_MB
  NS_INIT_ISUPPORTS();
#endif

  if (initialize() != MB_SUCCESS)
  {
    printf("Error initializing MB\n");
    exit(1);
  }
}

//! destructor
MBCore::~MBCore()
{
  if(mMBWriteUtil)
    delete mMBWriteUtil;
  if(mMBReadUtil) 
    delete mMBReadUtil;
  
  mMBWriteUtil = NULL;
  mMBReadUtil = NULL;

  deinitialize();
}


MBErrorCode MBCore::initialize()
{
  geometricDimension = 3;
  materialTag      = 0;
  neumannBCTag     = 0;
  dirichletBCTag   = 0;
  geomDimensionTag = 0;
  globalIdTag      = 0;

  tagServer = new TagServer();
  if (!tagServer)
    return MB_MEMORY_ALLOCATION_FAILED;
  
  sequenceManager = new EntitySequenceManager( procInfo );
  if (!sequenceManager)
    return MB_MEMORY_ALLOCATION_FAILED;

  aEntityFactory = new AEntityFactory(this);
  if (!aEntityFactory)
    return MB_MEMORY_ALLOCATION_FAILED;

  mError = new MBError;

  mMBWriteUtil = NULL;
  mMBReadUtil = NULL;
    
    // Readers and writers try to get pointers to above utils.
    // Do this after pointers are initialized. (Pointers should
    // really be initialized in constructor to avoid this kind
    // of thing -- j.kraftcheck.)
  readerWriterSet = new MBReaderWriterSet( this, mError );
  if (!readerWriterSet)
    return MB_MEMORY_ALLOCATION_FAILED;
  
  material_tag();
  neumannBC_tag();
  dirichletBC_tag();
  geom_dimension_tag();
  globalId_tag();

  return MB_SUCCESS;
}

MBEntityHandle MBCore::get_root_set() 
{
  return 0;
}

void MBCore::deinitialize()
{
  if (aEntityFactory)
    delete aEntityFactory;

  aEntityFactory = 0;

  if (tagServer)
    delete tagServer;

  tagServer = 0;
  
  if (sequenceManager)
    delete sequenceManager;

  sequenceManager = 0;
  
  delete readerWriterSet;
  readerWriterSet = 0;

  if(mError)
    delete mError;
  mError = 0;
}

MBErrorCode MBCore::query_interface(const std::string& iface_name, void** iface)
{
  if(iface_name == "MBReadUtilIface")
  {
    if(mMBReadUtil)
      *iface = (MBReadUtilIface*)mMBReadUtil;
    else
      *iface = (MBReadUtilIface*)(mMBReadUtil = new MBReadUtil(this, mError));
    return MB_SUCCESS;
  }
  else if(iface_name == "MBWriteUtilIface")
  {
    if(mMBWriteUtil)
      *iface = (MBWriteUtilIface*)mMBWriteUtil;
    else
      *iface = (MBWriteUtilIface*)(mMBWriteUtil = new MBWriteUtil(this, mError));
    return MB_SUCCESS;
  }
  else if(iface_name == "MBReaderWriterSet")
  {
    *iface = reader_writer_set();
    return MB_SUCCESS;
  }
  else if(iface_name == "ExoIIInterface")
  {
    *iface = (void*)(ExoIIInterface*) new ExoIIUtil(this);
    return MB_SUCCESS;
  }
  return MB_FAILURE;
}


MBErrorCode MBCore::release_interface(const std::string& iface_name, void* iface)
{
  if(iface == NULL)
    return MB_FAILURE;

  if(iface_name == "MBReadUtilIface")
  {
      // Is it possible to crash here?  We should fail gracefully instead.
    return MB_SUCCESS;
  }
  else if(iface_name == "MBWriteUtilIface")
  {
    return MB_SUCCESS;
  }
  else if(iface_name == "MBReaderWriterSet")
  {
    return MB_SUCCESS;
  }
  else if(iface_name == "ExoIIInterface")
  {
    delete (ExoIIInterface*)iface;
    return MB_SUCCESS;
  }
  
  return MB_FAILURE;
}


#ifdef XPCOM_MB
// provides basic implementation of nsISupports methods
NS_IMPL_ISUPPORTS1_CI(MBCore, MBInterface);
#endif

int MBCore::QueryInterface(const MBuuid& uuid, MBUnknownInterface** iface)
{
  *iface = 0;
  if(uuid == IDD_MBUnknown)
    *iface = this;
  if(uuid == IDD_MBCore)
    *iface = this;
  else
    return 0;
  return 1;
}

float MBCore::impl_version( std::string *version_string )
{
  if (version_string)
    *version_string = MB_VERSION_STRING;
  
  return MB_VERSION_MAJOR + MB_VERSION_MINOR / 100.0f;
}

//! get the type from a handle, returns type
MBEntityType MBCore::type_from_handle(const MBEntityHandle handle) const
{
  return TYPE_FROM_HANDLE(handle);
}
  
//! get the id from a handle, returns id
MBEntityID MBCore::id_from_handle(const MBEntityHandle handle) const
{
  return ID_FROM_HANDLE(handle);
}

//! get a handle from an id and type
MBErrorCode MBCore::handle_from_id( const MBEntityType type, 
                                    const MBEntityID id, 
                                    MBEntityHandle& handle) const
{
  int err;
  handle = CREATE_HANDLE(type, id, err);

    //check to see if handle exists 
  MBEntitySequence *dummy_seq = 0;
  MBErrorCode error_code = sequence_manager()->find(handle, dummy_seq);
  if(error_code == MB_SUCCESS && type != MBVERTEX)
    error_code = dummy_seq->is_valid_entity(handle) ? MB_SUCCESS : MB_ENTITY_NOT_FOUND;

  return error_code; 
}

int MBCore::dimension_from_handle(const MBEntityHandle handle) const
{
  return MBCN::Dimension(TYPE_FROM_HANDLE(handle));
}

//! load mesh from data in file
//! NOTE: if there is mesh already present, the new mesh will be added
MBErrorCode  MBCore::load_mesh( const char *file_name,
                                const int* block_id_list,
                                const int num_blocks )
{
  MBEntityHandle file_set;
  return load_file( file_name, file_set, 0, MATERIAL_SET_TAG_NAME, block_id_list, num_blocks );
}

MBErrorCode MBCore::load_file( const char* file_name,
                               MBEntityHandle& file_set,
                               const char* options,
                               const char* set_tag_name,
                               const int* set_tag_values,
                               int num_set_tag_values )
{
  if (num_set_tag_values < 0)
    return MB_INDEX_OUT_OF_RANGE;
   
  FileOptions opts(options);
  file_set = 0;
  
  MBErrorCode rval;
  const MBReaderWriterSet* set = reader_writer_set();
  
    // convert from to old block-based reader interface
  const int* block_id_list = 0;
  int num_blocks = 0;
  if (set_tag_name) {
    if (strcmp(set_tag_name, MATERIAL_SET_TAG_NAME ))
      return MB_NOT_IMPLEMENTED;
    block_id_list = set_tag_values;
    num_blocks = num_set_tag_values;
  }
  
    // if reading in parallel, call a different reader
  std::string parallel_opt;
  rval = opts.get_option( "PARALLEL", parallel_opt);
  if (MB_SUCCESS == rval && !parallel_opt.empty()) {
    return ReadParallel(this).load_file(file_name, file_set, opts,
                                        block_id_list, num_blocks);
  }

    // otherwise try using the file extension to select a reader
  MBReaderIface* reader = set->get_file_extension_reader( file_name );
  if (reader)
  { 
    rval = reader->load_file( file_name, file_set, opts, block_id_list, num_blocks );
    delete reader;
  }
  else
  {  
      // Try all the readers
    MBReaderWriterSet::iterator iter;
    for (iter = set->begin(); iter != set->end(); ++iter)
    {
      MBReaderIface* reader = iter->make_reader( this );
      if (NULL != reader)
      {
        rval = reader->load_file( file_name, file_set, opts, block_id_list, num_blocks );
        delete reader;
        if (MB_SUCCESS == rval)
          break;
      }
    }
  }
  
  return rval; 
}

MBErrorCode  MBCore::write_mesh(const char *file_name,
                                  const MBEntityHandle *output_list,
                                  const int num_sets)
{
  return write_file( file_name, 0, 0, output_list, num_sets );
}

MBErrorCode MBCore::write_file( const char* file_name,
                                const char* file_type,
                                const char* options,
                                const MBEntityHandle* output_sets,
                                int num_output_sets,
                                const MBTag* tag_list,
                                int num_tags )
{
  MBRange range;
  std::copy( output_sets, output_sets+num_output_sets, mb_range_inserter(range) );
  return write_file( file_name, file_type, options, range, tag_list, num_tags );
}

MBErrorCode MBCore::write_file( const char* file_name,
                                const char* file_type,
                                const char* options_string,
                                const MBRange& output_sets,
                                const MBTag* tag_list,
                                int num_tags )
{
    // convert range to vector
  std::vector<MBEntityHandle> list( output_sets.size() );
  std::copy( output_sets.begin(), output_sets.end(), list.begin() );
  
    // parse some options
  FileOptions opts( options_string );
  MBErrorCode rval;
  
  rval = opts.get_null_option( "CREATE" );
  if (rval == MB_TYPE_OUT_OF_RANGE) {
    mError->set_last_error( "Unexpected value for CREATE option\n" );
    return MB_FAILURE;
  }
  bool overwrite = (rval == MB_ENTITY_NOT_FOUND);

    // Get the file writer
  MBReaderWriterSet::iterator i;
  if (file_type) {
    i = reader_writer_set()->handler_by_name( file_type );
    if (i == reader_writer_set()->end()) {
      mError->set_last_error( "Unknown file type: %s\n", file_type );
      return MB_NOT_IMPLEMENTED;
    }
  }
  else {
    std::string ext = MBReaderWriterSet::extension_from_filename( file_name );
    i = reader_writer_set()->handler_from_extension( ext );
  }
  
  MBWriterIface* writer;
  if (i == reader_writer_set()->end())
    writer = new DefaultWriter(this);
  else
    writer = i->make_writer( this );
  
  if (!writer) {
    mError->set_last_error( "File format supported for reading only.\n" );
    return MB_NOT_IMPLEMENTED;
  }
  
    // write the file
  std::vector<std::string> qa_records;
  rval = writer->write_file(file_name, overwrite, opts, &list[0], list.size(), qa_records );
  delete writer;
  
  return rval;
}
   
  

//! deletes all mesh entities from this datastore
MBErrorCode MBCore::delete_mesh()
{

  MBErrorCode result = MB_SUCCESS;

    // perform all deinitialization procedures to clean up
  if (aEntityFactory)
    delete aEntityFactory;
  aEntityFactory = new AEntityFactory(this);

  tagServer->reset_all_data();
  
  if (sequenceManager)
    delete sequenceManager;
  sequenceManager = new EntitySequenceManager( procInfo );

  return result;
}

  //! get overall geometric dimension
MBErrorCode MBCore::get_dimension(int &dim) const
{
  dim = geometricDimension;
  return MB_SUCCESS;
}

  //! set overall geometric dimension
  /** Returns error if setting to 3 dimensions, mesh has been created, and 
   *  there are only 2 dimensions on that mesh
   */
MBErrorCode MBCore::set_dimension(const int dim) 
{
    // check to see if current dimension is smaller
  if (geometricDimension < dim) 
  {
      // need to check the number of entities 
    int num;
    /*MBErrorCode result = */ get_number_entities_by_dimension(0, geometricDimension, num);
    
      // test written to be more readable but possibly less efficient
      //if (MB_SUCCESS != result) return MB_FAILURE;
      //else if (0 != num && dim == 2 && ycoordTag == 0) return MB_FAILURE;
      //else if (0 != num && dim == 3 && (ycoordTag == 0 || zcoordTag == 0)) return MB_FAILURE;
      //TODO -- replace this with not using xcoordTag, etc...
  }
    
    // if we got here, it's ok to set dimension
  geometricDimension = dim;
  return MB_SUCCESS;
}

  //! get blocked vertex coordinates for all vertices
  /** Blocked = all x, then all y, etc. 
   */
MBErrorCode MBCore::get_vertex_coordinates(std::vector<double> &coords) const
{
    // INEFFICIENT implementation for now, until we get blocked tag access
  MBRange vertices;
  MBErrorCode result = get_entities_by_type(0, MBVERTEX, vertices);
  if (MB_SUCCESS != result) return result;
  
    // the least we can do is resize the vector and only go through the 
    // vertex list once
  int num_verts = vertices.size();
  int vec_pos = 0;
  double xyz[3];
  coords.resize(geometricDimension*num_verts);
  for (MBRange::iterator it = vertices.begin(); it != vertices.end(); it++) 
  {
    result = get_coords(&(*it), 1, xyz);
    if (MB_SUCCESS != result) return result;

    coords[vec_pos] = xyz[0];
    coords[num_verts+vec_pos] = xyz[1];
    coords[2*num_verts+vec_pos] = xyz[2];

    vec_pos++;
  }
  
  return result;
}

MBErrorCode  MBCore::get_coords(const MBRange& entities, double *coords) const
{

  int node_index = 0;

  // lets get ready to iterate the entity sequences and copy data
  // we'll iterate the map in the sequence manager and
  // the entities range at the same time

  MBRange::const_iterator range_iter = entities.begin();
  MBRange::const_iterator range_iter_end = entities.end();

  std::map<MBEntityHandle, MBEntitySequence*>::const_iterator seq_iter
    = sequence_manager()->entity_map(MBVERTEX)->begin();
  
  std::map<MBEntityHandle, MBEntitySequence*>::const_iterator seq_iter_end
    = sequence_manager()->entity_map(MBVERTEX)->end();

  // lets find the entity sequence which holds the first entity
  std::map<MBEntityHandle, MBEntitySequence*>::const_iterator seq_iter_lookahead = seq_iter;
  seq_iter_lookahead++;
  for( ; seq_iter_lookahead != seq_iter_end && 
      seq_iter_lookahead->second->get_start_handle() < *range_iter; )
  {
    ++seq_iter;
    ++seq_iter_lookahead;
  }

  // a look ahead iterator
  MBRange::const_iterator range_iter_lookahead = range_iter;

  // our main loop
  for(; range_iter != range_iter_end && seq_iter != seq_iter_end; /* ++ is handled in loop*/ )
  {
    // find a range that fits in the current entity sequence
    for(; range_iter_lookahead != range_iter_end && 
        *range_iter_lookahead <= seq_iter->second->get_end_handle(); 
        ++range_iter_lookahead)
    {}

    double* coord_array[3];
    static_cast<VertexEntitySequence*>(seq_iter->second)->get_coordinate_arrays(
        coord_array[0], coord_array[1], coord_array[2]);
    MBEntityHandle start_ent = seq_iter->second->get_start_handle();

    // for each of the entities in this entity sequence, copy data
    for(MBRange::const_iterator tmp_iter = range_iter; 
        tmp_iter != range_iter_lookahead;
        ++tmp_iter)
    {
      coords[node_index] = coord_array[0][*tmp_iter - start_ent];
      node_index++;
      coords[node_index] = coord_array[1][*tmp_iter - start_ent];
      node_index++;
      coords[node_index] = coord_array[2][*tmp_iter - start_ent];
      node_index++;
    }

    // go to the next entity sequence
    ++seq_iter;
    // start with the next entities
    range_iter = range_iter_lookahead;
  }
  return MB_SUCCESS;
}

/**\author Jason Kraftcheck <kraftche@cae.wisc.edu> - 2007-5-15 */
MBErrorCode MBCore::get_coords( const MBRange& entities, 
                                double *x_coords,
                                double *y_coords,
                                double *z_coords ) const
{
    // declare some iterators
  MBRange::const_iterator iter = entities.begin();
  const MBRange::const_iterator end = entities.end();
  std::map<MBEntityHandle, MBEntitySequence*>::const_iterator seq_iter
    = sequence_manager()->entity_map(MBVERTEX)->begin();
  const std::map<MBEntityHandle, MBEntitySequence*>::const_iterator seq_end
    = sequence_manager()->entity_map(MBVERTEX)->end();

    // if no vertices in database, give up now
  if (seq_iter == seq_end)
    return MB_ENTITY_NOT_FOUND;

    // loop over all entities in range
  while(iter != end) {
      // find sequence for vertex at *iter
    while (seq_iter->second->get_end_handle() < *iter)
       if (++seq_iter == seq_end)
        return MB_ENTITY_NOT_FOUND;
    if (*iter < seq_iter->second->get_start_handle())
      return MB_ENTITY_NOT_FOUND;

      // get data for sequence
    VertexEntitySequence* seq = static_cast<VertexEntitySequence*>(seq_iter->second);
    const MBEntityHandle start_ent = seq->get_start_handle();
    const MBEntityHandle end_ent = seq->get_end_handle();
    double *seq_x, *seq_y, *seq_z;
    seq->get_coordinate_arrays( seq_x, seq_y, seq_z );

      // for each entity in range and sequence 
    while (iter != end && *iter <= end_ent) {
    
        // get block of consecutive handles
      const MBEntityHandle range_start = *iter;
      iter = iter.end_of_block();
      if (*iter > end_ent) // trim block to those in the sequence
        iter -= *iter - end_ent;
       
      const MBEntityID offset = range_start - start_ent;
      const MBEntityID count = *iter - range_start + 1;
      ++iter; // done with iter for this iteration, advance to next handle
      
        // copy coordinate data
      if (x_coords) {
        memcpy( x_coords, seq_x + offset, sizeof(double)*count );
        x_coords += count;
      }
      if (y_coords) {
        memcpy( y_coords, seq_y + offset, sizeof(double)*count );
        y_coords += count;
      }
      if (z_coords) {
        memcpy( z_coords, seq_z + offset, sizeof(double)*count );
        z_coords += count;
      }
    }
  }

  return MB_SUCCESS;
}

MBErrorCode  MBCore::get_coords(const MBEntityHandle* entities, 
                                  const int num_entities, 
                                  double *coords) const
{
  MBErrorCode status;
  MBEntitySequence* seq;
  const MBEntityHandle* const end = entities + num_entities;

  for(const MBEntityHandle* iter = entities; iter != end; ++iter)
  {
    if(TYPE_FROM_HANDLE(*iter) != MBVERTEX)
      return MB_TYPE_OUT_OF_RANGE;

    status = sequence_manager()->find(*iter, seq);
    if(status != MB_SUCCESS || !seq->is_valid_entity(*iter) )
      return MB_ENTITY_NOT_FOUND;
    
    static_cast<VertexEntitySequence*>(seq)->get_coordinates(*iter, coords);
    coords += 3;
  }

  return MB_SUCCESS; 
}


MBErrorCode  MBCore::get_coords(const MBEntityHandle entity_handle, 
                                  const double *& x, const double *& y, const double *& z) const
{

  MBErrorCode status = MB_TYPE_OUT_OF_RANGE;

  if ( TYPE_FROM_HANDLE(entity_handle) == MBVERTEX )
  {
    MBEntitySequence* seq = 0;
    status = sequence_manager()->find(entity_handle, seq);

    if (seq == 0 || status != MB_SUCCESS || !seq->is_valid_entity(entity_handle) )
      return MB_ENTITY_NOT_FOUND;

    status = static_cast<VertexEntitySequence*>(seq)->get_coordinates_ref(entity_handle, 
                                                                          x, y, z);

  }

  return status; 

}

//! set the coordinate information for this handle if it is of type Vertex
//! otherwise, return an error
MBErrorCode  MBCore::set_coords(MBEntityHandle *entity_handles, const int num_entities, 
                                  const double *coords)
{

  MBErrorCode status = MB_SUCCESS;

  int i, j = 0;

  for (i = 0; i < num_entities; i++) {
    if ( TYPE_FROM_HANDLE(entity_handles[i]) == MBVERTEX )
    {
      MBEntitySequence* seq = 0;
      status = sequence_manager()->find(entity_handles[i], seq);

      if (seq != 0 && status == MB_SUCCESS) {
        status = static_cast<VertexEntitySequence*>(seq)->set_coordinates(entity_handles[i], coords[j], coords[j+1], coords[j+2]);
        j += 3;
      }
    }
    else if (status == MB_SUCCESS)
      status = MB_TYPE_OUT_OF_RANGE;
  }

  return status; 

}

//! set the coordinate information for this handle if it is of type Vertex
//! otherwise, return an error
MBErrorCode  MBCore::set_coords(MBRange entity_handles, const double *coords)
{

  MBErrorCode status = MB_SUCCESS;

  int j = 0;

  for (MBRange::iterator rit = entity_handles.begin(); rit != entity_handles.end(); rit++) {
    if ( TYPE_FROM_HANDLE(*rit) == MBVERTEX )
    {
      MBEntitySequence* seq = 0;
      status = sequence_manager()->find(*rit, seq);

      if (seq != 0 && status == MB_SUCCESS) {
        status = static_cast<VertexEntitySequence*>(seq)->set_coordinates(*rit, coords[j], coords[j+1], coords[j+2]);
        j += 3;
      }
    }
    else if (status == MB_SUCCESS)
      status = MB_TYPE_OUT_OF_RANGE;
  }

  return status; 

}

  //! get global connectivity array for specified entity type
  /**  Assumes just vertices, no higher order nodes
   */
MBErrorCode MBCore::get_connectivity_by_type(const MBEntityType type, 
                                               std::vector<MBEntityHandle> &connect) const
{
    // inefficient implementation until we get blocked tag access
  
    // get the range of entities of this type
  MBRange this_range;
  MBErrorCode result = get_entities_by_type(0, type, this_range);
  
  int num_ents = this_range.size();
  connect.reserve(num_ents*MBCN::VerticesPerEntity(type));
  
    // now loop over these entities, getting connectivity for each
  for (MBRange::iterator this_it = this_range.begin(); 
       this_it != this_range.end();
       this_it++)
  {
    const MBEntityHandle *connect_vec;
    result = get_connectivity(*this_it, connect_vec, num_ents, true);
    if (MB_SUCCESS != result) 
      return result;
    connect.insert(connect.end(), &connect_vec[0], &connect_vec[num_ents]); 
  }
  
  return MB_SUCCESS;
}
  

//! get the connectivity for element /handles.  For non-element handles, return an error
MBErrorCode  MBCore::get_connectivity(const MBEntityHandle *entity_handles, 
                                      const int num_handles,
                                      MBRange &connectivity,
                                      bool topological_connectivity) const
{
  std::vector<MBEntityHandle> tmp_connect;
  MBErrorCode result = get_connectivity(entity_handles, num_handles, tmp_connect,
                                        topological_connectivity);
  if (MB_SUCCESS != result) return result;
  
  std::copy(tmp_connect.begin(), tmp_connect.end(), mb_range_inserter(connectivity));
  return result;
}

//! get the connectivity for element /handles.  For non-element handles, return an error
MBErrorCode  MBCore::get_connectivity(const MBEntityHandle *entity_handles, 
                                      const int num_handles,
                                      std::vector<MBEntityHandle> &connectivity,
                                      bool topological_connectivity) const
{
  if (num_handles == 0) return MB_FAILURE;
  
    // Make sure the entity should have a connectivity.
  MBEntityType type = TYPE_FROM_HANDLE(entity_handles[0]);

  int i;
  connectivity.clear();

  // handle the special case where someone asks for the connectivity
  // of a vertex.  This actually kind of makes sense for a sphere element.
  if (type == MBVERTEX)
  {
    connectivity.reserve(num_handles);
    std::copy(entity_handles, entity_handles+num_handles, 
              std::back_inserter(connectivity));
    return MB_SUCCESS;
  }

    // WARNING: This is very dependent on the ordering of the MBEntityType enum
  if(type <= MBVERTEX || type >= MBENTITYSET)
    return MB_TYPE_OUT_OF_RANGE;

  MBErrorCode result = MB_SUCCESS, temp_result;
  for (i = 0; i < num_handles; i++) {
    MBEntitySequence* seq = 0;

      // We know that connectivity is stored in an EntitySequence so jump straight
      // to the entity sequence
    temp_result = sequence_manager()->find(entity_handles[i], seq);
    if (seq == NULL || !seq->is_valid_entity(entity_handles[i])) {
      result = MB_ENTITY_NOT_FOUND;
      continue;
    }
    else if (temp_result != MB_SUCCESS) {
      result = temp_result;
      continue;
    }

    ElementEntitySequence *elem_seq = static_cast<ElementEntitySequence*>(seq);
      // let's be smart about this...
    temp_result = elem_seq->get_connectivity(entity_handles[i], connectivity,
                                             topological_connectivity);
    if (MB_SUCCESS != temp_result) {
      result = temp_result;
      continue;
    }
  }
  
  return result;
}

//! get the connectivity for element handles.  For non-element handles, return an error
MBErrorCode MBCore::get_connectivity(const MBEntityHandle entity_handle, 
                                     const MBEntityHandle*& connectivity,
                                     int& number_nodes,
                                     bool topological_connectivity) const
{

  MBErrorCode status;

    // Make sure the entity should have a connectivity.
  MBEntityType type = TYPE_FROM_HANDLE(entity_handle);
  
    // WARNING: This is very dependent on the ordering of the MBEntityType enum
  if(type < MBVERTEX || type >= MBENTITYSET)
    return MB_TYPE_OUT_OF_RANGE;

  else if (type == MBVERTEX) {
    return MB_FAILURE;
  }
  
  MBEntitySequence* seq = 0;

    // We know that connectivity is stored in an EntitySequence so jump straight
    // to the entity sequence
  status = sequence_manager()->find(entity_handle, seq);
  if (seq == 0 || status != MB_SUCCESS || !seq->is_valid_entity(entity_handle)) 
    return MB_ENTITY_NOT_FOUND;

  return static_cast<ElementEntitySequence*>(seq)->get_connectivity(entity_handle, connectivity,
                                                                    number_nodes,
                                                                    topological_connectivity);
}

//! set the connectivity for element handles.  For non-element handles, return an error
MBErrorCode  MBCore::set_connectivity(const MBEntityHandle entity_handle, 
                                      MBEntityHandle *connect,
                                      const int num_connect)
{
  MBErrorCode status = MB_FAILURE;

    // Make sure the entity should have a connectivity.
    // WARNING: This is very dependent on the ordering of the MBEntityType enum
  MBEntityType type = TYPE_FROM_HANDLE(entity_handle);
  
  MBEntitySequence* seq = 0;

  if (type < MBVERTEX || type > MBENTITYSET)
    return MB_TYPE_OUT_OF_RANGE;
  
  status = sequence_manager()->find(entity_handle, seq);
  if (seq == 0 || status != MB_SUCCESS || !seq->is_valid_entity(entity_handle))
    return (status != MB_SUCCESS ? status : MB_ENTITY_NOT_FOUND);

  const MBEntityHandle* old_conn;
  int len;
  status = static_cast<ElementEntitySequence*>(seq)->get_connectivity(entity_handle, old_conn, len);
  if (status != MB_SUCCESS) return status;

  aEntityFactory->notify_change_connectivity(
    entity_handle, old_conn, connect, num_connect);
#ifdef MOAB_WITH_REFCOUNT
  decrement_reference_count( old_conn, num_connect );
#endif
  
  status = static_cast<ElementEntitySequence*>(seq)->set_connectivity(entity_handle, 
                                                                      connect, num_connect);
  if (status != MB_SUCCESS) {
    aEntityFactory->notify_change_connectivity(
      entity_handle, connect, old_conn, num_connect);
#ifdef MOAB_WITH_REFCOUNT
    increment_reference_count( old_conn, num_connect );
  }
  else {
    increment_reference_count( connect, num_connect );
#endif
  }
  
  return status;
}

MBErrorCode MBCore::get_adjacencies(const MBEntityHandle *from_entities,
                                      const int num_entities,
                                      const int to_dimension,
                                      const bool create_if_missing,
                                      std::vector<MBEntityHandle> &adj_entities,
                                      const int operation_type)
{
  MBErrorCode result;
  if (num_entities == 1) {
    if(to_dimension == 0 && TYPE_FROM_HANDLE(from_entities[0]) != MBPOLYHEDRON)
      result = get_connectivity(&from_entities[0], 1, adj_entities);
    else
      result = aEntityFactory->get_adjacencies(from_entities[0], to_dimension, 
                                               create_if_missing, adj_entities);
    std::vector<MBEntityHandle>::iterator iter = 
      std::remove(adj_entities.begin(), adj_entities.end(), MBEntityHandle(0));
    if(iter != adj_entities.end())
      adj_entities.erase(iter, adj_entities.end());
    return result;
  }
  
  MBRange temp_range, temp_range2;
  std::copy(from_entities, from_entities+num_entities,
            mb_range_inserter(temp_range));
  
  result = get_adjacencies(temp_range, to_dimension, create_if_missing, temp_range2,
                           operation_type);
  if (MB_SUCCESS != result) 
    return result;
  
  adj_entities.clear();
  adj_entities.reserve(temp_range2.size());
  for (MBRange::const_iterator it = temp_range2.begin();
       it != temp_range2.end();
       ++it)
    adj_entities.push_back(*it);
  
  return MB_SUCCESS;
}

MBErrorCode MBCore::get_adjacencies(const MBEntityHandle *from_entities,
                                    const int num_entities,
                                    const int to_dimension,
                                      const bool create_if_missing,
                                      MBRange &adj_entities,
                                      const int operation_type)
{
  MBRange tmp_from_entities;
  std::copy(from_entities, from_entities+num_entities, mb_range_inserter(tmp_from_entities));
  return get_adjacencies(tmp_from_entities, to_dimension, create_if_missing, 
                         adj_entities, operation_type);
}
  
MBErrorCode MBCore::get_adjacencies(const MBRange &from_entities,
                                      const int to_dimension,
                                      const bool create_if_missing,
                                      MBRange &adj_entities,
                                      const int operation_type)
{
  if (operation_type != MBInterface::INTERSECT &&
      operation_type != MBInterface::UNION) return MB_FAILURE;

  if(from_entities.size() == 0)
    return MB_SUCCESS;

  MBRange temp_range;
  std::vector<MBEntityHandle> temp_vec;
  MBErrorCode result = MB_SUCCESS, tmp_result;

  for (MBRange::const_iterator from_it = from_entities.begin(); 
       from_it != from_entities.end(); from_it++) 
  {
      // running results kept in adj_entities; clear temp_vec and temp_range, which are working space
    temp_vec.clear();
    temp_range.clear();

      // get the next set of adjacencies
    if(to_dimension == 0 && TYPE_FROM_HANDLE(*from_it) != MBPOLYHEDRON)
      tmp_result = get_connectivity(&(*from_it), 1, temp_vec);
    else
      tmp_result = aEntityFactory->get_adjacencies(*from_it, to_dimension, 
                                                   create_if_missing, temp_vec);
    
    if (MB_SUCCESS != tmp_result) result = tmp_result;
    
      // if we're on the first iteration and we didn't come in with entities,
      // just get the first results and move on
    if (adj_entities.empty() && from_it == from_entities.begin()) {
      std::copy(temp_vec.begin(), temp_vec.end(), mb_range_inserter(adj_entities));
      continue;
    }

      // operate on the vectors
    if (operation_type == MBInterface::INTERSECT) {
        // only have to sort if we're doing intersection
      std::sort(temp_vec.begin(), temp_vec.end());
      std::set_intersection(adj_entities.begin(), adj_entities.end(), 
                            temp_vec.begin(), temp_vec.end(),
                            mb_range_inserter(temp_range));
      adj_entities = temp_range;
    }
    else if (operation_type == MBInterface::UNION) {
      std::copy(temp_vec.begin(), temp_vec.end(), mb_range_inserter(adj_entities));
    }
  }

  return result;
}

MBErrorCode MBCore::add_adjacencies(const MBEntityHandle entity_handle, 
                                    const MBEntityHandle *adjacencies,
                                    const int num_handles,
                                    bool both_ways)
{
  MBErrorCode result = MB_SUCCESS, temp_result;
  
  for (const MBEntityHandle *it = adjacencies; 
       it != adjacencies+num_handles; it++) {
    temp_result = aEntityFactory->add_adjacency(entity_handle, *it, both_ways);
    if (MB_SUCCESS != temp_result) result = temp_result;
  }

  return result;
}

MBErrorCode MBCore::add_adjacencies(const MBEntityHandle entity_handle, 
                                    MBRange &adjacencies,
                                    bool both_ways)
{
  MBErrorCode result = MB_SUCCESS, temp_result;
  
  for (MBRange::iterator rit = adjacencies.begin(); rit != adjacencies.end(); rit++) {
    temp_result = aEntityFactory->add_adjacency(entity_handle, *rit, both_ways);
    if (MB_SUCCESS != temp_result) result = temp_result;
  }

  return result;
}

MBErrorCode MBCore::remove_adjacencies(const MBEntityHandle entity_handle,
                                       const MBEntityHandle *adjacencies,
                                       const int num_handles)
{
  MBErrorCode result = MB_SUCCESS, temp_result;
  
  for (const MBEntityHandle *it = adjacencies; 
       it != adjacencies+num_handles; it++) {
    temp_result = aEntityFactory->remove_adjacency(entity_handle, *it);
    if (MB_SUCCESS != temp_result) result = temp_result;
    temp_result = aEntityFactory->remove_adjacency(*it, entity_handle);
    if (MB_SUCCESS != temp_result) result = temp_result;
  }

  return result;
}

MBErrorCode MBCore::get_entities_by_dimension(const MBEntityHandle meshset,
                                                const int dimension, 
                                                MBRange &entities,
                                                const bool recursive) const
{
  MBErrorCode result = MB_SUCCESS;
  if (meshset) {
    MBEntitySequence* seq;
    result = sequence_manager()->find( meshset, seq );
    if (MB_SUCCESS != result)
      return result;
    MeshSetSequence* mseq = reinterpret_cast<MeshSetSequence*>(seq);
    result = mseq->get_dimension( meshset, dimension, entities, recursive );
  }
  else if (dimension > 3) {
    result = sequence_manager()->get_entities( MBENTITYSET, entities );
  } 
  else {
    for (MBEntityType this_type = MBCN::TypeDimensionMap[dimension].first;
         this_type <= MBCN::TypeDimensionMap[dimension].second;
         this_type++) {
      result = sequence_manager()->get_entities( this_type, entities );
      if (MB_SUCCESS != result)
        break;
    }
  }

  return result;
}

MBErrorCode MBCore::get_entities_by_type(const MBEntityHandle meshset,
                                           const MBEntityType type, 
                                           MBRange &entities,
                                           const bool recursive) const
{
  if (recursive && type == MBENTITYSET)  // will never return anything
    return MB_TYPE_OUT_OF_RANGE;
  
  MBErrorCode result = MB_SUCCESS;
  if (meshset) {
    MBEntitySequence* seq;
    result = sequence_manager()->find( meshset, seq );
    if (MB_SUCCESS != result)
      return result;
    MeshSetSequence* mseq = reinterpret_cast<MeshSetSequence*>(seq);
    result = mseq->get_type( meshset, type, entities, recursive );
  }  
  else
    result = sequence_manager()->get_entities( type, entities );

  return result;
}

MBErrorCode MBCore::get_entities_by_type_and_tag(const MBEntityHandle meshset,
                                                   const MBEntityType type,
                                                   const MBTag *tags,
                                                   const void* const* values,
                                                   const int num_tags,
                                                   MBRange &entities,
                                                   const int condition,
                                                   const bool recursive) const
{
  if (recursive && type == MBENTITYSET)  // will never return anything
    return MB_TYPE_OUT_OF_RANGE;

  MBErrorCode result;
  MBRange tmp_range;

  result = get_entities_by_type( meshset, type, tmp_range, recursive );
  if (MB_SUCCESS != result)
    return result;

    // if range is empty, return right away; if intersecting condition, 
    // empty the list too
  if (tmp_range.empty()) {
    if (MBInterface::INTERSECT == condition) entities.clear();
    return MB_SUCCESS;
  }
  else if (!entities.empty() && MBInterface::INTERSECT == condition) {
    entities = entities.intersect(tmp_range);
    if (entities.empty()) return MB_SUCCESS;
    tmp_range = entities;
  }
    
  result = tagServer->get_entities_with_tag_values(tmp_range, type, 
                                                   tags, values, num_tags, 
                                                   entities, condition); 
  
  return result;
}

MBErrorCode MBCore::get_entities_by_handle(const MBEntityHandle meshset,
                                             MBRange &entities,
                                             const bool recursive) const
{
  MBErrorCode result = MB_SUCCESS;
  if (meshset) {
    MBEntitySequence* seq;
    result = sequence_manager()->find( meshset, seq );
    if (MB_SUCCESS != result)
      return result;
    MeshSetSequence* mseq = reinterpret_cast<MeshSetSequence*>(seq);
    result = mseq->get_entities( meshset, entities, recursive );
  }  
  else {
    // iterate backards so range insertion is quicker
    for (MBEntityType type = MBENTITYSET; type >= MBVERTEX; --type) {
      result = sequence_manager()->get_entities( type, entities );
      if (MB_SUCCESS != result) 
        break;
    }
  }

  return result;
}


MBErrorCode MBCore::get_entities_by_handle(const MBEntityHandle meshset,
                                   std::vector<MBEntityHandle> &entities,
                                   const bool recursive) const
{
  MBErrorCode result;
  if (recursive || !meshset) {
    MBRange tmp_range;
    result = get_entities_by_handle( meshset, tmp_range, recursive);
    size_t offset = entities.size();
    entities.resize( offset + tmp_range.size() );
    std::copy( tmp_range.begin(), tmp_range.end(), entities.begin() + offset );
  }
  else {
    MBEntitySequence* seq;
    result = sequence_manager()->find( meshset, seq );
    if (MB_SUCCESS != result)
      return result;
    MeshSetSequence* mseq = reinterpret_cast<MeshSetSequence*>(seq);
    result = mseq->get_entities( meshset, entities );
  }  
  return result;
}

  //! get # entities of a given dimension
MBErrorCode MBCore::get_number_entities_by_dimension(const MBEntityHandle meshset,
                                                       const int dim, 
                                                       int &number,
                                                       const bool recursive) const
{
  MBErrorCode result = MB_SUCCESS;
 
  if (!meshset) {
    number = 0;
    for (MBEntityType this_type = MBCN::TypeDimensionMap[dim].first;
         this_type <= MBCN::TypeDimensionMap[dim].second;
         this_type++) {
      MBEntityID dummy = 0;
      result = sequence_manager()->get_number_entities( this_type, dummy );
      if (result != MB_SUCCESS) 
        break;
      number += dummy;
    }
  }
  else {
    MBEntitySequence* seq;
    result = sequence_manager()->find( meshset, seq );
    if (MB_SUCCESS != result)
      return result;
    MeshSetSequence* mseq = reinterpret_cast<MeshSetSequence*>(seq);
    result = mseq->num_dimension( meshset, dim, number, recursive );
  }  
  
  return result;
}

//! returns the number of entities with a given type and tag
MBErrorCode MBCore::get_number_entities_by_type(const MBEntityHandle meshset,
                                                  const MBEntityType type, 
                                                  int& num_ent,
                                                  const bool recursive) const
{
  MBErrorCode result = MB_SUCCESS;

  if (recursive && type == MBENTITYSET)  // will never return anything
    return MB_TYPE_OUT_OF_RANGE;
  
  if (meshset) {
    MBEntitySequence* seq;
    result = sequence_manager()->find( meshset, seq );
    if (MB_SUCCESS != result)
      return result;
    MeshSetSequence* mseq = reinterpret_cast<MeshSetSequence*>(seq);
    result = mseq->num_type( meshset, type, num_ent, recursive );
  }
  else {
    MBEntityID tmp;
    result = sequence_manager()->get_number_entities( type, tmp );
    num_ent = tmp;
  }
  
  return result;
}

MBErrorCode MBCore::get_number_entities_by_type_and_tag(const MBEntityHandle meshset,
                                                          const MBEntityType type,
                                                          const MBTag *tag_handles,
                                                          const void** values,
                                                          const int num_tags,
                                                          int &num_entities,
                                                          const bool recursive) const
{
  MBRange dum_ents;
  MBErrorCode result = get_entities_by_type_and_tag(meshset, type, tag_handles, values, num_tags, 
                                                     dum_ents, recursive);
  num_entities = dum_ents.size();
  return result;
}

MBErrorCode MBCore::get_number_entities_by_handle(const MBEntityHandle meshset,
                                          int& num_ent,
                                          const bool recursive) const
{
  MBErrorCode result;
  if (meshset) {
    MBEntitySequence* seq;
    result = sequence_manager()->find( meshset, seq );
    if (MB_SUCCESS != result)
      return result;
    MeshSetSequence* mseq = reinterpret_cast<MeshSetSequence*>(seq);
    return mseq->num_entities( meshset, num_ent, recursive );
  }

  num_ent = 0;
  for (MBEntityType this_type = MBVERTEX;
       this_type < MBMAXTYPE;
       this_type++) {
    int dummy = 0;
    result = get_number_entities_by_type(0, this_type, dummy);
    if (result != MB_SUCCESS) {
      num_ent = 0;
      return result;
    }
    num_ent += dummy;
  }

  return MB_SUCCESS;
}

//! return the tag data for a given EntityHandle and MBTag
MBErrorCode  MBCore::tag_get_data(const MBTag tag_handle, 
                                    const MBEntityHandle* entity_handles, 
                                    const int num_entities,
                                    void *tag_data) const
{
  if (NULL == entity_handles && 0 == num_entities)
    return tagServer->get_mesh_data(tag_handle, tag_data);

  else return tagServer->get_data(tag_handle, entity_handles, num_entities, tag_data);
}

//! return the tag data for a given EntityHandle and MBTag
MBErrorCode  MBCore::tag_get_data(const MBTag tag_handle, 
                                    const MBRange& entity_handles,
                                    void *tag_data) const
{
  return tagServer->get_data(tag_handle, entity_handles, tag_data);
}

//! set the data  for given EntityHandles and MBTag
MBErrorCode  MBCore::tag_set_data(const MBTag tag_handle, 
                                    const MBEntityHandle* entity_handles, 
                                    const int num_entities,
                                    const void *tag_data)
{
  MBErrorCode rval;
#ifdef MOAB_WITH_REFCOUNT
  std::vector<MBEntityHandle> old_data;
  std::vector<MBEntityHandle> unique_handles;
  const TagInfo* tag_info = tagServer->get_tag_info( tag_handle );
  if(!tag_info)
    return MB_TAG_NOT_FOUND;
  int count;
#endif

  if (NULL == entity_handles && 0 == num_entities) {
#ifdef MOAB_WITH_REFCOUNT   
    if (tag_info->get_data_type() == MB_TYPE_HANDLE) {
      count = tag_info->get_size() / sizeof(MBEntityHandle);
      old_data.resize( count );
      rval = tagServer->get_mesh_data( tag_handle, &old_data[0] );
      if (MB_SUCCESS == rval) {
        rval = decrement_reference_count( &old_data[0], count );
        if (MB_SUCCESS != rval)
          return rval;
      }
      else 
        old_data.clear();
    }
#endif  
    rval = tagServer->set_mesh_data(tag_handle, tag_data);
  }
  else {
  //verify handles
    MBEntitySequence* seq;
    const MBEntityHandle* iter;
    const MBEntityHandle* end = entity_handles + num_entities;
    for(iter = entity_handles; iter != end; ++iter)
    {
      if (TYPE_FROM_HANDLE(*iter) == MBENTITYSET) continue;

      else if(sequenceManager->find(*iter, seq) != MB_SUCCESS)
        return MB_ENTITY_NOT_FOUND;
    }
    
#ifdef MOAB_WITH_REFCOUNT   
    if (tag_info->get_data_type() == MB_TYPE_HANDLE) {
      unique_handles.resize( num_entities );
      std::copy( entity_handles, entity_handles+num_entities, unique_handles.begin() );
      std::sort( unique_handles.begin(), unique_handles.end() );
      unique_handles.erase( std::unique( unique_handles.begin(), unique_handles.end() ), unique_handles.end() );
      
      count = tag_info->get_size() * unique_handles.size() / sizeof(MBEntityHandle);
      old_data.resize( count );
      rval = tagServer->get_mesh_data( tag_handle, &old_data[0] );
      if (MB_SUCCESS == rval) {
        rval = decrement_reference_count( &old_data[0], count );
        if (MB_SUCCESS != rval)
          return rval;
      }
      else 
        old_data.clear();
    }
#endif

    rval = tagServer->set_data(tag_handle, entity_handles, num_entities, tag_data);
  }

#ifdef MOAB_WITH_REFCOUNT
  if (tag_info->get_data_type() == MB_TYPE_HANDLE) {
    if (MB_SUCCESS == rval) {
      if (unique_handles.size()) {
        tag_server()->get_data( tag_handle, &unique_handles[0], unique_handles.size(), &old_data[0] );
        rval = increment_reference_count( &old_data[0], count );
      }
      else {
        rval = increment_reference_count( (const MBEntityHandle*)tag_data, count );
      }
    }
    else if (!old_data.empty())
      increment_reference_count( &old_data[0], count );  
  }
#endif

  return rval;
}

//! set the data  for given EntityHandles and MBTag
MBErrorCode  MBCore::tag_set_data(const MBTag tag_handle, 
                                    const MBRange& entity_handles, 
                                    const void *tag_data)
{
  //verify handles
  MBRange::const_iterator iter;
  MBErrorCode result;
  for(iter = entity_handles.begin(); 
      iter != entity_handles.end() && TYPE_FROM_HANDLE(*iter) != MBENTITYSET;
      ++iter)
  {
    MBEntitySequence* seq = NULL;
    result = sequenceManager->find(*iter, seq);
    if(result != MB_SUCCESS)
      return result;
  }

#ifdef MOAB_WITH_REFCOUNT
  std::vector<MBEntityHandle> old_data;
  const TagInfo* tag_info = tagServer->get_tag_info( tag_handle );
  if(!tag_info)
    return MB_TAG_NOT_FOUND;
  int count;
  if (tag_info->get_data_type() == MB_TYPE_HANDLE) {
    count = tag_info->get_size() * entity_handles.size() / sizeof(MBEntityHandle);
    old_data.resize( count );
    result = tagServer->get_data( tag_handle, entity_handles, &old_data[0] );
    if (MB_SUCCESS != result)
      old_data.clear();
    else {
      result = decrement_reference_count( &old_data[0], count );
      if (MB_SUCCESS != result)
        return result;
    }
  }
#endif

  result = tagServer->set_data(tag_handle, entity_handles, tag_data);

#ifdef MOAB_WITH_REFCOUNT
  if (tag_info->get_data_type() == MB_TYPE_HANDLE) {
    if (MB_SUCCESS == result) 
      result = increment_reference_count( (const MBEntityHandle*)tag_data, count );
    else if (!old_data.empty())
      increment_reference_count( &old_data[0], count );  
  }
#endif

  return result;

}

//! adds a sparse tag for this specific EntityHandle/tag_name combination
MBErrorCode MBCore::tag_create(const char *tag_name,
                                 const int tag_size, 
                                 const MBTagType tag_type,
                                 MBTag &tag_handle, 
                                 const void *default_value)
{
  MBDataType data_type = (tag_type == MB_TAG_BIT) ? MB_TYPE_BIT : MB_TYPE_OPAQUE;
  return tag_create( tag_name, tag_size, tag_type, data_type, tag_handle, default_value, false );
}

MBErrorCode MBCore::tag_create( const char* name,
                                const int size,
                                const MBTagType storage,
                                const MBDataType data,
                                MBTag& handle,
                                const void* def_val,
                                bool use_existing )
{
#ifdef MOAB_WITH_REFCOUNT
    // don't allow a non-zero default value for handle-type tags
    // if reference counting is enabled.
  if (data == MB_TYPE_HANDLE && def_val) {
    const MBEntityHandle* handles = (const MBEntityHandle*)def_val;
    for (size_t i = 0; i < size / sizeof(MBEntityHandle); ++i)
      if (handles[i])
        return MB_FAILURE;
  }
#endif

  MBErrorCode rval = tagServer->add_tag( name, size, storage, data, handle, def_val );

    // If it is okay to use an existing tag of the same name, check that it 
    // matches the input values.  NOTE: we don't check the storage type for 
    // the tag because the choice of dense vs. sparse is a matter of optimi-
    // zation, not correctness.  Bit tags require a redundant MB_TYPE_BIT 
    // for the data type, so we catch those when we check the data type.
  if (rval == MB_ALREADY_ALLOCATED && use_existing) {
    handle = tagServer->get_handle( name );
    const TagInfo* info = tagServer->get_tag_info( handle );
    if (info->get_size() == size &&  info->get_data_type() == data)  {
        // If we were not passed a default value, the caller is presumably
        // OK with an arbitrary default value, so its OK if there is one
        // set for the tag.
      if (!def_val)
        rval = MB_SUCCESS;
        // If caller specified a default value, there MUST be an existing one
        // that matches.  We could just set a default value for the tag, but
        // given the dense tag representation, it isn't feasible to change
        // the default value once the tag has been created.  For now, don't 
        // bother because there isn't any mechanism to set it for sparse tags 
        // anyway.
      else if (info->default_value() && !memcmp(info->default_value(), def_val, size))
        rval = MB_SUCCESS;
    }
  }

  return rval;
}

//! removes the tag from the entity
MBErrorCode  MBCore::tag_delete_data(const MBTag tag_handle, 
                                       const MBEntityHandle *entity_handles,
                                       const int num_handles)
{
  if (PROP_FROM_TAG_HANDLE(tag_handle) == MB_TAG_DENSE)
    return MB_FAILURE;

  MBErrorCode status = MB_SUCCESS, temp_status;
#ifdef MOAB_WITH_REFCOUNT
  std::vector<MBEntityHandle> old_data;
  const TagInfo* tag_info = tagServer->get_tag_info( tag_handle );
  if(!tag_info)
    return MB_TAG_NOT_FOUND;
  if (tag_info->get_data_type() == MB_TYPE_HANDLE)
    old_data.resize( tag_info->get_size()/sizeof(MBEntityHandle) );
#endif

  for (int i = 0; i < num_handles; i++) {
    if (0 == entity_handles[i]) {
#ifdef MOAB_WITH_REFCOUNT
      if (tag_info->get_data_type() == MB_TYPE_HANDLE &&
          MB_SUCCESS == tagServer->get_mesh_data( tag_handle, &old_data[0] ))
        decrement_reference_count( &old_data[0], old_data.size() );
#endif
      temp_status = tagServer->remove_mesh_data(tag_handle);
    }
    else {
#ifdef MOAB_WITH_REFCOUNT
      if (tag_info->get_data_type() == MB_TYPE_HANDLE &&
          MB_SUCCESS == tagServer->get_data( tag_handle, entity_handles+i, 1, &old_data[0] ))
        decrement_reference_count( &old_data[0], old_data.size() );
#endif
      temp_status = tagServer->remove_data(tag_handle, entity_handles[i]);
    }
    if (temp_status != MB_SUCCESS) status = temp_status;
  }

  return status;
}

//! removes the tag from the entity
MBErrorCode  MBCore::tag_delete_data(const MBTag tag_handle, 
                                     const MBRange &entity_handles)
{
  if (PROP_FROM_TAG_HANDLE(tag_handle) == MB_TAG_DENSE)
    return MB_FAILURE;

#ifdef MOAB_WITH_REFCOUNT
  std::vector<MBEntityHandle> old_data;
  const TagInfo* tag_info = tagServer->get_tag_info( tag_handle );
  if(!tag_info)
    return MB_TAG_NOT_FOUND;
  if (tag_info->get_data_type() == MB_TYPE_HANDLE)
    old_data.resize( tag_info->get_size()/sizeof(MBEntityHandle) );
#endif

  MBErrorCode status = MB_SUCCESS, temp_status;
  for (MBRange::const_iterator it = entity_handles.begin(); it != entity_handles.end(); it++) {
#ifdef MOAB_WITH_REFCOUNT
    if (tag_info->get_data_type() == MB_TYPE_HANDLE &&
        MB_SUCCESS == tagServer->get_data( tag_handle, &*it, 1, &old_data[0] ))
      decrement_reference_count( &old_data[0], old_data.size() );
#endif
    temp_status = tagServer->remove_data(tag_handle, *it);
    if (temp_status != MB_SUCCESS) status = temp_status;
  }

  return status;
}

//! removes the tag from MB
MBErrorCode  MBCore::tag_delete(MBTag tag_handle)
{
#ifdef MOAB_WITH_REFCOUNT
  const TagInfo* tag_info = tagServer->get_tag_info( tag_handle );
  if(!tag_info)
    return MB_TAG_NOT_FOUND;
  if (tag_info->get_data_type() == MB_TYPE_HANDLE) {
    std::vector<MBEntityHandle> tag_data;
    for (MBEntityType t = MBVERTEX; t <= MBENTITYSET; ++t) {
      MBRange entities;
      tagServer->get_entities( tag_handle, t, entities );
      tag_data.resize( entities.size() );
      tagServer->get_data( tag_handle, entities, &tag_data[0] );
      decrement_reference_count( &tag_data[0], tag_data.size() );
    }
  }
#endif
  return tag_server()->remove_tag(tag_handle);
}

//! gets the tag name string for the tag_handle
MBErrorCode  MBCore::tag_get_name(const MBTag tag_handle, 
                                    std::string& tag_name) const
{
  const TagInfo* tag_info = tagServer->get_tag_info( tag_handle );
  if(!tag_info)
    return MB_TAG_NOT_FOUND;
  
  tag_name = tag_info->get_name();
  return MB_SUCCESS;

}

//! gets tag handle from its name.
//! the type must be specified because the same name is valid for multiple types
MBErrorCode  MBCore::tag_get_handle(const char *tag_name, 
                                      MBTag &tag_handle) const
{
  MBErrorCode status = MB_TAG_NOT_FOUND;

  tag_handle = tagServer->get_handle( tag_name );
  
  if (tag_handle != 0)
  {
    status = MB_SUCCESS;
  }

  return status;
}

  //! get size of tag in bytes
MBErrorCode MBCore::tag_get_size(const MBTag tag_handle, int &tag_size) const
{
  const TagInfo* tag_info = tagServer->get_tag_info( tag_handle );
  if(!tag_info)
    return MB_TAG_NOT_FOUND;
  
  tag_size = tag_info->get_size();
  return MB_SUCCESS;
}

MBErrorCode MBCore::tag_get_data_type( const MBTag handle, 
                                       MBDataType& type ) const
{
  const TagInfo* info = tagServer->get_tag_info( handle );
  if (!info)
    return MB_TAG_NOT_FOUND;
  
  type = info->get_data_type();
  return MB_SUCCESS;
}

  //! get default value of the tag
MBErrorCode MBCore::tag_get_default_value(const MBTag tag_handle, void *def_value) const
{
  const TagInfo* tag_info = tagServer->get_tag_info( tag_handle );
  if(!tag_info)
    return MB_TAG_NOT_FOUND;

  if (NULL == def_value) return MB_FAILURE;
  
  if (tag_info->default_value() == NULL)
    return MB_ENTITY_NOT_FOUND;
  
  memcpy(def_value, tag_info->default_value(), tag_info->get_size());

  return MB_SUCCESS;
}

  //! get type of tag (sparse, dense, etc.; 0 = dense, 1 = sparse, 2 = bit, 3 = static)
MBErrorCode MBCore::tag_get_type(const MBTag tag_handle, MBTagType &tag_type) const
{
  tag_type = PROP_FROM_TAG_HANDLE(tag_handle);
  return tag_type < MB_TAG_LAST+1 ? MB_SUCCESS : MB_TAG_NOT_FOUND;
}

  //! get handles for all tags defined
MBErrorCode MBCore::tag_get_tags(std::vector<MBTag> &tag_handles) const
{
  return tagServer->get_tags(tag_handles);
}

  //! Get handles for all tags defined on this entity
MBErrorCode MBCore::tag_get_tags_on_entity(const MBEntityHandle entity,
                                            std::vector<MBTag> &tag_handles) const 
{
  if (0 == entity)
    return tagServer->get_mesh_tags(tag_handles);
  else return tagServer->get_tags(entity, tag_handles);
}

MBTag MBCore::material_tag()
{
  if (0 == materialTag)
    tagServer->add_tag(MATERIAL_SET_TAG_NAME, sizeof(int), 
                       MB_TAG_SPARSE, MB_TYPE_INTEGER, materialTag);
  return materialTag;
}

MBTag MBCore::neumannBC_tag()
{
  if (0 == neumannBCTag)
    tagServer->add_tag(NEUMANN_SET_TAG_NAME, sizeof(int), 
                       MB_TAG_SPARSE, MB_TYPE_INTEGER, neumannBCTag);
  return neumannBCTag;
}

MBTag MBCore::dirichletBC_tag()
{
  if (0 == dirichletBCTag)
    tagServer->add_tag(DIRICHLET_SET_TAG_NAME, sizeof(int), 
                       MB_TAG_SPARSE, MB_TYPE_INTEGER, dirichletBCTag);
  return dirichletBCTag;
}

MBTag MBCore::globalId_tag()
{
  if (0 == globalIdTag)
    tagServer->add_tag(GLOBAL_ID_TAG_NAME, sizeof(int), 
                       MB_TAG_DENSE, MB_TYPE_INTEGER, globalIdTag);
  return globalIdTag;
}

MBTag MBCore::geom_dimension_tag()
{
  if (0 == geomDimensionTag)
    tagServer->add_tag(GEOM_DIMENSION_TAG_NAME, sizeof(int), 
                       MB_TAG_SPARSE, MB_TYPE_INTEGER, geomDimensionTag);
  return geomDimensionTag;
}

//! creates an element based on the type and connectivity.  returns a handle and error code
MBErrorCode MBCore::create_element(const MBEntityType type, 
                                   const MBEntityHandle *connectivity,
                                   const int num_nodes, 
                                   MBEntityHandle &handle)
{
  return create_element( type, procInfo.rank(), connectivity, num_nodes, handle );
}

MBErrorCode MBCore::create_element( const MBEntityType type,
                                    const unsigned processor_id,
                                    const MBEntityHandle* connectivity,
                                    const int num_nodes,
                                    MBEntityHandle& handle )
{
    // make sure we have enough vertices for this entity type
  if(num_nodes < MBCN::VerticesPerEntity(type))
    return MB_FAILURE;
  
  if (processor_id >= procInfo.size())
    return MB_INDEX_OUT_OF_RANGE;
  
  MBErrorCode status = sequence_manager()->create_element(type, processor_id, connectivity, num_nodes, handle);
  if (MB_SUCCESS == status)
    status = aEntityFactory->notify_create_entity( handle, connectivity, num_nodes); 

  return status;
}  

//! creates a vertex based on coordinates, returns a handle and error code
MBErrorCode MBCore::create_vertex(const double coords[3], MBEntityHandle &handle )
{
  return create_vertex( procInfo.rank(), coords, handle );
}

MBErrorCode MBCore::create_vertex( const unsigned processor_id, const double* coords, MBEntityHandle& handle )
{
  if (processor_id >= procInfo.size())
    return MB_INDEX_OUT_OF_RANGE;
    
    // get an available vertex handle
  return sequence_manager()->create_vertex( processor_id, coords, handle );
}

//! merges two  entities
MBErrorCode MBCore::merge_entities( MBEntityHandle entity_to_keep, 
                                      MBEntityHandle entity_to_remove,
                                      bool auto_merge,
                                      bool delete_removed_entity)
{
  if (auto_merge) return MB_FAILURE;
  
    // The two entities to merge must be of the same type
  MBEntityType type_to_keep = TYPE_FROM_HANDLE(entity_to_keep);

  if (type_to_keep != TYPE_FROM_HANDLE(entity_to_remove))
    return MB_TYPE_OUT_OF_RANGE;

    // Make sure both entities exist before trying to merge.
  MBEntitySequence* seq = 0;
  MBErrorCode result, status;
  status = sequence_manager()->find(entity_to_keep, seq);
  if(seq == 0 || status != MB_SUCCESS ||
     !seq->is_valid_entity(entity_to_keep))
    return MB_ENTITY_NOT_FOUND;
  status = sequence_manager()->find(entity_to_remove, seq);
  if(seq == 0 || status != MB_SUCCESS ||
     !seq->is_valid_entity(entity_to_remove))
    return MB_ENTITY_NOT_FOUND;
  
    // If auto_merge is not set, all sub-entities should
    // be merged if the entities are to be merged.
  int ent_dim = MBCN::Dimension(type_to_keep);
  if(ent_dim > 0)
  {
    std::vector<MBEntityHandle> conn, conn2;

    result = get_connectivity(&entity_to_keep, 1, conn);
    if(result != MB_SUCCESS)
      return result;
    result = get_connectivity(&entity_to_remove, 1, conn2);
    if(result != MB_SUCCESS)
      return result;

      // Check to see if we can merge before pulling adjacencies.
    int dum1, dum2;
    if(!auto_merge && 
       (conn.size() != conn2.size() ||
        !MBCN::ConnectivityMatch(&conn[0], &conn2[0], conn.size(), dum1, dum2)))
      return MB_FAILURE;
  }

  result = aEntityFactory->merge_adjust_adjacencies(entity_to_keep, entity_to_remove);
  
  if (delete_removed_entity) 
    result = delete_entities(&entity_to_remove, 1);

  return result;
}


MBErrorCode MBCore::delete_entity( MBEntityHandle handle )
{
  MBErrorCode result;

    // tell AEntityFactory that this element is going away
  result = aEntityFactory->notify_delete_entity(handle);
  if (MB_SUCCESS != result) 
    return result;

    // reset and/or clean out data associated with this entity handle
  result = tagServer->reset_data(handle);
  if (MB_SUCCESS != result) 
    return result;

  if (TYPE_FROM_HANDLE(handle) == MBENTITYSET) {
    if (MBMeshSet* ptr = get_mesh_set( sequence_manager(), handle )) {
      int j, count;
      const MBEntityHandle* rel;
      ptr->clear( handle, a_entity_factory() );
      rel = ptr->get_parents( count );
      for (j = 0; j < count; ++j)
        remove_child_meshset( rel[j], handle );
      rel = ptr->get_children( count );
      for (j = 0; j < count; ++j)
        remove_parent_meshset( rel[j],handle );
    }
  }

    // now delete the entity
  result = sequence_manager()->delete_entity(handle);
  if (MB_SUCCESS != result) 
    return result;
  
  return MB_SUCCESS;
}

//! deletes an entity vector
MBErrorCode MBCore::delete_entities(const MBEntityHandle *entities,
                                      const int num_entities)
{
  // If reference counting is enabled, make sure all entities 
  // are unreferenced, assuming all other entities in the list
  // are also deleted.
#ifdef MOAB_WITH_REFCOUNT
  std::vector<MBTag> handle_tags;
  tag_server()->get_tags( MB_TYPE_HANDLE, handle_tags );
  
  int ii, jj;
    // before deleting each entity, decrement the reference count on
    // everything that entity references.
  for (ii = 0; ii < num_entities; ++ii) 
    if (MB_SUCCESS != decrement_all_referenced_entities( entities[ii], handle_tags )) 
      break;
    // make sure everything that is to be deleted has a reference count 
    // of 1 (only referenced by the containing sequence).
  if (ii == num_entities) 
    for (jj = 0; jj < num_entities; ++jj) 
      if (get_reference_count( entities[jj] ) != 1)
        break;
    // if we can't delete the entities, but the reference counts back
    // to what they were before and return failure.
  if (ii < num_entities || jj < num_entities) {
    while (ii > 0)
      increment_all_referenced_entities( entities[--ii], handle_tags );
    return MB_FAILURE;
  }
  
#endif  

  MBErrorCode result = MB_SUCCESS, temp_result;
  for (int i = 0; i < num_entities; i++) {
    temp_result = delete_entity( entities[i] );
    if (MB_SUCCESS != temp_result)
      result = temp_result;
  }
    
  return result;
}


//! deletes an entity range
MBErrorCode MBCore::delete_entities(const MBRange &range)
{
  // If reference counting is enabled, make sure all entities 
  // are unreferenced, assuming all other entities in the list
  // are also deleted.
#ifdef MOAB_WITH_REFCOUNT
  std::vector<MBTag> handle_tags;
  tag_server()->get_tags( MB_TYPE_HANDLE, handle_tags );
  
  MBRange::const_iterator ii, jj;
    // before deleting each entity, decrement the reference count on
    // everything that entity references.
  for (ii = range.begin(); ii != range.end(); ++ii) 
    if (MB_SUCCESS != decrement_all_referenced_entities( *ii, handle_tags )) 
      break;
    // make sure everything that is to be deleted has a reference count 
    // of 1 (only referenced by the containing sequence).
  if (ii == range.end()) 
    for (jj = range.begin(); jj != range.end(); ++jj) 
      if (get_reference_count( *jj ) != 1)
        break;
    // if we can't delete the entities, but the reference counts back
    // to what they were before and return failure.
  if (ii != range.end() || jj != range.end()) {
    while (ii != range.begin())
      increment_all_referenced_entities( *--ii, handle_tags );
    return MB_FAILURE;
  }
#endif  
  MBErrorCode result = MB_SUCCESS, rval;
  for (MBRange::iterator i = range.begin(); i != range.end(); ++i)
    if (MB_SUCCESS != (rval = delete_entity( *i )))
      result = rval;
  return rval;
}

MBErrorCode MBCore::list_entities(const MBEntityHandle *entities,
                                    const int num_entities) const
{
  MBRange temp_range;
  MBErrorCode result = MB_SUCCESS;
  if (NULL == entities && num_entities <= 0) {
      // just list the numbers of each entity type
    int num_ents;
    std::cout << std::endl;
    std::cout << "Number of entities per type: " << std::endl;
    for (MBEntityType this_type = MBVERTEX; this_type < MBMAXTYPE; this_type++) {
      result = get_number_entities_by_type(0, this_type, num_ents);
      std::cout << MBCN::EntityTypeName(this_type) << ": " << num_ents << std::endl;
    }
    std::cout << std::endl;

      // if negative num_entities, list the set hierarchy too
    if (0 > num_entities) {
      MBRange sets;
      result = this->get_entities_by_type(0, MBENTITYSET, sets);
      if (MB_SUCCESS != result) return result;
      for (MBRange::iterator rit = sets.begin(); rit != sets.end(); rit++) {
        this->print(*rit, "", false);
        result = this->get_number_entities_by_handle(*rit, num_ents);
        std::cout << "(" << num_ents << " total entities)" << std::endl;
      }
    }
    
    return MB_SUCCESS;
  }
      
  else if (NULL == entities) {

      // list all entities of all types
    std::cout << std::endl;
    for (MBEntityType this_type = MBVERTEX; this_type < MBMAXTYPE; this_type++) {
      result = get_entities_by_type(0, this_type, temp_range);
    }

    return list_entities(temp_range);
  }

  else {
    MBErrorCode tmp_result;
    for (int i = 0; i < num_entities; i++) {
      MBEntityType this_type = TYPE_FROM_HANDLE(entities[i]);
      std::cout << MBCN::EntityTypeName(this_type) << " " 
                << ID_FROM_HANDLE(entities[i]) << ":" << endl;

      tmp_result = (const_cast<MBCore*>(this))->list_entity(entities[i]);
      if (MB_SUCCESS != tmp_result) result = tmp_result;
    }
  }

  return result;
}

MBErrorCode MBCore::list_entities(const MBRange &temp_range) const
{
  MBErrorCode result = MB_SUCCESS, tmp_result;
  
  for (MBRange::const_iterator rit = temp_range.begin(); rit != temp_range.end(); rit++) {
    MBEntityType this_type = TYPE_FROM_HANDLE(*rit);
    std::cout << MBCN::EntityTypeName(this_type) << " " << ID_FROM_HANDLE(*rit) << ":" << endl;

    tmp_result = (const_cast<MBCore*>(this))->list_entity(*rit);
    if (MB_SUCCESS != tmp_result) result = tmp_result;
  }
    
  return result;
}
  
MBErrorCode MBCore::list_entity(const MBEntityHandle entity) const
{
  MBErrorCode result;
  MBHandleVec adj_vec;

  if (!is_valid(entity)) {
    std::cout << "(invalid)" << std::endl;
    return MB_SUCCESS;
  }

  if (0 != globalIdTag) {
    int dum;
    result = tag_get_data(globalIdTag, &entity, 1, &dum);
    if (MB_SUCCESS == result)
      std::cout << "Global id = " << dum << std::endl;
  }
  
    // list entity
  MBEntityType this_type = TYPE_FROM_HANDLE(entity);
  if (this_type == MBVERTEX) {
    double coords[3];
    result = get_coords(&(entity), 1, coords);
    if (MB_SUCCESS != result) return result;
    std::cout << "Coordinates: (" << coords[0] << ", " << coords[1] << ", " << coords[2] 
              << ")" << std::endl;
  }
  else if (this_type == MBENTITYSET)
    this->print(entity, "");
    
  std::cout << "  Adjacencies:" << std::endl;
  bool some = false;
  int multiple = 0;
  for (int dim = 0; dim <= 3; dim++) {
    if (dim == MBCN::Dimension(this_type)) continue;
    adj_vec.clear();
      // use const_cast here 'cuz we're in a const function and we're passing 'false' for
      // create_if_missing, so we know we won't change anything
    result = (const_cast<MBCore*>(this))->get_adjacencies(&(entity), 1, dim, false, adj_vec);
    if (MB_FAILURE == result) continue;
    for (MBHandleVec::iterator adj_it = adj_vec.begin(); adj_it != adj_vec.end(); adj_it++) {
      if (adj_it != adj_vec.begin()) std::cout << ", ";
      else std::cout << "   ";
      std::cout << MBCN::EntityTypeName(TYPE_FROM_HANDLE(*adj_it)) << " " << ID_FROM_HANDLE(*adj_it);
    }
    if (!adj_vec.empty()) {
      std::cout << std::endl;
      some = true;
    }
    if (MB_MULTIPLE_ENTITIES_FOUND == result)
      multiple += dim;
  }
  if (!some) std::cout << "(none)" << std::endl;
  const MBEntityHandle *explicit_adjs;
  int num_exp;
  aEntityFactory->get_adjacencies(entity, explicit_adjs, num_exp);
  if (NULL != explicit_adjs && 0 != num_exp) {
    std::cout << "  Explicit adjacencies: ";
    for (int i = 0; i < num_exp; i++) {
      if (i != 0) std::cout << ", ";
      std::cout << MBCN::EntityTypeName(TYPE_FROM_HANDLE(explicit_adjs[i])) << " " 
                << ID_FROM_HANDLE(explicit_adjs[i]);
    }
    std::cout << std::endl;
  }
  if (multiple != 0)
    std::cout << "   (MULTIPLE = " << multiple << ")" << std::endl;

  std::cout << std::endl;

  return MB_SUCCESS;
}

MBErrorCode MBCore::convert_entities( const MBEntityHandle meshset, 
                                        const bool mid_side, const bool mid_face, const bool mid_volume,
                                        MBInterface::HONodeAddedRemoved* function_object )
{
  HigherOrderFactory fact(this, function_object);
  return fact.convert(meshset, mid_side, mid_face, mid_volume);
}

  //! function to get the side number given two elements; returns
  //! MB_FAILURE if child not related to parent; does *not* create adjacencies
  //! between parent and child
MBErrorCode MBCore::side_number(const MBEntityHandle parent,
                                  const MBEntityHandle child,
                                  int &side_number,
                                  int &sense,
                                  int &offset) const
{
    // get the connectivity of parent and child
  const MBEntityHandle *parent_conn, *child_conn;
  int num_parent_vertices, num_child_vertices;
  MBErrorCode result = get_connectivity(parent, parent_conn, num_parent_vertices, true);
  if (MB_SUCCESS != result) return result;
  result = get_connectivity(child, child_conn, num_child_vertices, true);
  if (MB_SUCCESS != result) return result;

    // call handle vector-based function
  if (TYPE_FROM_HANDLE(parent) != MBPOLYGON &&
      TYPE_FROM_HANDLE(parent) != MBPOLYHEDRON) {

      // find indices into parent_conn for each entry in child_conn
    int child_conn_indices[10];
    assert((unsigned)num_child_vertices <= sizeof(child_conn_indices)/sizeof(child_conn_indices[0]));
    for (int i = 0; i < num_child_vertices; ++i) {
      child_conn_indices[i] = std::find( parent_conn,
        parent_conn + num_parent_vertices, child_conn[i] ) - parent_conn;
      if (child_conn_indices[i] >= num_parent_vertices) {
        side_number = -1;
        return MB_SUCCESS;
      }
    }
    
    int temp_result = MBCN::SideNumber(TYPE_FROM_HANDLE(parent),
                                       child_conn_indices, num_child_vertices, 
                                       MBCN::Dimension(TYPE_FROM_HANDLE(child)), 
                                       side_number, sense, offset);
    return (0 == temp_result ? MB_SUCCESS : MB_FAILURE);
  }
  else if (TYPE_FROM_HANDLE(parent) == MBPOLYGON) {
      // find location of 1st vertex
    const MBEntityHandle *first_v = std::find(parent_conn, parent_conn+num_parent_vertices,
                                              child_conn[0]);
    if (first_v == parent_conn+num_parent_vertices) return MB_ENTITY_NOT_FOUND;
    side_number = first_v - parent_conn;
    offset = side_number;
    if (TYPE_FROM_HANDLE(child) == MBVERTEX) {
      sense = 0;
      return MB_SUCCESS;
    }
    else if (TYPE_FROM_HANDLE(child) == MBPOLYGON) {
      bool match = MBCN::ConnectivityMatch(parent_conn, child_conn,
                                           num_parent_vertices,
                                           sense, offset);
      side_number = 0;
      if (match) return MB_SUCCESS;
      else return MB_ENTITY_NOT_FOUND;
    }
    else if (TYPE_FROM_HANDLE(child) == MBEDGE) {
      if (parent_conn[(side_number+1)%num_parent_vertices] == child_conn[1])
        sense = 1;
      else if (parent_conn[(side_number+num_parent_vertices-1)%num_parent_vertices] ==
               child_conn[1])
        sense = -1;
      return MB_SUCCESS;
    }
  }
  
  return MB_FAILURE;
}

  //! given an entity and the connectivity and type of one of its subfacets, find the
  //! high order node on that subfacet, if any
MBErrorCode MBCore::high_order_node(const MBEntityHandle parent_handle,
                                      const MBEntityHandle *subfacet_conn,
                                      const MBEntityType subfacet_type,
                                      MBEntityHandle &high_order_node) const
{
  high_order_node = 0;

  MBEntityType parent_type = TYPE_FROM_HANDLE(parent_handle);

    // get the parent's connectivity
  const MBEntityHandle *parent_conn;
  int num_parent_vertices;
  MBErrorCode result = get_connectivity(parent_handle, parent_conn, 
                                         num_parent_vertices, false);
  if (result != MB_SUCCESS) return result;

    // find whether this entity has ho nodes
  int mid_nodes[4];
  MBCN::HasMidNodes(parent_type, num_parent_vertices, mid_nodes);

    // check whether this entity has mid nodes on this dimension subfacet; 
    // use dimension-1 because vertices don't have mid nodes
  if (!mid_nodes[MBCN::Dimension(subfacet_type)]) return MB_SUCCESS;

    // ok, we have mid nodes; now must compute expected index in connectivity array; 
    // ho nodes stored for edges, faces then entity

    // offset starts with # corner vertices
  int offset = MBCN::VerticesPerEntity(parent_type);
  int i;

  for (i = 0; i < MBCN::Dimension(subfacet_type)-1; i++)
      // for each dimension lower than that of the subfacet we're looking for, 
      // if this entity has midnodes in that dimension, increment offset by # 
      // of subfacets of that dimension; use dimension-1 in loop because 
      // canon numbering table only has 2 positions, for edges and faces;
    if (mid_nodes[i+1]) offset += MBCN::mConnectivityMap[parent_type][i].num_sub_elements;

    // now add the index of this subfacet; only need to if it's not the highest dimension
  if (subfacet_type != parent_type) {

      // find indices into parent_conn for each entry in child_conn
    unsigned subfacet_size = MBCN::VerticesPerEntity(subfacet_type);
    int subfacet_indices[10];
    assert(subfacet_size <= sizeof(subfacet_indices)/sizeof(subfacet_indices[0]));
    for (unsigned i = 0; i < subfacet_size; ++i) {
      subfacet_indices[i] = std::find( parent_conn,
        parent_conn + num_parent_vertices, subfacet_conn[i] ) - parent_conn;
      if (subfacet_indices[i] >= num_parent_vertices) {
        return MB_FAILURE;
      }
    }

    int dum, side_no, temp_offset;
    int temp_result = 
      MBCN::SideNumber(  parent_type, subfacet_indices, 
                         subfacet_size, subfacet_type,
                         side_no, dum, temp_offset);
    if(temp_result != 0) return MB_FAILURE;

    offset += side_no;
  }

    // offset shouldn't be off the end of the connectivity vector
  if (offset >= num_parent_vertices) return MB_INDEX_OUT_OF_RANGE;

  high_order_node = parent_conn[offset];

  return MB_SUCCESS;
}

  //! given an entity and a target dimension & side number, get that entity
MBErrorCode MBCore::side_element(const MBEntityHandle source_entity,
                                   const int dim, 
                                   const int side_number,
                                   MBEntityHandle &target_entity) const
{
    // get a handle on the connectivity
  const MBEntityHandle *verts;
  int num_verts;
  MBErrorCode result = get_connectivity(source_entity, verts, num_verts);
  if (MB_SUCCESS != result) return result;
  
    // get the vertices comprising the target entity
  MBRange side_verts, target_ents;
  const MBEntityType source_type = TYPE_FROM_HANDLE(source_entity);
    // first get the indices
  std::vector<int> vertex_indices;
  int temp_result = 
    MBCN::AdjacentSubEntities(source_type, &side_number, 1, dim, 0, vertex_indices);
  if (0 != temp_result) return MB_FAILURE;
    // now get the actual vertices
  for (unsigned int i = 0; i < vertex_indices.size(); i++)
    side_verts.insert(verts[vertex_indices[i]]);
  
    // now look for an entity of the correct type
    // use const_cast here 'cuz we're in a const function and we're passing 'false' for
    // create_if_missing, so we know we won't change anything
  result = (const_cast<MBCore*>(this))->get_adjacencies(side_verts, dim, false, target_ents);
  if (MB_SUCCESS != result && MB_MULTIPLE_ENTITIES_FOUND != result) return result;
  
  if (!target_ents.empty() &&
         TYPE_FROM_HANDLE(*(target_ents.begin())) != 
         MBCN::mConnectivityMap[source_type][dim-1].target_type[side_number])
    return MB_ENTITY_NOT_FOUND;

  if (!target_ents.empty()) target_entity = *(target_ents.begin());
  
  return result;
}

//-------------------------MBSet Functions---------------------//

MBErrorCode MBCore::create_meshset(const unsigned int options, 
                                   MBEntityHandle &ms_handle,
                                   int ,
                                   int start_proc)
{
  if (-1 == start_proc) start_proc = procInfo.rank();
  return sequence_manager()->create_mesh_set( start_proc, options, ms_handle );
}

MBErrorCode MBCore::get_meshset_options( const MBEntityHandle ms_handle, 
                                          unsigned int& options) const
{
  MBMeshSet* set = get_mesh_set( sequence_manager(), ms_handle );
  if (!set)
    return MB_ENTITY_NOT_FOUND;
  
  options = set->flags();
  return MB_SUCCESS;
}

MBErrorCode MBCore::clear_meshset( const MBEntityHandle *ms_handles,
                                    const int num_meshsets)
{
  MBErrorCode result = MB_SUCCESS;
  for (int i = 0; i < num_meshsets; ++i) {
    MBMeshSet* set = get_mesh_set( sequence_manager(), ms_handles[i]);
    if (set)
      set->clear(ms_handles[i], a_entity_factory());
    else
      result = MB_ENTITY_NOT_FOUND;
  }

  return result;
}

MBErrorCode MBCore::clear_meshset(const MBRange &ms_handles)
{
  MBErrorCode result = MB_SUCCESS;
  for (MBRange::iterator i = ms_handles.begin(); i != ms_handles.end(); ++i) {
    MBMeshSet* set = get_mesh_set( sequence_manager(), *i);
    if (set)
      set->clear(*i, a_entity_factory());
    else
      result = MB_ENTITY_NOT_FOUND;
  }

  return result;
}

MBErrorCode MBCore::subtract_meshset(MBEntityHandle meshset1, const MBEntityHandle meshset2)
{ 
  MBMeshSet *set1 = get_mesh_set( sequence_manager(), meshset1 );
  MBMeshSet *set2 = get_mesh_set( sequence_manager(), meshset2 );
  if (!set1 || !set2)
    return MB_ENTITY_NOT_FOUND;
  
  return set1->subtract( set2, meshset1, a_entity_factory() );
}


MBErrorCode MBCore::intersect_meshset(MBEntityHandle meshset1, const MBEntityHandle meshset2)
{
  MBMeshSet *set1 = get_mesh_set( sequence_manager(), meshset1 );
  MBMeshSet *set2 = get_mesh_set( sequence_manager(), meshset2 );
  if (!set1 || !set2)
    return MB_ENTITY_NOT_FOUND;
  
  return set1->intersect( set2, meshset1, a_entity_factory() );
}

MBErrorCode MBCore::unite_meshset(MBEntityHandle meshset1, const MBEntityHandle meshset2)
{
  MBMeshSet *set1 = get_mesh_set( sequence_manager(), meshset1 );
  MBMeshSet *set2 = get_mesh_set( sequence_manager(), meshset2 );
  if (!set1 || !set2)
    return MB_ENTITY_NOT_FOUND;
  
  return set1->unite( set2, meshset1, a_entity_factory() );
}

MBErrorCode MBCore::add_entities(MBEntityHandle meshset, 
                                   const MBRange &entities)
{
  MBMeshSet* set = get_mesh_set( sequence_manager(), meshset );
  if (set)
    return set->add_entities( entities, meshset, a_entity_factory() );
  else
    return MB_ENTITY_NOT_FOUND;
}

MBErrorCode MBCore::add_entities(MBEntityHandle meshset, 
                                   const MBEntityHandle *entities,
                                   const int num_entities)
{
  MBMeshSet* set = get_mesh_set( sequence_manager(), meshset );
  if (set)
    return set->add_entities( entities, num_entities, meshset, a_entity_factory() );
  else
    return MB_ENTITY_NOT_FOUND;
}


//! remove a range of entities from a meshset
MBErrorCode MBCore::remove_entities(MBEntityHandle meshset, 
                                      const MBRange &entities)
{
  MBMeshSet* set = get_mesh_set( sequence_manager(), meshset );
  if (set)
    return set->remove_entities( entities, meshset, a_entity_factory() );
  else
    return MB_ENTITY_NOT_FOUND;
}

//! remove a vector of entities from a meshset
MBErrorCode MBCore::remove_entities( MBEntityHandle meshset, 
                                       const MBEntityHandle *entities,
                                       const int num_entities)
{
  MBMeshSet* set = get_mesh_set( sequence_manager(), meshset );
  if (set)
    return set->remove_entities( entities, num_entities, meshset, a_entity_factory() );
  else
    return MB_ENTITY_NOT_FOUND;
}

MBErrorCode MBCore::get_parent_meshsets(const MBEntityHandle meshset,
                                          std::vector<MBEntityHandle> &parents,
                                          const int num_hops) const
{
  if (0 == meshset) return MB_SUCCESS;

  MBEntitySequence *seq;
  MBErrorCode rval = sequence_manager()->find( meshset, seq );
  if (MB_SUCCESS != rval || !seq->is_valid_entity(meshset))
    return MB_ENTITY_NOT_FOUND;
  MeshSetSequence* mseq = reinterpret_cast<MeshSetSequence*>(seq);

  return mseq->get_parents( meshset, parents, num_hops );
}

MBErrorCode MBCore::get_parent_meshsets(const MBEntityHandle meshset,
                                        MBRange &parents,
                                          const int num_hops) const
{
  if (0 == meshset) return MB_SUCCESS;

  std::vector<MBEntityHandle> parent_vec;
  MBErrorCode result = get_parent_meshsets(meshset, parent_vec, num_hops);
  if (MB_SUCCESS != result) return result;
  std::sort( parent_vec.begin(), parent_vec.end() );
  std::copy(parent_vec.rbegin(), parent_vec.rend(), mb_range_inserter(parents));
  return MB_SUCCESS;
}

MBErrorCode MBCore::get_child_meshsets(const MBEntityHandle meshset,
                                         std::vector<MBEntityHandle> &children,
                                         const int num_hops) const
{
  if (0 == meshset) return MB_SUCCESS;

  MBEntitySequence *seq;
  MBErrorCode rval = sequence_manager()->find( meshset, seq );
  if (MB_SUCCESS != rval || !seq->is_valid_entity(meshset))
    return MB_ENTITY_NOT_FOUND;
  MeshSetSequence* mseq = reinterpret_cast<MeshSetSequence*>(seq);

  return mseq->get_children( meshset, children, num_hops );
}

MBErrorCode MBCore::get_child_meshsets(const MBEntityHandle meshset,
                                        MBRange &children,
                                          const int num_hops) const
{
  if (0 == meshset) return MB_SUCCESS;

  std::vector<MBEntityHandle> child_vec;
  MBErrorCode result = get_child_meshsets(meshset, child_vec, num_hops);
  if (MB_SUCCESS != result) return result;
  std::sort( child_vec.begin(), child_vec.end() );
  std::copy(child_vec.rbegin(), child_vec.rend(), mb_range_inserter(children));
  return MB_SUCCESS;
}

MBErrorCode MBCore::num_parent_meshsets(const MBEntityHandle meshset, int* number,
                                        const int num_hops) const
{
  if (0 == meshset) {
    *number = 0;
    return MB_SUCCESS;
  }

  MBEntitySequence *seq;
  MBErrorCode rval = sequence_manager()->find( meshset, seq );
  if (MB_SUCCESS != rval || !seq->is_valid_entity(meshset))
    return MB_ENTITY_NOT_FOUND;
  MeshSetSequence* mseq = reinterpret_cast<MeshSetSequence*>(seq);

  return mseq->num_parents( meshset, *number, num_hops );
}

MBErrorCode MBCore::num_child_meshsets(const MBEntityHandle meshset, int* number,
                                       const int num_hops) const
{
  if (0 == meshset) {
    *number = 0;
    return MB_SUCCESS;
  }
  
  MBEntitySequence *seq;
  MBErrorCode rval = sequence_manager()->find( meshset, seq );
  if (MB_SUCCESS != rval || !seq->is_valid_entity(meshset))
    return MB_ENTITY_NOT_FOUND;
  MeshSetSequence* mseq = reinterpret_cast<MeshSetSequence*>(seq);

  return mseq->num_children( meshset, *number, num_hops );
}


MBErrorCode MBCore::add_parent_meshset( MBEntityHandle meshset, 
                                        const MBEntityHandle parent_meshset)
{
  return add_parent_meshsets( meshset, &parent_meshset, 1 );
}

MBErrorCode MBCore::add_parent_meshsets( MBEntityHandle meshset, 
                                         const MBEntityHandle* parents,
                                         int count )
{
  MBMeshSet* set_ptr = get_mesh_set( sequence_manager(), meshset );
  if (!set_ptr)
    return MB_ENTITY_NOT_FOUND;

  for (int i = 0; i < count; ++i)
    if (!get_mesh_set( sequence_manager(), parents[i] ))
      return MB_ENTITY_NOT_FOUND;
    
#ifdef MOAB_WITH_REFCOUNT
  for (int i = 0; i < count; ++i)
    if (set_ptr->add_parent( parents[i] ))
      increment_reference_count( parents + i, 1 );
#else
  for (int i = 0; i < count; ++i)
    set_ptr->add_parent( parents[i] );
#endif
  return MB_SUCCESS;
}

MBErrorCode MBCore::add_child_meshset(MBEntityHandle meshset, 
                                        const MBEntityHandle child_meshset)
{
  return add_child_meshsets( meshset, &child_meshset, 1 );
}

MBErrorCode MBCore::add_child_meshsets( MBEntityHandle meshset, 
                                        const MBEntityHandle* children,
                                        int count )
{
  MBMeshSet* set_ptr = get_mesh_set( sequence_manager(), meshset );
  if (!set_ptr)
    return MB_ENTITY_NOT_FOUND;

  for (int i = 0; i < count; ++i)
    if (!get_mesh_set( sequence_manager(), children[i] ))
      return MB_ENTITY_NOT_FOUND;
    
#ifdef MOAB_WITH_REFCOUNT
  for (int i = 0; i < count; ++i)
    if (set_ptr->add_child( children[i] ))
      increment_reference_count( children+i, 1 );
#else
  for (int i = 0; i < count; ++i)
    set_ptr->add_child( children[i] );
#endif
  return MB_SUCCESS;
}


MBErrorCode MBCore::add_parent_child(MBEntityHandle parent, 
                                       MBEntityHandle child)
{
  MBMeshSet* parent_ptr = get_mesh_set( sequence_manager(), parent );
  MBMeshSet* child_ptr = get_mesh_set( sequence_manager(), child );
  if (!parent_ptr || !child_ptr)
    return MB_ENTITY_NOT_FOUND;

#ifdef MOAB_WITH_REFCOUNT
  if (parent_ptr->add_child( child ))
    increment_reference_count( &child, 1 );
  if (child_ptr->add_parent( parent ))
    increment_reference_count( &parent, 1 );
#else
  parent_ptr->add_child( child );
  child_ptr->add_parent( parent );
#endif
  return MB_SUCCESS;
}

MBErrorCode MBCore::remove_parent_child(MBEntityHandle parent, 
                                          MBEntityHandle child)
{
  MBMeshSet* parent_ptr = get_mesh_set( sequence_manager(), parent );
  MBMeshSet* child_ptr = get_mesh_set( sequence_manager(), child );
  if (!parent_ptr || !child_ptr)
    return MB_ENTITY_NOT_FOUND;
  
#ifdef MOAB_WITH_REFCOUNT
  if (parent_ptr->remove_child( child ))
    decrement_reference_count( &child, 1 );
  if (child_ptr->remove_parent( parent ))
    decrement_reference_count( &parent, 1 );
#else
  parent_ptr->remove_child( child );
  child_ptr->remove_parent( parent );
#endif
  return MB_SUCCESS;
}


MBErrorCode MBCore::remove_parent_meshset(MBEntityHandle meshset, 
                                            const MBEntityHandle parent_meshset)
{
  MBMeshSet* set_ptr = get_mesh_set( sequence_manager(), meshset );
  if (!set_ptr)
    return MB_ENTITY_NOT_FOUND;
#ifdef MOAB_WITH_REFCOUNT
  if (set_ptr->remove_parent( parent_meshset ))
    decrement_reference_count( &parent_meshset, 1 );
#else
  set_ptr->remove_parent( parent_meshset );
#endif
  return MB_SUCCESS;
}

MBErrorCode MBCore::remove_child_meshset(MBEntityHandle meshset, 
                                           const MBEntityHandle child_meshset)
{
  MBMeshSet* set_ptr = get_mesh_set( sequence_manager(), meshset );
  if (!set_ptr)
    return MB_ENTITY_NOT_FOUND;
#ifdef MOAB_WITH_REFCOUNT
  if (set_ptr->remove_child( child_meshset ))
    decrement_reference_count( &child_meshset, 1 );
#else
  set_ptr->remove_child( child_meshset );
#endif
  return MB_SUCCESS;
}


MBErrorCode MBCore::get_last_error(std::string& info) const
{
  return mError->get_last_error(info);
}

std::string MBCore::get_error_string(const MBErrorCode code) const 
{
  return errorStrings[code];
}

void MBCore::print(const MBEntityHandle ms_handle, const char *prefix,
                   bool first_call) const
{
    // get the entities
  MBRange entities;
  
  if (0 != ms_handle) {
    get_entities_by_handle( ms_handle, entities );
    if (!first_call)
      std::cout << prefix << "MBENTITYSET " << ID_FROM_HANDLE(ms_handle) 
                << std::endl;
  }
  else {
    get_entities_by_dimension(0, 3, entities);
    if (entities.empty()) get_entities_by_dimension(0, 2, entities);
    if (entities.empty()) get_entities_by_dimension(0, 1, entities);
    get_entities_by_dimension(0, 0, entities);
    get_entities_by_type(0, MBENTITYSET, entities);
    std::cout << prefix << "--MB: " << std::endl;
  }
    
  std::string indent_prefix = prefix;
  indent_prefix += "  ";
  entities.print(indent_prefix.c_str());

  if (!first_call || !ms_handle) return;
  
    // print parent/children
  MBRange temp;
  this->get_parent_meshsets(ms_handle, temp);
  std::cout << "  Parent sets: ";
  if (temp.empty()) std::cout << "(none)" << std::endl;
  else {
    for (MBRange::iterator rit = temp.begin(); rit != temp.end(); rit++) {
      if (rit != temp.begin()) std::cout << ", ";
      std::cout << ID_FROM_HANDLE(*rit);
    }
    std::cout << std::endl;
  }

  temp.clear();
  this->get_child_meshsets(ms_handle, temp);
  std::cout << "  Child sets: ";
  if (temp.empty()) std::cout << "(none)" << std::endl;
  else {
    for (MBRange::iterator rit = temp.begin(); rit != temp.end(); rit++) {
      if (rit != temp.begin()) std::cout << ", ";
      std::cout << ID_FROM_HANDLE(*rit);
    }
    std::cout << std::endl;
  }

    // print all sparse tags
  std::vector<MBTag> set_tags;
  MBErrorCode result = this->tag_get_tags_on_entity(ms_handle, set_tags);
  std::cout << indent_prefix << "Sparse tags:" << std::endl;
  indent_prefix += "  ";
  
  for (std::vector<MBTag>::iterator vit = set_tags.begin(); 
       vit != set_tags.end(); vit++) {
    MBTagType this_type;
    result = this->tag_get_type(*vit, this_type);
    if (MB_SUCCESS != result || MB_TAG_SPARSE != this_type) continue;
    MBDataType this_data_type;
    result = this->tag_get_data_type(*vit, this_data_type);
    int this_size;
    result = this->tag_get_size(*vit, this_size);
    if (MB_SUCCESS != result || (int) sizeof(double) < this_size) continue;
      // use double since this is largest single-valued tag
    double this_val;
    result = this->tag_get_data(*vit, &ms_handle, 1, &this_val);
    if (MB_SUCCESS != result) continue;
    std::string tag_name;
    result = this->tag_get_name(*vit, tag_name);
    if (MB_SUCCESS != result) continue;
    switch (this_data_type) {
      case MB_TYPE_INTEGER:
        std::cout << indent_prefix << tag_name << " = " 
                  << *((int*)&this_val) << std::endl;
        break;
      case MB_TYPE_DOUBLE:
        std::cout << indent_prefix << tag_name << " = " 
                  << this_val << std::endl;
        break;
      case MB_TYPE_HANDLE:
        std::cout << indent_prefix << tag_name << " = " 
                  << *((MBEntityID*)&this_val) << std::endl;
        break;
      case MB_TYPE_BIT:
      case MB_TYPE_OPAQUE:
        break;
    }
  }
}

MBErrorCode MBCore::check_adjacencies() 
{
    // run through all entities, checking adjacencies and reverse-evaluating them
  MBRange all_ents;
  MBErrorCode result = get_entities_by_handle(0, all_ents);
  if (MB_SUCCESS != result) return result;
  
  MBErrorCode tmp_result;
  for (MBRange::iterator rit = all_ents.begin(); rit != all_ents.end(); rit++) {
    tmp_result = check_adjacencies(&(*rit), 1);
    if (MB_SUCCESS != tmp_result) result = tmp_result;
  }
  
  return result;
}

MBErrorCode MBCore::check_adjacencies(const MBEntityHandle *ents, int num_ents) 
{

  MBErrorCode result = MB_SUCCESS, tmp_result;
  std::ostringstream oss;
  
  for (int i = 0; i < num_ents; i++) {
    MBEntityHandle this_ent = ents[i];
    std::ostringstream ent_str;
    ent_str << MBCN::EntityTypeName(TYPE_FROM_HANDLE(this_ent)) << " "
            << ID_FROM_HANDLE(this_ent) << ": ";
    int this_dim = dimension_from_handle(this_ent);

    if (!is_valid(this_ent)) {
      std::cerr << ent_str.str()
                << "Not a valid entity." << std::endl;
      result = MB_FAILURE;
    }

    else {
      if (TYPE_FROM_HANDLE(this_ent) == MBENTITYSET) continue;
      
        // get adjacencies for this entity
      MBRange adjs;
      for (int dim = 0; dim <= 3; dim++) {
        if (dim == this_dim) continue;
        tmp_result = get_adjacencies(&this_ent, 1, dim, false, adjs, MBInterface::UNION);
        if (MB_SUCCESS != tmp_result) {
          oss << ent_str.str()
              << "Failed to get adjacencies for dimension " << dim << "." << std::endl;
          result = tmp_result;
        }
      }
      if (!oss.str().empty()) {
        std::cerr << oss.str();
        oss.str("");
      }

        // now check and reverse-evaluate them
      for (MBRange::iterator rit = adjs.begin(); rit != adjs.end(); rit++) {
        MBEntitySequence* seq = 0;
        tmp_result = sequence_manager()->find(*rit, seq);
        if(seq == 0 || tmp_result != MB_SUCCESS || !seq->is_valid_entity(*rit)) {
          oss << ent_str.str() << 
            "Adjacent entity " << MBCN::EntityTypeName(TYPE_FROM_HANDLE(*rit)) << " "
              << ID_FROM_HANDLE(*rit) << " is invalid." << std::endl;
          result = tmp_result;
        }
        else {
          MBRange rev_adjs;
          tmp_result = get_adjacencies(&(*rit), 1, this_dim, false, rev_adjs);
          if (MB_SUCCESS != tmp_result) {
            oss << ent_str.str() 
                << "Failed to get reverse adjacency from " 
                << MBCN::EntityTypeName(TYPE_FROM_HANDLE(*rit)) << " "
                << ID_FROM_HANDLE(*rit);
            if (MB_MULTIPLE_ENTITIES_FOUND == tmp_result)
              oss << " (MULTIPLE)" << std::endl;
            else oss << " (" << tmp_result << ")" << std::endl;
            result = tmp_result;
          }
          else if (rev_adjs.find(this_ent) == rev_adjs.end()) {
            oss << ent_str.str() 
                << "Failed to find adjacency to this entity from " 
                << MBCN::EntityTypeName(TYPE_FROM_HANDLE(*rit)) << " "
                << ID_FROM_HANDLE(*rit) << "." << std::endl;
            result = tmp_result;
          }
        }
        if (!oss.str().empty()) {
          std::cerr << oss.str();
          oss.str("");
        }
      }
    }
  }
  
  return MB_SUCCESS;
}

bool MBCore::is_valid(const MBEntityHandle this_ent) const
{
  MBEntitySequence* seq = 0;
  MBErrorCode result = sequence_manager()->find(this_ent, seq);
  if(seq == 0 || result != MB_SUCCESS || !seq->is_valid_entity(this_ent))
    return false;
  else return true;
}

static unsigned long get_num_entities_with_tag( TagServer* ts, 
                                                MBTag tag,
                                                const MBRange& entities )
{
  if (entities.empty())
    return 0;
  
  int tmp;
  unsigned long total = 0;
  MBEntityType t = TYPE_FROM_HANDLE( entities.front() );
  MBEntityType e = TYPE_FROM_HANDLE( entities.back() );
  ++e;
  for (; t != e; ++t) {
    tmp = 0;
    if (MB_SUCCESS == ts->get_number_entities( entities, tag, t, tmp ))
      total += tmp;
  }
  
  return total;
}

void MBCore::estimated_memory_use_internal( const MBRange* ents,
                                  unsigned long* total_storage,
                                  unsigned long* total_amortized_storage,
                                  unsigned long* entity_storage,
                                  unsigned long* amortized_entity_storage,
                                  unsigned long* adjacency_storage,
                                  unsigned long* amortized_adjacency_storage,
                                  const MBTag* tag_array,
                                  unsigned num_tags,
                                  unsigned long* tag_storage,
                                  unsigned long* amortized_tag_storage )
{
    // Figure out which values we need to calulate
  unsigned long i_entity_storage,    ia_entity_storage, 
                i_adjacency_storage, ia_adjacency_storage, 
                i_tag_storage,       ia_tag_storage;
  unsigned long *total_tag_storage = 0, 
                *amortized_total_tag_storage =0;
  if (!tag_array) {
    total_tag_storage = tag_storage;
    amortized_total_tag_storage = amortized_tag_storage;
  }
  if (total_storage || total_amortized_storage) {
    if (!entity_storage)
      entity_storage = &i_entity_storage;
    if (!amortized_entity_storage)
      amortized_entity_storage = &ia_entity_storage;
    if (!adjacency_storage)
      adjacency_storage = &i_adjacency_storage;
    if (!amortized_adjacency_storage)
      amortized_adjacency_storage = &ia_adjacency_storage;
  }
  else {
    if (entity_storage || amortized_entity_storage) {
      if (!amortized_entity_storage)
        amortized_entity_storage = &ia_entity_storage;
      else if (!entity_storage)
        entity_storage = &i_entity_storage;
    }
    if (adjacency_storage || amortized_adjacency_storage) {
      if (!amortized_adjacency_storage)
        amortized_adjacency_storage = &ia_adjacency_storage;
      else if (!adjacency_storage)
        adjacency_storage = &i_adjacency_storage;
    }
  }
  if (!total_tag_storage && total_storage)
    total_tag_storage = &i_tag_storage;
  if (!amortized_total_tag_storage && total_amortized_storage)
    amortized_total_tag_storage = &ia_tag_storage;
    
    // get entity storage
  if (amortized_entity_storage) {
    if (ents)
      sequenceManager->get_memory_use( *ents, *entity_storage, *amortized_entity_storage );
    else
      sequenceManager->get_memory_use( *entity_storage, *amortized_entity_storage );
  }
  
    // get adjacency storage
  if (amortized_adjacency_storage) {
    if (ents)
      aEntityFactory->get_memory_use( *ents, *adjacency_storage, *amortized_adjacency_storage );
    else
      aEntityFactory->get_memory_use( *adjacency_storage, *amortized_adjacency_storage );
  }
  
    // get storage for requested list of tags
  if (tag_array) {
    for (unsigned i = 0; i < num_tags; ++i) {
      unsigned long total, per_ent, count;
      tagServer->get_memory_use( tag_array[i], total, per_ent );
      
      if (ents) {
        count = get_num_entities_with_tag( tagServer, tag_array[i], *ents );
        if (tag_storage)
          tag_storage[i] = count * per_ent;
        if (amortized_tag_storage) {
          tagServer->get_number_entities( tag_array[i], per_ent );
          if (per_ent)
            amortized_tag_storage[i] = (unsigned long)((double)total * count / per_ent);
        }
      }
      else {
        if (tag_storage) {
          tagServer->get_number_entities( tag_array[i], count );
          tag_storage[i] = count * per_ent;
        }
        if (amortized_tag_storage)
          amortized_tag_storage[i] = total;
      }
    }
  }
  
    // get storage for all tags
  if (total_tag_storage || amortized_total_tag_storage) {
    if (amortized_total_tag_storage)
      *amortized_total_tag_storage = 0;
    if (total_tag_storage)
      *total_tag_storage =0;
      
    std::vector<MBTag> tags;
    tag_get_tags( tags );
    for (unsigned i = 0; i < tags.size(); ++i) {
      unsigned long total, per_ent, count;
      tagServer->get_memory_use( tags[i], total, per_ent );
      
      if (ents) {
        count = get_num_entities_with_tag( tagServer, tags[i], *ents );
        if (total_tag_storage)
          *total_tag_storage += count * per_ent;
        if (amortized_total_tag_storage) {
          tagServer->get_number_entities( tags[i], per_ent );
          if (per_ent)
            *amortized_total_tag_storage += (unsigned long)((double)total * count / per_ent);
        }
      }
      else {
        if (total_tag_storage) {
          tagServer->get_number_entities( tags[i], count );
          *total_tag_storage += count * per_ent;
        }
        if (amortized_total_tag_storage)
          *amortized_total_tag_storage += total;
      }
    }
  }
  
    // calculate totals
  if (total_storage)
    *total_storage = *entity_storage + *adjacency_storage + *total_tag_storage;
  
  if (total_amortized_storage)
    *total_amortized_storage = *amortized_entity_storage 
                             + *amortized_adjacency_storage
                             + *amortized_total_tag_storage;
}


void  MBCore::estimated_memory_use( const MBEntityHandle* ent_array,
                                    unsigned long num_ents,
                                    unsigned long* total_storage,
                                    unsigned long* total_amortized_storage,
                                    unsigned long* entity_storage,
                                    unsigned long* amortized_entity_storage,
                                    unsigned long* adjacency_storage,
                                    unsigned long* amortized_adjacency_storage,
                                    const MBTag* tag_array,
                                    unsigned num_tags,
                                    unsigned long* tag_storage,
                                    unsigned long* amortized_tag_storage ) 
{
  MBRange range;
  
    // If non-empty entity list, call range version of function
  if (ent_array) {
    if (num_ents > 20) {
      std::vector<MBEntityHandle> list(num_ents);
      std::copy(ent_array, ent_array+num_ents, list.begin());
      std::sort( list.begin(), list.end() );
      MBRange::iterator j = range.begin();
      for (std::vector<MBEntityHandle>::reverse_iterator i = list.rbegin(); i != list.rend(); ++i)
        j = range.insert( j, *i, *i );
    }
    else {
      std::copy( ent_array, ent_array + num_ents, mb_range_inserter(range) );
    }
  }
  
  estimated_memory_use_internal( ent_array ? &range : 0,
                         total_storage,     total_amortized_storage,
                         entity_storage,    amortized_entity_storage,
                         adjacency_storage, amortized_adjacency_storage,
                         tag_array,         num_tags,
                         tag_storage,       amortized_tag_storage );
}

void MBCore::estimated_memory_use( const MBRange& ents,
                                   unsigned long* total_storage,
                                   unsigned long* total_amortized_storage,
                                   unsigned long* entity_storage,
                                   unsigned long* amortized_entity_storage,
                                   unsigned long* adjacency_storage,
                                   unsigned long* amortized_adjacency_storage,
                                   const MBTag* tag_array,
                                   unsigned num_tags,
                                   unsigned long* tag_storage,
                                   unsigned long* amortized_tag_storage )
{
  estimated_memory_use_internal( &ents,
                         total_storage,     total_amortized_storage,
                         entity_storage,    amortized_entity_storage,
                         adjacency_storage, amortized_adjacency_storage,
                         tag_array,         num_tags,
                         tag_storage,       amortized_tag_storage );
}

#ifdef MOAB_WITH_REFCOUNT

MBErrorCode MBCore::increment_reference_count( const MBEntityHandle* array, size_t count )
{
  MBErrorCode result = MB_SUCCESS, tmp;
  for (size_t i = 0; i < count; ++i) {
    MBEntitySequence* seq;
    tmp = sequence_manager()->find( array[i], seq );
    if (MB_SUCCESS != tmp)
      result = tmp;
    else if (!seq->is_valid_entity(array[i]))
      result = MB_ENTITY_NOT_FOUND;
    else
      seq->increment_reference_count( array[i] );
  }
  return result;
}

MBErrorCode MBCore::decrement_reference_count( const MBEntityHandle* array, size_t count )
{
  MBErrorCode result = MB_SUCCESS, tmp;
  for (size_t i = 0; i < count; ++i) {
    MBEntitySequence* seq;
    tmp = sequence_manager()->find( array[i], seq );
    if (MB_SUCCESS != tmp)
      result = tmp;
    else if (!seq->is_valid_entity(array[i]))
      result = MB_ENTITY_NOT_FOUND;
    else
      seq->decrement_reference_count( array[i] );
  }
  return result;
}

unsigned MBCore::get_reference_count( MBEntityHandle handle )
{
  MBEntitySequence* seq;
  MBErrorCode rval = sequence_manager()->find( handle, seq );
  return MB_SUCCESS == rval ? seq->get_reference_count(handle) : 0;
}

MBErrorCode MBCore::decrement_all_referenced_entities( MBEntityHandle handle,
                                   const std::vector<MBTag>& handle_tags )
{
  MBEntitySequence* seq;
  MBErrorCode rval = sequence_manager()->find( handle, seq );
  if (MB_SUCCESS != rval)
    return rval;
  
  seq->decrement_all_referenced_entities( handle, a_entity_factory() );
    
  rval = a_entity_factory()->decrement_referenced_entities( handle );
  if (MB_SUCCESS != rval) {
    seq->increment_all_referenced_entities( handle, a_entity_factory() );
    return rval;
  }
  
  std::vector<MBEntityHandle> tag_data;
  for (std::vector<MBTag>::const_iterator i = handle_tags.begin(); 
       i != handle_tags.end(); ++i) {
    int size;
    tag_get_size( *i, size );
    size /= sizeof(MBEntityHandle);
    tag_data.resize( size );
    rval = tag_server()->get_data( *i, &handle, 1, &tag_data[0] );
    if (MB_SUCCESS == rval) 
      decrement_reference_count( &tag_data[0], size );
  }
  return MB_SUCCESS;
}


MBErrorCode MBCore::increment_all_referenced_entities( MBEntityHandle handle,
                                   const std::vector<MBTag>& handle_tags )
{
  MBEntitySequence* seq;
  MBErrorCode rval = sequence_manager()->find( handle, seq );
  if (MB_SUCCESS != rval)
    return rval;
  
  seq->increment_all_referenced_entities( handle, a_entity_factory() );
    
  rval = a_entity_factory()->increment_referenced_entities( handle );
  if (MB_SUCCESS != rval) {
    seq->decrement_all_referenced_entities( handle, a_entity_factory() );
    return rval;
  }
  
  std::vector<MBEntityHandle> tag_data;
  for (std::vector<MBTag>::const_iterator i = handle_tags.begin(); 
       i != handle_tags.end(); ++i) {
    int size;
    tag_get_size( *i, size );
    size /= sizeof(MBEntityHandle);
    tag_data.resize( size );
    rval = tag_server()->get_data( *i, &handle, 1, &tag_data[0] );
    if (MB_SUCCESS == rval) 
      increment_reference_count( &tag_data[0], size );
  }
  return MB_SUCCESS;
}

#else

MBErrorCode MBCore::increment_reference_count( const MBEntityHandle* , size_t  )
{
  abort();
  return MB_FAILURE;
}

MBErrorCode MBCore::decrement_reference_count( const MBEntityHandle* , size_t  )
{
  abort();
  return MB_FAILURE;
}

unsigned MBCore::get_reference_count( MBEntityHandle  )
{
  abort();
  return (unsigned)-1;
}

MBErrorCode MBCore::decrement_all_referenced_entities( MBEntityHandle,
                                   const std::vector<MBTag>&  )
{
  abort();
  return MB_FAILURE;
}


MBErrorCode MBCore::increment_all_referenced_entities( MBEntityHandle ,
                                   const std::vector<MBTag>&  )
{
  abort();
  return MB_FAILURE;
}

#endif


MBErrorCode MBCore::find_all_referencing_entities( MBEntityHandle entity,
                                                   MBRange& results )
{
  MBErrorCode rval;
  MBRange entities;

    // if entity is a vertex, check all connectivity lists
  if (TYPE_FROM_HANDLE(entity) == MBVERTEX) {
    for (MBEntityType t = MBEDGE; t < MBENTITYSET; ++t) {
      if (t == MBPOLYHEDRON)
        continue;
      entities.clear();
      rval = get_entities_by_type( 0, t, entities );
      if (MB_SUCCESS != rval)
        return rval;
      for (MBRange::iterator i = entities.begin(); i != entities.end(); ++i) {
        const MBEntityHandle* conn;
        int len;
        rval = get_connectivity( *i, conn, len );
        if (MB_SUCCESS != rval)
          return rval;
        if (std::find( conn, conn+len, entity ) - conn < len)
          results.insert( *i );
      }
    }
  }
  
    // if entity type is a face, check polyhedron connectivity
  if (MBCN::Dimension( TYPE_FROM_HANDLE(entity) ) == 2) {
    entities.clear();
    rval = get_entities_by_type( 0, MBPOLYHEDRON, entities );
    if (MB_SUCCESS != rval)
      return rval;
    for (MBRange::iterator i = entities.begin(); i != entities.end(); ++i) {
      const MBEntityHandle* conn;
      int len;
      rval = get_connectivity( *i, conn, len );
      if (MB_SUCCESS != rval)
        return rval;
      if (std::find( conn, conn+len, entity ) - conn < len)
        results.insert( *i );
    }
  }
  
    // check all adjacnecy data
  entities.clear();
  rval = get_entities_by_handle( 0, entities );
  if (MB_SUCCESS != rval)
    return rval;
  
  for (MBRange::iterator i = entities.begin(); i != entities.end(); ++i) {
    const MBEntityHandle* adj;
    int len;
    rval = a_entity_factory()->get_adjacencies( *i, adj, len );
    if (MB_SUCCESS != rval)
      return rval;
    
    if (std::find( adj, adj+len, entity ) - adj < len)
      results.insert( *i );
  }
  
    // check set contents
  entities.clear();
  rval = get_entities_by_type( 0, MBENTITYSET, entities );
  if (MB_SUCCESS != rval)
    return rval;
  for (MBRange::iterator i = entities.begin(); i != entities.end(); ++i) {
    MBRange contents;
    rval = get_entities_by_handle( *i, contents );
    if (contents.find( entity ) != contents.end())
      results.insert( *i );
  }
  
    // if entity is a set, check set parent/child links
  if (TYPE_FROM_HANDLE(entity) == MBENTITYSET) {
    for (MBRange::iterator i = entities.begin(); i != entities.end(); ++i) {
      MBRange parents, children;
      rval = get_parent_meshsets( *i, parents );
      if (MB_SUCCESS != rval)
        return rval;
      rval = get_child_meshsets( *i, children );
      if (MB_SUCCESS != rval)
        return rval;
      if (parents.find(entity) != parents.end() || children.find(entity) != children.end())
        results.insert( *i );
    }
  }
  
    // check all handle tags
  std::vector<MBTag> tags;
  rval = tag_server()->get_tags( MB_TYPE_HANDLE, tags );
  while (!tags.empty()) {
    MBTag tag = tags.back();
    tags.pop_back();
    int size;
    rval = tag_get_size( tag, size );
    if (MB_SUCCESS != rval)
      return rval;
    std::vector<MBEntityHandle> tag_values(size/sizeof(MBEntityHandle));;
    
    for (MBEntityType t = MBVERTEX; t <= MBENTITYSET; ++t) {
      entities.clear();
      rval = get_entities_by_type_and_tag( 0, t, &tag, 0, 1, entities );
      if (MB_SUCCESS != rval)
        return rval;
  
      for (MBRange::iterator i = entities.begin(); i != entities.end(); ++i) {
        rval = tag_get_data( tag, &*i, 1, &tag_values[0] );
        if (MB_SUCCESS != rval)
          return rval;
        if (std::find(tag_values.begin(), tag_values.end(), entity) != tag_values.end())
          results.insert( *i );
      }
    }
  }
  
  return MB_SUCCESS;
}

    
    
