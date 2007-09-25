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

/**
 * \class ReadGmsh
 * \brief Gmsh (http://www.geuz.org/gmsh) file reader
 * \author Jason Kraftcheck
 */

#include "ReadGmsh.hpp"
#include "FileTokenizer.hpp" // for file tokenizer
#include "MBInternals.hpp"
#include "MBInterface.hpp"
#include "MBReadUtilIface.hpp"
#include "MBRange.hpp"
#include "MBTagConventions.hpp"
#include "MBParallelConventions.h"
#include "MBCN.hpp"

#include <errno.h>
#include <string.h>
#include <map>
#include <set>

MBReaderIface* ReadGmsh::factory( MBInterface* iface )
  { return new ReadGmsh(iface); }

ReadGmsh::ReadGmsh(MBInterface* impl)
    : mdbImpl(impl)
{
  mdbImpl->query_interface("MBReadUtilIface", reinterpret_cast<void**>(&readMeshIface));
}

ReadGmsh::~ReadGmsh()
{
  if (readMeshIface)
   mdbImpl->release_interface("MBReadUtilIface", reinterpret_cast<void**>(&readMeshIface));
}

// Type info indexed by type id used in file format.
const int hex_27_node_order[] =  {  
    0,  1,  2,  3,  4,  5,  6,  7,                 // corners
    8, 11, 12,  9, 13, 10, 14, 15, 16, 19, 17, 18, // edges
   24, 20, 23, 21, 22, 25,                         // faces
   26 };                                           // volume
const ReadGmsh::ElementType typemap[] = {
  { MBMAXTYPE, 0, 0 },
  { MBEDGE,    2, 0 },
  { MBTRI,     3, 0 },
  { MBQUAD,    4, 0 },
  { MBTET,     4, 0 },
  { MBHEX,     8, 0 },
  { MBPRISM,   6, 0 },
  { MBPYRAMID, 5, 0 },
  { MBEDGE,    3, 0 },
  { MBTRI,     6, 0 },
  { MBQUAD,    8, 0 },
  { MBTET,    10, 0 },
  { MBHEX,    27, hex_27_node_order },
  { MBMAXTYPE,18, 0 }, // prism w/ mid-face nodes on quads but not tris
  { MBMAXTYPE,14, 0 }, // pyramid w/ mid-face nodes on quad but not tris
  { MBMAXTYPE, 1, 0 }  // point element (0-rad sphere element?)
};
const int max_type_int = sizeof(typemap) / sizeof(typemap[0]) - 1;


MBErrorCode ReadGmsh::load_file( const char* filename, 
                                 MBEntityHandle& file_set,
                                 const FileOptions& ,
                                 const int* blocks,
                                 const int num_blocks )
{
  mCurrentMeshHandle = 0;
  const MBErrorCode result = load_file_impl( filename, blocks, num_blocks );
  
    // If file read has failed, destroy anything that was
    // created during the read.
  if (MB_SUCCESS != result && mCurrentMeshHandle)
  {
    MBRange entities;
    mdbImpl->get_entities_by_handle( mCurrentMeshHandle, entities );
    entities.insert( mCurrentMeshHandle );
    mdbImpl->delete_entities( entities );
    mCurrentMeshHandle = 0;
  }
  
  file_set = mCurrentMeshHandle;
  return result;
}


MBErrorCode ReadGmsh::load_file_impl( const char* filename, 
                                      const int* material_set_list,
                                      const int num_material_sets )
{
  geomSets.clear();
  MBErrorCode result = mdbImpl->tag_get_handle( GLOBAL_ID_TAG_NAME, globalId );
  if (MB_TAG_NOT_FOUND == result)
    result = mdbImpl->tag_create( GLOBAL_ID_TAG_NAME,
                                  sizeof(int), MB_TAG_SPARSE,
                                  MB_TYPE_INTEGER, globalId, 0 );
  if (MB_SUCCESS != result)
    return result;
    
    // Create set for more convienient check for material set ids
  std::set<int> blocks;
  for (const int* mat_set_end = material_set_list + num_material_sets;
       material_set_list != mat_set_end; ++material_set_list)
    blocks.insert( *material_set_list );
  
    // Map of ID->handle for nodes
  std::map<long,MBEntityHandle> node_id_map;
  int data_size = 8;
  
    // Open file and hand off pointer to tokenizer
  FILE* file_ptr = fopen( filename, "r" );
  if (!file_ptr)
  {
    readMeshIface->report_error( "%s: %s\n", filename, strerror(errno) );
    return MB_FILE_DOES_NOT_EXIST;
  }
  FileTokenizer tokens( file_ptr, readMeshIface );

    // Determine file format version
  const char* const start_tokens[] = { "$NOD", "$MeshFormat", 0 };
  int format_version = tokens.match_token( start_tokens );
  if (!format_version)
    return MB_FILE_DOES_NOT_EXIST;
  
    // If version 2.0, read additional header info
  if (2 == format_version)
  {
    double version;
    if (!tokens.get_doubles( 1, &version ))
      return MB_FILE_WRITE_ERROR;
    
    if (version != 2.0) {
      readMeshIface->report_error( "%s: unknown format version: %f\n", filename, version );
      return MB_FILE_DOES_NOT_EXIST;
    }
    
    int file_format;
    if (!tokens.get_integers( 1, &file_format ) ||
        !tokens.get_integers( 1, &data_size )   ||
        !tokens.match_token( "$EndMeshFormat" ) ||
        !tokens.match_token( "$Nodes") ) 
      return MB_FILE_WRITE_ERROR;
  }
  
    // read number of nodes
  long num_nodes;
  if (!tokens.get_long_ints( 1, &num_nodes ))
    return MB_FILE_WRITE_ERROR;

    // make a meshset for this mesh
  result = mdbImpl->create_meshset(MESHSET_SET, mCurrentMeshHandle);
  if (MB_SUCCESS != result) return result;
  
    // allocate nodes
  std::vector<double*> coord_arrays;
  MBEntityHandle handle = 0;
  result = readMeshIface->get_node_arrays( 3, num_nodes, MB_START_ID, 
                                           readMeshIface->parallel_rank(),
                                           handle, coord_arrays );
  if (MB_SUCCESS != result)
    return result;
  
    // put nodes in set of all loaded entities
  MBRange node_handles( handle, handle + num_nodes - 1 );
  result = mdbImpl->add_entities( mCurrentMeshHandle, node_handles );
  if (MB_SUCCESS != result)
    return result;

    // read nodes
  double *x = coord_arrays[0], 
         *y = coord_arrays[1],
         *z = coord_arrays[2];
  for( long i = 0; i < num_nodes; ++i )
  {
    long id;
    if (!tokens.get_long_ints( 1, &id ) ||
        !tokens.get_doubles( 1, x++ ) ||
        !tokens.get_doubles( 1, y++ ) ||
        !tokens.get_doubles( 1, z++ ))
      return MB_FILE_WRITE_ERROR;
    
    if (!node_id_map.insert( std::pair<long,MBEntityHandle>( id, handle++ ) ).second)
    {
      readMeshIface->report_error( "Dulicate node ID at line %d\n",
                                   tokens.line_number() );
      return MB_FILE_WRITE_ERROR;
    }
  }
  
    // create reverse map from handle to id
  std::vector<int> ids( num_nodes );
  std::vector<int>::iterator id_iter = ids.begin();
  std::vector<MBEntityHandle> handles( num_nodes );
  std::vector<MBEntityHandle>::iterator h_iter = handles.begin();
  for (std::map<long,MBEntityHandle>::iterator i = node_id_map.begin();
        i != node_id_map.end(); ++i, ++id_iter, ++h_iter)
  {
    *id_iter = i->first;
    * h_iter = i->second;
  }
    // store IDs in tags
  result = mdbImpl->tag_set_data( globalId, &handles[0], num_nodes, &ids[0] );
  if (MB_SUCCESS != result)
    return result;
  ids.clear();
  handles.clear();
  
    // get tokens signifying end of node data and start of elements
  if (!tokens.match_token( format_version == 1 ? "$ENDNOD" : "$EndNodes" ) ||
      !tokens.match_token( format_version == 1 ? "$ELM" : "$Elements" ))
    return MB_FILE_WRITE_ERROR;
  
    // get element count
  long num_elem;
  if (!tokens.get_long_ints( 1, &num_elem ))
    return MB_FILE_WRITE_ERROR;
  
    // lists of data accumulated for elements
  std::vector<MBEntityHandle> connectivity;
  std::vector<int> mat_set_list, geom_set_list, part_set_list, id_list;
    // temporary, per-element data
  std::vector<int> int_data(5), tag_data(2);
  std::vector<long> tmp_conn;
  int curr_elem_type = 0;
  for (long i = 0; i < num_elem; ++i)
  {
      // Read element description
      // File format 1.0
    if (1 == format_version)
    {
      if (!tokens.get_integers( 5, &int_data[0] ))
        return MB_FILE_WRITE_ERROR;
      tag_data[0] = int_data[2];
      tag_data[1] = int_data[3];
      if (int_data[4] != typemap[int_data[1]].nodes)
      {
        readMeshIface->report_error( "Invalid node count for element type at line %d\n",
                                     tokens.line_number() );
        return MB_FILE_WRITE_ERROR;
      }
    }
      // File format 2.0
    else
    {
      if (!tokens.get_integers( 3, &int_data[0] ))
        return MB_FILE_WRITE_ERROR;
      tag_data.resize( int_data[2] );
      if (!tokens.get_integers( tag_data.size(), &tag_data[0] ))
        return MB_FILE_WRITE_ERROR;
    }

      // If a list of material sets was specified in the
      // argument list, skip any elements for which the
      // material set is not specified or is not in the
      // passed list.
    if (!blocks.empty() && (tag_data.empty() ||  
        blocks.find( tag_data[0] ) != blocks.end()))
      continue;
    
    
      // If the next element is not the same type as the last one,
      // create a sequence for the block of elements we've read
      // to this point (all of the same type), and clear accumulated
      // data.
    if (int_data[1] != curr_elem_type)
    {
      if (!id_list.empty())  // first iteration
      {
        result = create_elements( typemap[curr_elem_type],
                                  id_list,
                                  mat_set_list,
                                  geom_set_list,
                                  part_set_list,
                                  connectivity ) ;
        if (MB_SUCCESS != result)
          return result;
      }
      
      id_list.clear();
      mat_set_list.clear();
      geom_set_list.clear();
      part_set_list.clear();
      connectivity.clear();
      curr_elem_type = int_data[1];
      if (curr_elem_type > max_type_int ||
          typemap[curr_elem_type].mbtype == MBMAXTYPE)
      {
        readMeshIface->report_error( "Unsupported element type %d at line %d\n",
                                     curr_elem_type, tokens.line_number() );
        return MB_FILE_WRITE_ERROR;
      }
      tmp_conn.resize( typemap[curr_elem_type].nodes );
    }
    
      // Store data from element description
    id_list.push_back( int_data[0] );
    part_set_list.push_back( tag_data.size() > 2 ? tag_data[2] : 0 );
    geom_set_list.push_back( tag_data.size() > 1 ? tag_data[1] : 0 );
     mat_set_list.push_back( tag_data.size() > 0 ? tag_data[0] : 0 );

      // Get element connectivity
    if (!tokens.get_long_ints( tmp_conn.size(), &tmp_conn[0] ))
      return MB_FILE_WRITE_ERROR;

      // Convert conectivity from IDs to handles
    for (unsigned j = 0; j < tmp_conn.size(); ++j)
    {
      std::map<long,MBEntityHandle>::iterator k = node_id_map.find( tmp_conn[j] );
      if (k == node_id_map.end()) {
        readMeshIface->report_error( "Invalid node ID at line %d\n",
                                     tokens.line_number() );
        return MB_FILE_WRITE_ERROR;
      }
      connectivity.push_back( k->second );
    }
  } // for (num_nodes)
 
    // Create entity sequence for last element(s).
  if (!id_list.empty())
  {
    result = create_elements( typemap[curr_elem_type],
                              id_list,
                              mat_set_list,
                              geom_set_list,
                              part_set_list,
                              connectivity ) ;
    if (MB_SUCCESS != result)
      return result;
  }
  
    // Construct parent-child relations for geometric sets.
    // Note:  At the time this comment was written, the following
    //        function was not impelemented.
  result = create_geometric_topology();
  geomSets.clear();
  return result;
}

//! Create an element sequence
MBErrorCode ReadGmsh::create_elements( const ElementType& type,
                               const std::vector<int>& elem_ids,
                               const std::vector<int>& matl_ids,
                               const std::vector<int>& geom_ids,
                               const std::vector<int>& prtn_ids,
                               const std::vector<MBEntityHandle>& connectivity )
{
  MBErrorCode result;
  
    // Make sure input is consistent
  const unsigned long num_elem = elem_ids.size();
  const int node_per_elem = type.nodes;
  if (matl_ids.size() != num_elem ||
      geom_ids.size() != num_elem ||
      prtn_ids.size() != num_elem ||
      connectivity.size() != num_elem*node_per_elem)
    return MB_FAILURE;
  
    // Create the element sequence
  MBEntityHandle handle = 0;
  MBEntityHandle* conn_array;
  result = readMeshIface->get_element_array( num_elem, node_per_elem, type.mbtype,
                                             MB_START_ID, 
                                             readMeshIface->parallel_rank(), 
                                             handle, conn_array );
  if (MB_SUCCESS != result)
    return result;
  
    // Put newly created elements in set of all entities read from file.
  MBRange elements( handle, handle + num_elem - 1 );
  result = mdbImpl->add_entities( mCurrentMeshHandle, elements );
  if (MB_SUCCESS != result)
    return result;
  
    // Copy passed element connectivity into entity sequence data.
  if (type.node_order)
  {
    for (unsigned long i = 0; i < num_elem; ++i)
      for (int j = 0; j < node_per_elem; ++j)
        conn_array[i*node_per_elem+type.node_order[j]] = connectivity[i*node_per_elem+j];
  }
  else
  {
    memcpy( conn_array, &connectivity[0], connectivity.size() * sizeof(MBEntityHandle) );
  }

    // Store element IDs
  result = mdbImpl->tag_set_data( globalId, elements, &elem_ids[0] );
  if (MB_SUCCESS != result)
    return result;
  
    // Add elements to material sets
  result = create_sets( type.mbtype, elements, matl_ids, 0 );
  if (MB_SUCCESS != result)
    return result;
    // Add elements to geometric sets
  result = create_sets( type.mbtype, elements, geom_ids, 1 );
  if (MB_SUCCESS != result)
    return result;
    // Add elements to parallel partitions
  result = create_sets( type.mbtype, elements, prtn_ids, 2 );
  if (MB_SUCCESS != result)
    return result;
  
  return MB_SUCCESS;
}

//! Add elements to sets as dictated by grouping ID in file.
MBErrorCode ReadGmsh::create_sets( MBEntityType type,
                                   const MBRange& elements,
                                   const std::vector<int>& set_ids,
                                   int set_type )
{ 
  MBErrorCode result;
  
    // Get a unque list of set IDs
  std::set<int> ids;
  for (std::vector<int>::const_iterator i = set_ids.begin(); i != set_ids.end(); ++i)
    ids.insert( *i );
  
    // No Sets?
  if (ids.empty() || (ids.size() == 1 && *ids.begin() == 0))
    return MB_SUCCESS; // no sets (all ids are zero)
  

    // Get/create tag handles
  int num_tags;
  MBTag tag_handles[2];
  int tag_val;
  const void* tag_values[2] = { &tag_val, NULL };
  
  switch( set_type ) 
  {
    default:
      return MB_FAILURE;
    case 0:
    case 2:
    {
      const char* name = set_type ? PARALLEL_PARTITION_TAG_NAME : MATERIAL_SET_TAG_NAME;
      result = mdbImpl->tag_get_handle( name, tag_handles[0] );
      if (result == MB_TAG_NOT_FOUND)
        result = mdbImpl->tag_create( name, 
                                      sizeof(int),
                                      MB_TAG_SPARSE,
                                      MB_TYPE_INTEGER,
                                      tag_handles[0], 0 );
      if (MB_SUCCESS != result)
        return result;
      num_tags = 1;
      break;
    }
    case 1:
    {
      const char* name = GEOM_DIMENSION_TAG_NAME;
      result = mdbImpl->tag_get_handle( name, tag_handles[1] );
      if (result == MB_TAG_NOT_FOUND)
        result = mdbImpl->tag_create( name, 
                                      sizeof(int),
                                      MB_TAG_SPARSE,
                                      MB_TYPE_INTEGER,
                                      tag_handles[1], 0 );
      if (MB_SUCCESS != result)
        return result;
      tag_values[1] = NULL;
      tag_handles[0]= globalId;
      num_tags = 2;
      break;
    }
  } // switch
  
    // For each unique set ID...
  for (std::set<int>::iterator i = ids.begin(); i != ids.end(); ++i)
  {
      // Skip "null" set ID
    if (*i == 0)
      continue;
    
      // Get all entities with the current set ID
    MBRange entities, sets;
    std::vector<int>::const_iterator j = set_ids.begin();
    for (MBRange::iterator k = elements.begin(); k != elements.end(); ++j, ++k)
      if (*i == *j)
        entities.insert( *k );
  
      // Get set by ID
    tag_val = *i;
    result = mdbImpl->get_entities_by_type_and_tag( 0, MBENTITYSET,
                                                    tag_handles, tag_values, num_tags,
                                                    sets );
    if (MB_SUCCESS != result && MB_ENTITY_NOT_FOUND != result) 
      return result;
    
      // Don't use existing geometry sets (from some other file)
    if (1 == set_type) // geometry
      sets = sets.intersect( geomSets );
    
      // Get set handle
    MBEntityHandle set;
      // If no sets with ID, create one
    if (sets.empty())
    {
      result = mdbImpl->create_meshset( MESHSET_SET, set );
      if (MB_SUCCESS != result)
        return result;
         
      result = mdbImpl->add_entities( mCurrentMeshHandle, &set, 1 );
      if (MB_SUCCESS != result)
        return result;
     
      result = mdbImpl->tag_set_data( tag_handles[0], &set, 1, &*i );
      if (MB_SUCCESS != result)
        return result;
      
      if (1 == set_type) // geometry
      {
        int dim = MBCN::Dimension( type );
        result = mdbImpl->tag_set_data( tag_handles[1], &set, 1, &dim );
        if (MB_SUCCESS != result)
          return result;
        geomSets.insert( set );
      }
    }
    else
    {
      set = *sets.begin();
      if (1 == set_type) // geometry
      {
        int dim = MBCN::Dimension( type );
          // Get dimension of set
        int dim2;
        result = mdbImpl->tag_get_data( tag_handles[1], &set, 1, &dim2 );
        if (MB_SUCCESS != result) 
          return result;
          // If we're putting geometry of a higher dimension into the
          // set, increase the dimension of the set.
        if (dim > dim2) {
          result = mdbImpl->tag_set_data( tag_handles[1], &set, 1, &dim );
          if (MB_SUCCESS != result)
            return result;
        }
      }
    }
    
      // Put the mesh entities into the set
    result = mdbImpl->add_entities( set, entities );
    if (MB_SUCCESS != result)
      return result;
  } // for (ids)
  
  return MB_SUCCESS;
}


//! NOT IMPLEMENTED
//! Reconstruct parent-child relations for geometry sets from
//! mesh connectivity.  
MBErrorCode ReadGmsh::create_geometric_topology()
{
  if (geomSets.empty())
    return MB_SUCCESS;
  
    // not implemented yet
  geomSets.clear();
  return MB_SUCCESS;
}
