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
 * \class ReadVtk
 * \brief VTK reader from Mesquite
 * \author Jason Kraftcheck
 */

#include <assert.h>

#include "ReadVtk.hpp"
#include "MBRange.hpp"
#include "MBInternals.hpp"
#include "MBInterface.hpp"
#include "MBReadUtilIface.hpp"
#include "FileTokenizer.hpp"
#include "VtkUtil.hpp"

MBReaderIface* ReadVtk::factory( MBInterface* iface )
  { return new ReadVtk( iface ); }

ReadVtk::ReadVtk(MBInterface* impl)
    : mdbImpl(impl)
{
  impl->query_interface("MBReadUtilIface", reinterpret_cast<void**>(&readMeshIface));
}

const char* const vtk_type_names[] = { "bit",
                                       "char",
                                       "unsigned_char",
                                       "short",
                                       "unsigned_short",
                                       "int",
                                       "unsigned_int",
                                       "long",
                                       "unsigned_long",
                                       "float",
                                       "double",
                                       0 };


MBErrorCode ReadVtk::load_file(const char *filename,
                               MBEntityHandle& file_set,
                               const FileOptions&,
                               const int*, const int) 
{
  MBErrorCode result;
  file_set = 0;

  int major, minor;
  char vendor_string[257];
  std::vector<MBRange> element_list;
  MBRange vertices;
  
  FILE* file = fopen( filename, "r" );
  if (!file)
  {
    return MB_FILE_DOES_NOT_EXIST;
  }
  
    // Read file header
    
  if (!fgets( vendor_string, sizeof(vendor_string), file ))
  {
    fclose( file );
    return MB_FAILURE;
  }
  
  if (!strchr( vendor_string, '\n' ) ||
      2 != sscanf( vendor_string, "# vtk DataFile Version %d.%d", &major, &minor ))
  {
    fclose( file );
    return MB_FAILURE;
  }
  
  if (!fgets( vendor_string, sizeof(vendor_string), file )) 
  {
    fclose( file );
    return MB_FAILURE;
  }
  
    // VTK spec says this should not exceed 256 chars.
  if (!strchr( vendor_string, '\n' ))
  {
    readMeshIface->report_error( "Vendor string (line 2) exceeds 256 characters." );
    fclose( file );
    return MB_FAILURE;
  }
  
  
    // Check file type
  
  FileTokenizer tokens( file, readMeshIface );
  const char* const file_type_names[] = { "ASCII", "BINARY", 0 };
  int filetype = tokens.match_token( file_type_names );
  switch (filetype) {
    case 2:  // BINARY
      readMeshIface->report_error( "Cannot read BINARY VTK files" );
    default: // ERROR 
      return MB_FAILURE;
    case 1:  // ASCII
      break;
  }

    // make a meshset for this mesh
  result = mdbImpl->create_meshset(MESHSET_SET, mCurrentMeshHandle);
  if (MB_SUCCESS != result) return result;

    // Read the mesh
  if (!tokens.match_token( "DATASET" ))
    return MB_FAILURE;
  result = vtk_read_dataset( tokens, vertices, element_list );
  if (MB_SUCCESS != result) 
    return result;
  
    // Count the number of elements read
  long elem_count = 0;
  for (std::vector<MBRange>::iterator it = element_list.begin(); it != element_list.end(); ++it )
    elem_count += it->size();
  
    // Read attribute data until end of file.
  const char* const block_type_names[] = { "POINT_DATA", "CELL_DATA", 0 };
  std::vector<MBRange> vertex_list(1);
  vertex_list[0] = vertices;
  int blocktype = 0;
  while (!tokens.eof())
  {
      // get POINT_DATA or CELL_DATA
    int new_block_type = tokens.match_token( block_type_names, false );
    if (tokens.eof())
      break;

    if (!new_block_type)
    {
        // If next token was neither POINT_DATA nor CELL_DATA,
        // then there's another attribute under the current one.
      if (blocktype)
        tokens.unget_token();
      else
        break;
    }
    else
    {
      blocktype = new_block_type;
      long count;
      if (!tokens.get_long_ints( 1, &count ))
        return MB_FAILURE;
      
      if (blocktype == 1 && (unsigned long)count != vertices.size())
      {
        readMeshIface->report_error(
              "Count inconsistent with number of vertices at line %d.", 
              tokens.line_number());
        return MB_FAILURE;
      }
      else if (blocktype == 2 && count != elem_count)
      {
        readMeshIface->report_error(
              "Count inconsistent with number of elements at line %d.", 
              tokens.line_number());
        return MB_FAILURE;
      }
    }
      
   
    if (blocktype == 1)
      result = vtk_read_attrib_data( tokens, vertex_list );
    else
      result = vtk_read_attrib_data( tokens, element_list );
    
    if (MB_SUCCESS != result)
      return result;
  }
  
  file_set = mCurrentMeshHandle;
  return MB_SUCCESS;
}

MBErrorCode ReadVtk::allocate_vertices( long num_verts,
                                        MBEntityHandle& start_handle_out,
                                        double*& x_coord_array_out,
                                        double*& y_coord_array_out,
                                        double*& z_coord_array_out )
{
  MBErrorCode result;
  
    // Create vertices
  std::vector<double*> arrays;
  start_handle_out = 0;
  result = readMeshIface->get_node_arrays( 3, num_verts, MB_START_ID, 0,
                                           start_handle_out, arrays );
  if (MB_SUCCESS != result)
    return result;
    
  x_coord_array_out = arrays[0];
  y_coord_array_out = arrays[1];
  z_coord_array_out = arrays[2];
  
    // Add vertices to meshset
  MBRange vert_range(start_handle_out, start_handle_out+num_verts-1);
  return mdbImpl->add_entities(mCurrentMeshHandle, vert_range);
}

MBErrorCode ReadVtk::read_vertices( FileTokenizer& tokens,
                                    long num_verts, 
                                    MBEntityHandle& start_handle_out )
{
  MBErrorCode result;
  double *x, *y, *z;
  
  result = allocate_vertices( num_verts, start_handle_out, x, y, z );
  if (MB_SUCCESS != result)
    return result;

    // Read vertex coordinates
  for (long vtx = 0; vtx < num_verts; ++vtx)
  {
    if (!tokens.get_doubles( 1, x++ ) ||
        !tokens.get_doubles( 1, y++ ) ||
        !tokens.get_doubles( 1, z++ ))
      return MB_FAILURE;
  }
  
  return MB_SUCCESS;
}

MBErrorCode ReadVtk::allocate_elements( long num_elements,
                                        int vert_per_element,
                                        MBEntityType type,
                                        MBEntityHandle& start_handle_out,
                                        MBEntityHandle*& conn_array_out,
                                        std::vector<MBRange>& append_to_this )
{
  MBErrorCode result;
  
  start_handle_out = 0;
  result = readMeshIface->get_element_array( num_elements,
                                             vert_per_element,
                                             type,
                                             MB_START_ID, 0,
                                             start_handle_out,
                                             conn_array_out );
  if (MB_SUCCESS != result)
    return result;
  
  MBRange range(start_handle_out, start_handle_out+num_elements-1);
  append_to_this.push_back( range );
  return mdbImpl->add_entities(mCurrentMeshHandle, range);
}

MBErrorCode ReadVtk::allocate_poly_elems( long num_elements,
                                          int connectivity_length,
                                          MBEntityType type,
                                          MBEntityHandle& start_handle_out,
                                          MBEntityHandle*& conn_array_out,
                                          int*& index_array_out,
                                          std::vector<MBRange>& append_to_this )
{
  MBErrorCode result;
  
  start_handle_out = 0;
  result = readMeshIface->get_poly_element_array( num_elements,
                                                  connectivity_length,
                                                  type,
                                                  0, 0,
                                                  start_handle_out,
                                                  index_array_out,
                                                  conn_array_out );
  if (MB_SUCCESS != result)
    return result;
  
  MBRange range(start_handle_out, start_handle_out+num_elements-1);
  append_to_this.push_back( range );
  return mdbImpl->add_entities(mCurrentMeshHandle, range);
}

MBErrorCode ReadVtk::vtk_read_dataset( FileTokenizer& tokens,
                                       MBRange& vertex_list,
                                       std::vector<MBRange>& element_list )
{
  const char* const data_type_names[] = { "STRUCTURED_POINTS",
                                          "STRUCTURED_GRID",
                                          "UNSTRUCTURED_GRID",
                                          "POLYDATA",
                                          "RECTILINEAR_GRID",
                                          "FIELD",
                                          0 };
  int datatype = tokens.match_token( data_type_names );
  switch( datatype )
  {
    case 1:  return vtk_read_structured_points( tokens, vertex_list, element_list );
    case 2:  return vtk_read_structured_grid  ( tokens, vertex_list, element_list );
    case 3:  return vtk_read_unstructured_grid( tokens, vertex_list, element_list );
    case 4:  return vtk_read_polydata         ( tokens, vertex_list, element_list );
    case 5:  return vtk_read_rectilinear_grid ( tokens, vertex_list, element_list );
    case 6:  return vtk_read_field            ( tokens );
    default: return MB_FAILURE;
  }
  
}


MBErrorCode ReadVtk::vtk_read_structured_points( FileTokenizer& tokens, 
                                                 MBRange& vertex_list,
                                                 std::vector<MBRange>& elem_list )
{
  long i, j, k;
  long dims[3];
  double origin[3], space[3];
  MBErrorCode result;
 
  if (!tokens.match_token( "DIMENSIONS" ) || 
      !tokens.get_long_ints( 3, dims ) ||
      !tokens.get_newline( ))
    return MB_FAILURE; 
  
  if (dims[0] < 1 || dims[1] < 1 || dims[2] < 1)
  {
    readMeshIface->report_error( "Invalid dimension at line %d", 
                                 tokens.line_number());
    return MB_FAILURE;
  }
  
  if (!tokens.match_token( "ORIGIN" )  ||
      !tokens.get_doubles( 3, origin ) ||
      !tokens.get_newline( ))
    return MB_FAILURE;           
  
  const char* const spacing_names[] = { "SPACING", "ASPECT_RATIO", 0 };
  if (!tokens.match_token( spacing_names ) ||
      !tokens.get_doubles( 3, space )      ||
      !tokens.get_newline())
    return MB_FAILURE;
  
    // Create vertices
  double *x, *y, *z;
  MBEntityHandle start_handle = 0;
  long num_verts = dims[0]*dims[1]*dims[2];
  result = allocate_vertices( num_verts, start_handle, x, y, z );
  if (MB_SUCCESS != result)
    return result;
  vertex_list.insert( start_handle, start_handle + num_verts - 1 );
  
  for (k = 0; k < dims[2]; ++k)
    for (j = 0; j < dims[1]; ++j)
      for (i = 0; i < dims[0]; ++i)
      {
        *x = origin[0] + i*space[0]; ++x;
        *y = origin[1] + j*space[1]; ++y;
        *z = origin[2] + k*space[2]; ++z;
      }

  return vtk_create_structured_elems( dims, start_handle, elem_list );
}

MBErrorCode ReadVtk::vtk_read_structured_grid( FileTokenizer& tokens, 
                                               MBRange& vertex_list,
                                               std::vector<MBRange>& elem_list )
{
  long num_verts, dims[3];
  MBErrorCode result;
 
  if (!tokens.match_token( "DIMENSIONS" ) ||
      !tokens.get_long_ints( 3, dims )    ||
      !tokens.get_newline( ) )
    return MB_FAILURE;            
  
  if (dims[0] < 1 || dims[1] < 1 || dims[2] < 1)
  {
    readMeshIface->report_error( "Invalid dimension at line %d", 
                                 tokens.line_number());
    return MB_FAILURE;
  }
  
  if (!tokens.match_token( "POINTS" )        ||
      !tokens.get_long_ints( 1, &num_verts ) ||
      !tokens.match_token( vtk_type_names )  ||
      !tokens.get_newline())
    return MB_FAILURE;                
  
  if (num_verts != (dims[0] * dims[1] * dims[2]))
  {
    readMeshIface->report_error(                     
      "Point count not consistent with dimensions at line %d", 
      tokens.line_number() );
    return MB_FAILURE;
  }
  
    // Create and read vertices
  MBEntityHandle start_handle = 0;
  result = read_vertices( tokens, num_verts, start_handle );
  if (MB_SUCCESS != result)
    return result;
  vertex_list.insert( start_handle, start_handle + num_verts - 1 );
 
  return vtk_create_structured_elems( dims, start_handle, elem_list );
}

MBErrorCode ReadVtk::vtk_read_rectilinear_grid( FileTokenizer& tokens, 
                                                MBRange& vertex_list,
                                                std::vector<MBRange>& elem_list )
{
  int i, j, k;
  long dims[3];
  const char* labels[] = { "X_COORDINATES", "Y_COORDINATES", "Z_COORDINATES" };
  std::vector<double> coords[3];
  MBErrorCode result;
  
  if (!tokens.match_token( "DIMENSIONS" )||
      !tokens.get_long_ints( 3, dims )   ||
      !tokens.get_newline( ))
    return MB_FAILURE;
  
  if (dims[0] < 1 || dims[1] < 1 || dims[2] < 1)
  {
    readMeshIface->report_error( "Invalid dimension at line %d", 
                                 tokens.line_number());
    return MB_FAILURE;
  }

  for (i = 0; i < 3; i++)
  {
    long count;
    if (!tokens.match_token( labels[i] )      ||
        !tokens.get_long_ints( 1, &count )    ||
        !tokens.match_token( vtk_type_names ))
      return MB_FAILURE;
    
    if (count != dims[i])
    {
      readMeshIface->report_error(                     
        "Coordinate count inconsistent with dimensions at line %d", 
        tokens.line_number() );
      return MB_FAILURE;
    }
    
    coords[i].resize(count);
    if (!tokens.get_doubles( count, &coords[i][0] ))
      return MB_FAILURE;
  }
  
    // Create vertices
  double *x, *y, *z;
  MBEntityHandle start_handle = 0;
  long num_verts = dims[0]*dims[1]*dims[2];
  result = allocate_vertices( num_verts, start_handle, x, y, z );
  if (MB_SUCCESS != result)
    return result;
  vertex_list.insert( start_handle, start_handle + num_verts - 1 );
  
    // Calculate vertex coordinates
  for (k = 0; k < dims[2]; ++k)
    for (j = 0; j < dims[1]; ++j)
      for (i = 0; i < dims[0]; ++i)
      {
        *x = coords[0][i]; ++x;
        *y = coords[1][j]; ++y;
        *z = coords[2][k]; ++z;
      }
  
  
  return vtk_create_structured_elems( dims, start_handle, elem_list );
}

MBErrorCode ReadVtk::vtk_read_polydata( FileTokenizer& tokens, 
                                        MBRange& vertex_list,
                                        std::vector<MBRange>& elem_list )
{
  MBErrorCode result;
  long num_verts;
  std::vector<int> connectivity;
  const char* const poly_data_names[] = { "VERTICES",
                                          "LINES",
                                          "POLYGONS",
                                          "TRIANGLE_STRIPS", 
                                          0 };
  
  if (!tokens.match_token( "POINTS" )        ||
      !tokens.get_long_ints( 1, &num_verts ) ||
      !tokens.match_token( vtk_type_names )  ||
      !tokens.get_newline( ))
    return MB_FAILURE;
  
  if (num_verts < 1)
  {
    readMeshIface->report_error( "Invalid point count at line %d", 
                                 tokens.line_number());
    return MB_FAILURE;
  }
  
    // Create vertices and read coordinates
  MBEntityHandle start_handle = 0;
  result = read_vertices( tokens, num_verts, start_handle );
  if (MB_SUCCESS != result)
    return result;
  vertex_list.insert( start_handle, start_handle + num_verts - 1 );

  int poly_type = tokens.match_token( poly_data_names );
  switch (poly_type)
  {
    case 0:
      result = MB_FAILURE;
    case 1:
      readMeshIface->report_error( 
                       "Vertex element type at line %d",
                       tokens.line_number() );
      result = MB_FAILURE;
    case 2:
      readMeshIface->report_error( 
               "Unsupported type: polylines at line %d",
               tokens.line_number() );
      result = MB_FAILURE;
    case 3:
      result = vtk_read_polygons( tokens, start_handle, elem_list );
    case 4:
      readMeshIface->report_error( 
               "Unsupported type: triangle strips at line %d",
               tokens.line_number() );
      result = MB_FAILURE;
  }
  
  return result;
}

MBErrorCode ReadVtk::vtk_read_polygons( FileTokenizer& tokens,
                                        MBEntityHandle first_vtx, 
                                        std::vector<MBRange>& elem_list )
{
  MBErrorCode result;
  long size[2];
  
  if (!tokens.get_long_ints( 2, size ) ||
      !tokens.get_newline( ))
    return MB_FAILURE;

    // Create polygons
  MBEntityHandle first_poly = 0;
  MBEntityHandle* conn_array = 0;
  int* index_array = 0;
  result = allocate_poly_elems( size[0], size[1] - size[0], MBPOLYGON,
                                first_poly, conn_array, index_array,
                                elem_list );
  if (MB_SUCCESS != result)
    return result;
                                                  
    // Read connectivity
  int prev_index = -1;
  for (int i = 0; i < size[0]; ++i)
  {
    long count;
    if (!tokens.get_long_ints( 1, &count ) ||
        !tokens.get_long_ints( count, (long*)conn_array ))
      return MB_FAILURE;
    
    *index_array = ( prev_index += count );
    for (long j = 0; j < count; ++j)
    {
      *conn_array += first_vtx;
      ++conn_array;
    }
  } 
  
  return MB_SUCCESS;
}



MBErrorCode ReadVtk::vtk_read_unstructured_grid( FileTokenizer& tokens,
                                                 MBRange& vertex_list,
                                                 std::vector<MBRange>& elem_list  )
{
  MBErrorCode result;
  long i, num_verts, num_elems[2];
  MBEntityHandle tmp_conn_list[27];

  if (!tokens.match_token( "POINTS" )        ||
      !tokens.get_long_ints( 1, &num_verts ) ||
      !tokens.match_token( vtk_type_names)   ||
      !tokens.get_newline( ))
    return MB_FAILURE;
  
  if (num_verts < 1)
  {
    readMeshIface->report_error( "Invalid point count at line %d", tokens.line_number());
    return MB_FAILURE;
  }
  
    // Create vertices and read coordinates
  MBEntityHandle first_vertex = 0;
  result = read_vertices( tokens, num_verts, first_vertex );
  if (MB_SUCCESS != result)
    return result;
  vertex_list.insert( first_vertex, first_vertex + num_verts - 1 );
  
  if (!tokens.match_token( "CELLS" )        ||
      !tokens.get_long_ints( 2, num_elems ) ||
      !tokens.get_newline( ))
    return MB_FAILURE;
  
    // Read element connectivity for all elements
  std::vector<long> connectivity( num_elems[1] );
  if (!tokens.get_long_ints( num_elems[1], &connectivity[0] ))
    return MB_FAILURE;
  
  if (!tokens.match_token( "CELL_TYPES" )      ||
      !tokens.get_long_ints( 1, &num_elems[1]) ||
      !tokens.get_newline( ))
    return MB_FAILURE;

    // Read element types
  std::vector<long> types( num_elems[0] );
  if (!tokens.get_long_ints( num_elems[0], &types[0] ))
    return MB_FAILURE;
  
    // Create elements in blocks of the same type
    // It is important to preserve the order in 
    // which the elements were read for later reading
    // attribute data.
  long id = 0;
  std::vector<long>::iterator conn_iter = connectivity.begin();
  while (id < num_elems[0])
  {
    unsigned vtk_type = types[id];
    if (vtk_type >= VtkUtil::numVtkElemType)
      return MB_FAILURE;
      
    MBEntityType type = VtkUtil::vtkElemTypes[vtk_type].mb_type;
    int num_vtx = VtkUtil::vtkElemTypes[vtk_type].num_nodes;
    
    if (type == MBMAXTYPE) {
      readMeshIface->report_error( "Unsupported VTK element type: %s (%d)\n",
                                   VtkUtil::vtkElemTypes[vtk_type].name, vtk_type );
      return MB_FAILURE;
    }
    
      // Find any subsequent elements of the same type
    long end_id = id + 1;
    while ( end_id < num_elems[0] && (unsigned)types[end_id] == vtk_type)
      ++end_id;
    
      // Allocate element block
    long num_elem = end_id - id;
    MBEntityHandle start_handle = 0;
    MBEntityHandle* conn_array;
    int* index_array = 0;
    
    if (type == MBPOLYGON)
    {
        // Calculate total length of connectivity list
      std::vector<long>::iterator conn_iter2 = conn_iter;
      long conn_len = 0;
      for (i = 0; i < num_elem; ++i)
      {
        conn_len += *conn_iter2;
        conn_iter2 += *conn_iter2 + 1;
      }
      
        // Allocate elements
      result = allocate_poly_elems( num_elem, conn_len, type, start_handle,
                                    conn_array, index_array, elem_list );
    }
    else
    {
      result = allocate_elements( num_elem, num_vtx, type, start_handle,
                                  conn_array, elem_list );
    }
    
    if (MB_SUCCESS != result)
      return result;

      // Store element connectivity
    int index_array_prev = -1;
    for ( ; id < end_id; ++id)
    {
      if (conn_iter == connectivity.end()) 
      {
        readMeshIface->report_error( 
           "Connectivity data truncated at cell %ld\n", id );
        return MB_FAILURE;
      }

        // For poly elements, store end index
      if (index_array)
      {
        index_array_prev += num_vtx = *conn_iter;
        *index_array = index_array_prev;
        ++index_array;
      }
        // For non-poly elements, make sure connectivity length is correct.
      else if (*conn_iter != num_vtx)
      {
        readMeshIface->report_error(
          "Cell %ld is of type '%s' but has %u vertices.\n",
          id, VtkUtil::vtkElemTypes[vtk_type].name, num_vtx );
        return MB_FAILURE;
      }
      ++conn_iter;

      for (i = 0; i < num_vtx; ++i, ++conn_iter)
      {
        if (conn_iter == connectivity.end()) 
        {
          readMeshIface->report_error( 
             "Connectivity data truncated at cell %ld\n", id );
          return MB_FAILURE;
        }

        conn_array[i] = *conn_iter + first_vertex;
      }

      const unsigned* order = VtkUtil::vtkElemTypes[vtk_type].node_order;
      if ( order )
      {
        assert( num_vtx * sizeof(MBEntityHandle) <= sizeof(tmp_conn_list) );
        memcpy( tmp_conn_list, conn_array, num_vtx * sizeof(MBEntityHandle) );
        for (int j = 0; j < num_vtx; ++j)
          conn_array[order[j]] = tmp_conn_list[j];
      }       

      conn_array += num_vtx;
    }
  }      

  return MB_SUCCESS;
}

MBErrorCode ReadVtk::vtk_create_structured_elems( const long* dims, 
                                            MBEntityHandle first_vtx,
                                            std::vector<MBRange>& elem_list )
{
  MBErrorCode result;
  int non_zero[3] = {0,0,0};  // True if dim > 0 for x, y, z respectively
  long elem_dim = 0;          // Element dimension (2->quad, 3->hex)
  long num_elems = 1;         // Total number of elements
  long vert_per_elem;         // Element connectivity length
  long edims[3] = { 1, 1, 1 };// Number of elements in each grid direction
  
    // Populate above data
  for (int d = 0; d < 3; d++) 
    if (dims[d] > 1)
    {
      non_zero[elem_dim++] = d;
      edims[d] = dims[d] - 1;
      num_elems *= edims[d];
    }
  vert_per_elem = 1 << elem_dim;
  
    // Get element type from element dimension
  MBEntityType type;
  switch( elem_dim )
  {
    case 1: type = MBEDGE; break;
    case 2: type = MBQUAD; break;
    case 3: type = MBHEX;  break;
    default:
      readMeshIface->report_error( "Invalid dimension for structured elements: %ld\n",
                                   elem_dim );
      return MB_FAILURE;
  }

    // Allocate storage for elements
  MBEntityHandle start_handle = 0;
  MBEntityHandle* conn_array;
  result = allocate_elements( num_elems, vert_per_elem, type, start_handle,
                              conn_array, elem_list );
  if (MB_SUCCESS != result)
    return MB_FAILURE;

    // Offsets of element vertices in grid relative to corner closest to origin 
  long k = dims[0]*dims[1];
  const long corners[8] = { 0, 1, 1+dims[0], dims[0], k, k+1, k+1+dims[0], k+dims[0] };
                             
    // Populate element list
  for (long z = 0; z < edims[2]; ++z)
    for (long y = 0; y < edims[1]; ++y)
      for (long x = 0; x < edims[0]; ++x)
      {
        const long index = x + y*dims[0] + z*(dims[0]*dims[1]);
        for (long j = 0; j < vert_per_elem; ++j, ++conn_array)
          *conn_array = index + corners[j] + first_vtx;
      }
  
  return MB_SUCCESS;
}

MBErrorCode ReadVtk::vtk_read_field( FileTokenizer& tokens )
{
    // This is not supported yet.
    // Parse the data but throw it away because
    // Mesquite has no internal representation for it.
  
    // Could save this in tags, but the only useful thing that
    // could be done with the data is to write it back out
    // with the modified mesh.  As there's no way to save the
    // type of a tag in Mesquite, it cannot be written back
    // out correctly either.
    // FIXME: Don't know what to do with this data.
    // For now, read it and throw it out.

  long num_arrays;
  if (!tokens.get_long_ints( 1, &num_arrays ))
    return MB_FAILURE;
  
  for (long i = 0; i < num_arrays; ++i)
  {
    /*const char* name =*/ tokens.get_string( );
    
    long dims[2];
    if (!tokens.get_long_ints( 2, dims ) ||
        !tokens.match_token( vtk_type_names ))
      return MB_FAILURE;
    
    long num_vals = dims[0] * dims[1];
    
    for (long j = 0; j < num_vals; j++)
    {
      double junk;
      if (!tokens.get_doubles( 1, &junk ))
        return MB_FAILURE;
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode ReadVtk::vtk_read_attrib_data( FileTokenizer& tokens, 
                                           std::vector<MBRange>& entities )
{
  const char* const type_names[] = { "SCALARS",
                                     "COLOR_SCALARS",
                                     "VECTORS",
                                     "NORMALS",
                                     "TEXTURE_COORDINATES",
                                     "TENSORS",
                                     "FIELD",
                                     0 };

  int type = tokens.match_token( type_names );
  const char* tmp_name = tokens.get_string( );
  if (!type || !tmp_name)
    return MB_FAILURE;
  
  std::string name_alloc( tmp_name );
  const char* name = name_alloc.c_str();
  switch( type )
  {
    case 1: return vtk_read_scalar_attrib ( tokens, entities, name ); 
    case 2: return vtk_read_color_attrib  ( tokens, entities, name ); 
    case 3: return vtk_read_vector_attrib ( tokens, entities, name );
    case 4: return vtk_read_vector_attrib ( tokens, entities, name ); 
    case 5: return vtk_read_texture_attrib( tokens, entities, name ); 
    case 6: return vtk_read_tensor_attrib ( tokens, entities, name ); 
    case 7: // Can't handle field data yet.
      readMeshIface->report_error( "Error at line %d: field data not implemented.\n",
                                   tokens.line_number() );
    default:
      return MB_FAILURE;
  }

  return MB_SUCCESS;
}

MBErrorCode ReadVtk::vtk_read_tag_data( FileTokenizer& tokens, 
                                        int type, 
                                        size_t per_elem, 
                                        std::vector<MBRange>& entities,
                                        const char* name )
{
  MBErrorCode result;
  
  MBDataType mb_type;
  size_t size;
  if (type == 1)
  {
    mb_type = MB_TYPE_BIT;
    size = sizeof(bool);
  }
  else if (type >= 2 && type <= 9)
  {
    mb_type = MB_TYPE_INTEGER;
    size = sizeof(int);
  }
  else if (type == 10 || type == 11)
  {
    mb_type = MB_TYPE_DOUBLE;
    size = sizeof(double);
  }
  else
  {
    return MB_FAILURE;
  }
  
    // get/create tag
  MBTag handle;
  result = mdbImpl->tag_get_handle ( name, handle );
  if (result == MB_TAG_NOT_FOUND)
  {
    if (mb_type == MB_TYPE_BIT)
      result = mdbImpl->tag_create( name, per_elem, MB_TAG_BIT, mb_type, handle, 0 );
    else
      result = mdbImpl->tag_create( name, size*per_elem, MB_TAG_DENSE, mb_type, handle, 0 );
    if (MB_SUCCESS != result)
      return result;
  }
  else if (result == MB_SUCCESS)
  {
    int existing_size;
    MBDataType existing_type;
    
    result = mdbImpl->tag_get_size( handle, existing_size );
    if (MB_SUCCESS != result)
      return result;
    result = mdbImpl->tag_get_data_type( handle, existing_type );
    if (MB_SUCCESS != result)
      return result;
    
    if ( (mb_type == MB_TYPE_BIT && (size_t)existing_size != per_elem) ||
         (mb_type != MB_TYPE_BIT && (size_t)existing_size != size*per_elem) || 
         (existing_type != mb_type) )
    {
      readMeshIface->report_error( 
        "Tag name conflict for attribute \"%s\" at line %d\n",
        name, tokens.line_number() );
      return MB_FAILURE;
    }
  }
  else
  {
    return result;
  }
 
    // Count number of entities
  long count = 0;
  std::vector<MBRange>::iterator iter;
  for (iter = entities.begin(); iter != entities.end(); ++iter)
    count += iter->size();

  if (type == 1)
  {
    for (iter = entities.begin(); iter != entities.end(); ++iter)
    {
      bool *data = new bool[iter->size() * per_elem];
      if (!tokens.get_booleans( per_elem*iter->size(), data ))
      {
        delete [] data;
        return MB_FAILURE;
      }
      
      bool* data_iter = data;
      MBRange::iterator ent_iter = iter->begin();
      for ( ; ent_iter != iter->end(); ++ent_iter)
      {
        unsigned char bits = 0;
        for (unsigned j = 0; j < per_elem; ++j, ++data_iter)
          bits |= (unsigned char)(*data_iter << j);
        result = mdbImpl->tag_set_data( handle, &*ent_iter, 1, &bits );
        if (MB_SUCCESS != result)
        {
          delete [] data;
          return result;
        }
      } 
      delete [] data;
    }
  }
  else if (type >= 2 && type <= 9)
  {
    std::vector<int> data;
    for (iter = entities.begin(); iter != entities.end(); ++iter)
    {
      data.resize( iter->size() * per_elem );
      if (!tokens.get_integers( iter->size() * per_elem, &data[0] ))
        return MB_FAILURE;
      result = mdbImpl->tag_set_data( handle, *iter, &data[0] );
      if (MB_SUCCESS != result)
        return result;
    }
  }
  else if (type == 10 || type == 11)
  {
    std::vector<double> data;
    for (iter = entities.begin(); iter != entities.end(); ++iter)
    {
      data.resize( iter->size() * per_elem );
      if (!tokens.get_doubles( iter->size() * per_elem, &data[0] ))
        return MB_FAILURE;
      result = mdbImpl->tag_set_data( handle, *iter, &data[0] );
      if (MB_SUCCESS != result)
        return result;
    }
  }
  else
  {
    return MB_FAILURE;
  }
  
  return MB_SUCCESS;
}
  
      
MBErrorCode ReadVtk::vtk_read_scalar_attrib( FileTokenizer& tokens, 
                                             std::vector<MBRange>& entities,
                                             const char* name)
{
  int type = tokens.match_token( vtk_type_names );
  if (!type)
    return MB_FAILURE;
    
  long size;
  const char* tok = tokens.get_string( );
  if (!tok)
    return MB_FAILURE;
    
  const char* end = 0;
  size = strtol( tok, (char**)&end, 0 );
  if (*end)
  {
    size = 1;
    tokens.unget_token();
  }
  
    // VTK spec says cannot be greater than 4--do we care?
  if (size < 1 || size > 4)
  {
    readMeshIface->report_error(
                    "Scalar count out of range [1,4]" 
                    " at line %d", tokens.line_number());
    return MB_FAILURE;
  }
  
  if (!tokens.match_token("LOOKUP_TABLE") ||
      !tokens.match_token("default"))
    return MB_FAILURE;
  
  return vtk_read_tag_data( tokens, type, size, entities, name );
}


MBErrorCode ReadVtk::vtk_read_color_attrib( FileTokenizer& tokens, 
                                            std::vector<MBRange>& entities,
                                            const char* name)
{
  long size;
  if (!tokens.get_long_ints( 1, &size ) || size < 1)
    return MB_FAILURE;

  return vtk_read_tag_data( tokens, 10, size, entities, name );
}

MBErrorCode ReadVtk::vtk_read_vector_attrib( FileTokenizer& tokens, 
                                             std::vector<MBRange>& entities,
                                             const char* name)
{
  int type = tokens.match_token( vtk_type_names );
  if (!type)
    return MB_FAILURE;
    
  return vtk_read_tag_data( tokens, type, 3, entities, name );
}

MBErrorCode ReadVtk::vtk_read_texture_attrib( FileTokenizer& tokens,
                                              std::vector<MBRange>& entities,
                                              const char* name)
{
  int type, dim;
  if (!tokens.get_integers( 1, &dim ) ||
      !(type = tokens.match_token( vtk_type_names )))
    return MB_FAILURE;
    
  if (dim < 1 || dim > 3)
  {
    readMeshIface->report_error( "Invalid dimension (%d) at line %d.",
                     dim, tokens.line_number() );
    return MB_FAILURE;
  }
  
  return vtk_read_tag_data( tokens, type, dim, entities, name );
}

MBErrorCode ReadVtk::vtk_read_tensor_attrib( FileTokenizer& tokens,
                                             std::vector<MBRange>& entities,
                                             const char* name ) 
{
  int type = tokens.match_token( vtk_type_names );
  if (!type)
    return MB_FAILURE;
    
  return vtk_read_tag_data( tokens, type, 9, entities, name );
}  

