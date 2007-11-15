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
 * \class ReadSTL
 * \brief ASCII and Binary Stereo Lithography File readers.
 * \author Jason Kraftcheck
 */

#include "ReadSTL.hpp"
#include "FileTokenizer.hpp" // for FileTokenizer
#include "MBInternals.hpp"
#include "MBInterface.hpp"
#include "MBReadUtilIface.hpp"
#include "MBRange.hpp"
#include "MBEntityHandle.h"

#ifdef MOAB_HAVE_INTTYPES_H
#include <inttypes.h>
#endif
#ifdef MOAB_HAVE_STDDEF_H
#include <stddef.h>
#endif
#ifdef MOAB_HAVE_STDINT_H
#include <stdint.h>
#endif

#ifdef _MSC_VER /* windows */
#  include <BaseTsd.h>
typedef ULONG32 uint32_t;
#endif

#include <errno.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <map>


inline static uint32_t byte_swap( uint32_t value )
{
  return ((value & 0xFF000000) >> 24) |
         ((value & 0x00FF0000) >>  8) |
         ((value & 0x0000FF00) <<  8) |
         ((value & 0X000000FF) << 24);
}


inline static float byte_swap( float value )
{
  uint32_t bytes = byte_swap( *(uint32_t*)&value );
  return *(float*)&bytes;
}


inline static bool is_platform_little_endian()
{
  static const unsigned int one = 1;
  static const bool little = !*((char*)&one);
  return little;
}



ReadSTL::ReadSTL(MBInterface* impl)
    : mdbImpl(impl)
{
  mdbImpl->query_interface("MBReadUtilIface", reinterpret_cast<void**>(&readMeshIface));
}

ReadSTL::~ReadSTL()
{
  if (readMeshIface)
   mdbImpl->release_interface("MBReadUtilIface", reinterpret_cast<void**>(&readMeshIface));
}

// Used to put points in an STL tree-based container
bool ReadSTL::Point::operator<( const ReadSTL::Point& other ) const
{
  return 0 > memcmp( this, &other, sizeof(ReadSTL::Point) );
}


MBErrorCode ReadSTL::load_file( const char* filename, 
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
  
  return result;
}

// Generic load function for both ASCII and binary.  Calls
// pure-virtual function implemented in subclasses to read
// the data from the file.
MBErrorCode ReadSTL::load_file_impl(const char *filename,
                                    const int*, const int) 
{
  MBErrorCode result;

  std::vector<ReadSTL::Triangle> triangles;
 
  result = this->read_triangles( filename, triangles );
  if (MB_SUCCESS != result) return result; 

    // make a meshset for this mesh
  result = mdbImpl->create_meshset(MESHSET_SET, mCurrentMeshHandle);
  if (MB_SUCCESS != result) return result;

    // Create a std::map from position->handle, and such
    // that all positions are specified, and handles are zero.
  std::map<Point,MBEntityHandle> vertex_map;
  for (std::vector<Triangle>::iterator i = triangles.begin(); i != triangles.end(); ++i)
  {
    vertex_map[i->points[0]] = 0;
    vertex_map[i->points[1]] = 0;
    vertex_map[i->points[2]] = 0;
  }
  
    // Create vertices 
  std::vector<double*> coord_arrays;
  MBEntityHandle handle = 0;
  result = readMeshIface->get_node_arrays( 3, vertex_map.size(), MB_START_ID, -1,
                                           handle, coord_arrays );
  if (MB_SUCCESS != result)
    return result;

    // Add vertices to entity set
  MBRange range(handle, handle+vertex_map.size()-1);
  result = mdbImpl->add_entities(mCurrentMeshHandle, range);
  if (MB_SUCCESS != result)
    return result;
  
    // Copy vertex coordinates into entity sequence coordinate arrays
    // and copy handle into vertex_map.
  double *x = coord_arrays[0], *y = coord_arrays[1], *z = coord_arrays[2];
  for (std::map<Point,MBEntityHandle>::iterator i = vertex_map.begin();
       i != vertex_map.end(); ++i)
  {
    i->second = handle; ++handle;
    *x = i->first.coords[0]; ++x;
    *y = i->first.coords[1]; ++y;
    *z = i->first.coords[2]; ++z;
  }
  
    // Allocate triangles
  handle = 0;
  MBEntityHandle* connectivity;
  result = readMeshIface->get_element_array( triangles.size(),
                                             3,
                                             MBTRI,
                                             MB_START_ID, -1,
                                             handle,
                                             connectivity );
  if (MB_SUCCESS != result)
    return result;

    // Add triangles to entity set
  MBRange range2(handle, handle+triangles.size()-1);
  result = mdbImpl->add_entities(mCurrentMeshHandle, range2);
  if (MB_SUCCESS != result)
    return result;
  
    // Use vertex_map to reconver triangle connectivity from
    // vertex coordinates.
  for (std::vector<Triangle>::iterator i = triangles.begin(); i != triangles.end(); ++i)
  {
    *connectivity = vertex_map[i->points[0]]; ++connectivity;
    *connectivity = vertex_map[i->points[1]]; ++connectivity;
    *connectivity = vertex_map[i->points[2]]; ++connectivity;
  }
  
  return MB_SUCCESS;
}


long ReadSTL::get_file_size( FILE* file )
{
  long curr_pos = ftell( file );
  if (fseek( file, 0, SEEK_END ))
    return -1;
  
  long length = ftell( file );
  if (fseek( file, curr_pos, SEEK_SET))
  {
    assert(0); 
    return -2;
  }
  
  return length;
}


// Read ASCII file
MBErrorCode ReadASCIISTL::read_triangles( const char* name,
                                          std::vector<ReadSTL::Triangle>& tris )
{
  FILE* file = fopen( name, "r" );
  if (!file)
  {
    readMeshIface->report_error( "%s: %s\n", name, strerror(errno) );
    return MB_FILE_DOES_NOT_EXIST;
  }
  
  char header[82];
  if (!fgets( header, sizeof(header), file ) ||  // read header line
      strlen( header ) < 6                   ||  // must be at least 6 chars
      header[strlen(header) - 1] != '\n'     ||  // cannot exceed 80 chars
      memcmp( header, "solid", 5 )           ||  // must begin with "solid"
      !isspace( header[5] ) )                    // followed by a whitespace char
  {
    readMeshIface->report_error( "%s: %s\n", name, strerror(errno) );
    fclose( file );
    return MB_FILE_WRITE_ERROR;
  }
  
    // Use tokenizer for remainder of parsing
  FileTokenizer tokens( file, readMeshIface );
  
  Triangle tri;
  float norm[3];
  
    // Read until end of file.  If we reach "endsolid", read
    // was successful.  If EOF before "endsolid", return error.
  for (;;)
  {
      // Check for either another facet or the end of the list.
    const char* const expected[] = { "facet", "endsolid", 0 };
    switch( tokens.match_token( expected ) )
    {
      case 1:  break;                        // found another facet
      case 2:  return MB_SUCCESS;            // found "endsolid" -- done
      default: return MB_FILE_WRITE_ERROR;   // found something else, or EOF
    }
    
    if (!tokens.match_token( "normal" ) ||    // expect "normal" keyword
        !tokens.get_floats( 3, norm )   ||    // followed by normal vector
        !tokens.match_token( "outer" )  ||    // followed by "outer loop"
        !tokens.match_token( "loop" )   )
      return MB_FILE_WRITE_ERROR;
    
      // for each of three triangle vertices
    for (int i = 0; i < 3; i++)
    {
      if (!tokens.match_token("vertex") ||
          !tokens.get_floats( 3, tri.points[i].coords ))
        return MB_FILE_WRITE_ERROR;
    }
    
    if (!tokens.match_token( "endloop" ) ||   // facet ends with "endloop"
        !tokens.match_token( "endfacet" ) )   // and then "endfacet"
      return MB_FILE_WRITE_ERROR;
    
    tris.push_back( tri );
  }
}


// Header block from binary STL file (84 bytes long)
struct BinaryHeader {
  char comment[80]; // 80 byte comment string (null terminated?)
  uint32_t count;   // number of triangles - 4 byte integer
}; 

// Triangle spec from file (50 bytes)
struct BinaryTri {
  float normal[3];  // Normal as 3 4-byte little-endian IEEE floats
  float coords[9];  // Vertex coords as 9 4-byte little-endian IEEE floats
  char pad[2];
};

// Read a binary STL file
MBErrorCode ReadBinarySTL::read_triangles( const char* name,
                                           std::vector<ReadSTL::Triangle>& tris )
{
  FILE* file = fopen( name, "rb" );
  if (!file)
  {
    readMeshIface->report_error( "%s: %s\n", name, strerror(errno) );
    return MB_FILE_DOES_NOT_EXIST;
  }
  
    // Read header block
  BinaryHeader header;
  if (fread( &header, 84, 1, file ) != 1)
  {
    fclose ( file );
    readMeshIface->report_error( "%s: %s\n", name, strerror(errno) );
    return MB_FILE_WRITE_ERROR;
  }
  
  bool swap_bytes = !is_platform_little_endian();  // default to little endian

    // Check for tag specifying file byte order
  MBTag bo_tag = 0;
  MBErrorCode rval = mdbImpl->tag_get_handle( "__STL_BYTE_ORDER", bo_tag );
  if (MB_SUCCESS == rval)
  {
    int value;
    rval = mdbImpl->tag_get_data( bo_tag, 0, 1, &value );
    if (MB_SUCCESS != rval) 
      return rval;
    bool is_file_little_endian = (0 == value);
    swap_bytes = (is_platform_little_endian() != is_file_little_endian);
  } 
  else if (MB_TAG_NOT_FOUND != rval)
    return rval;
  
    // Compare the number of triangles to the length of the file.  
    // The file must contain an 80-byte description, a 4-byte 
    // triangle count and 50 bytes per triangle.
    //
    // The triangle count *may* allow us to determine the byte order
    // of the file, if it is not an endian-symetric value.  
    //
    // We need to compare the expected size calculated from the triangle
    // count with the file size anyway, as an invalid file or a byte-
    // swapping issue could result in a very large (incorrect) value for
    // num_tri, resulting in a SEGFAULT. 
  
    // Get expected number of triangles
  unsigned long num_tri = swap_bytes ? byte_swap(header.count) : header.count;
  
    // Get the file length
  long filesize = get_file_size( file );
  if (filesize >= 0) // -1 indicates could not determine file size (e.g. reading from FIFO)
  {
      // Check file size, but be careful of numeric overflow
    if (ULONG_MAX / 50 - 84 < num_tri ||   // next calc would have oveflow
        84 + 50 * num_tri != (unsigned long)filesize)
    {
        // Unless the byte order was specified explicitly in the 
        // tag, try the opposite byte order.
      unsigned long num_tri_swap = byte_swap( (uint32_t)num_tri );
      if (bo_tag || // If byte order was specified in tag, fail now
          ULONG_MAX / 50 - 84 < num_tri_swap  || // watch for overflow in next line
          84 + 50 * num_tri_swap != (unsigned long)filesize)
      {
        fclose( file );
        readMeshIface->report_error( "%s: not a binary STL file\n", name );
        return MB_FILE_DOES_NOT_EXIST;
      }
      swap_bytes = !swap_bytes;
      num_tri = num_tri_swap;
    }
  }
    
    // Allocate storage for triangles
  tris.resize( num_tri );
  
    // Read each triangle
  BinaryTri tri;  // binary block read from file
  for (std::vector<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i)
  {
    if (fread( &tri, 50, 1, file ) != 1)
    {
      fclose( file );
      readMeshIface->report_error( "%s: %s\n", name, strerror(errno) );
      return MB_FILE_WRITE_ERROR;
    }
    
      // Swap bytes if necessary
    for (unsigned j = 0; j < 9; ++j)
      if (swap_bytes)
        i->points[j/3].coords[j%3] = byte_swap( tri.coords[j] );
      else
        i->points[j/3].coords[j%3] = tri.coords[j];
  }
  
  fclose( file );
  return MB_SUCCESS;
}


MBReaderIface* ReadSTL::ascii_instance( MBInterface* iface )
  { return new ReadASCIISTL(iface); }

MBReaderIface* ReadSTL::binary_instance( MBInterface* iface )
  { return new ReadBinarySTL(iface); }
