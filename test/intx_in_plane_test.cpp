/*
 * intx_in_plane_test.cpp
 *
 *  Created on: Oct 4, 2012
 */
#include <iostream>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/IntxMesh/Intx2MeshInPlane.hpp"
#include "moab/IntxMesh/IntxUtils.hpp"
#include "TestUtil.hpp"
#include <cmath>

using namespace moab;

int main( int argc, char* argv[] )
{
    // check command line arg
    const char* filename_mesh1 = STRINGIFY( MESHDIR ) "/mbcslam/m1.vtk";
    const char* filename_mesh2 = STRINGIFY( MESHDIR ) "/mbcslam/m2.vtk";
    const char* newFile        = "intx1.vtk";
    const char* edgesFile      = "polyWithEdges.vtk";
    if( argc == 4 )
    {
        filename_mesh1 = argv[1];
        filename_mesh2 = argv[2];
        newFile        = argv[3];
    }
    else
    {
        printf( "Usage: %s <mesh_filename1> <mesh_filename2>  <newFile>\n", argv[0] );
        if( argc != 1 ) return 1;
        printf( "No files specified.  Defaulting to: %s  %s  %s\n", filename_mesh1, filename_mesh2, newFile );
    }

    // read meshes in 2 file sets
    Core moab;
    Interface* mb = &moab;  // global
    EntityHandle sf1, sf2;
    ErrorCode rval = mb->create_meshset( MESHSET_SET, sf1 );
    if( MB_SUCCESS != rval ) return 1;
    rval = mb->create_meshset( MESHSET_SET, sf2 );
    if( MB_SUCCESS != rval ) return 1;
    rval = mb->load_file( filename_mesh1, &sf1 );
    if( MB_SUCCESS != rval ) return 1;
    rval = mb->load_file( filename_mesh2, &sf2 );
    if( MB_SUCCESS != rval ) return 1;

    EntityHandle outputSet;
    rval = mb->create_meshset( MESHSET_SET, outputSet );
    if( MB_SUCCESS != rval ) return 1;

    Intx2MeshInPlane worker( mb );
    IntxAreaUtils areaAdaptor;
    rval = areaAdaptor.positive_orientation( mb, sf1, -1 );MB_CHK_ERR( rval );
    rval = areaAdaptor.positive_orientation( mb, sf2, -1 );MB_CHK_ERR( rval );

    worker.set_error_tolerance( 1.e-5 );
    rval = worker.FindMaxEdges( sf1, sf2 );MB_CHK_ERR( rval );
    // worker.enable_debug();
    rval = worker.intersect_meshes( sf1, sf2, outputSet );
    if( MB_SUCCESS != rval ) return 1;
    rval = mb->write_mesh( newFile, &outputSet, 1 );
    if( MB_SUCCESS != rval ) return 1;

    // retrieve polygons and get adjacent edges
    Range polygons;
    rval = mb->get_entities_by_type( outputSet, MBPOLYGON, polygons );
    if( MB_SUCCESS != rval ) return 1;

    Range edges;
    rval = mb->get_adjacencies( polygons, 1, true, edges, Interface::UNION );
    if( MB_SUCCESS != rval ) return 1;

    std::cout << "number of edges:" << edges.size() << "\n";
    // add edges to the output set
    rval = mb->add_entities( outputSet, edges );
    if( MB_SUCCESS != rval ) return 1;
    rval = mb->write_mesh( edgesFile, &outputSet, 1 );
    if( MB_SUCCESS != rval ) return 1;
    return 0;
}
