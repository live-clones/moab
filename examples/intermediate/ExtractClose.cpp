/** \brief This test shows how to extract mesh from a model, based on distance.
 *
 * MOAB's It is needed to extract from large mesh files cells close to some point, where we suspect
 * errors It would be useful for 9Gb input file that we cannot visualize, but we have overlapped
 * elements.
 */

#include <iostream>
#include <cstdlib>
#include <cstdio>

#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/CartVect.hpp"
#include "moab/IntxMesh/IntxUtils.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#endif

using namespace moab;
using namespace std;

#ifdef MOAB_HAVE_HDF5
string test_file_name = string( MESH_DIR ) + string( "/64bricks_512hex_256part.h5m" );
#endif
int main( int argc, char** argv )
{
    ProgOptions opts;

    // vertex 25 in filt out  (-0.784768, 0.594952, -0.173698)
    double x = -0.784768, y = 0.594952, z = -0.173698;
    double lat = 28.6551, lon = 61.0204;  // from Charlie Zender's page
    // https://acme-climate.atlassian.net/wiki/spaces/ED/pages/1347092547/Assessing+and+Addressing+Regridding+Weights+in+Map-Files

    string inputFile = test_file_name;
    opts.addOpt< string >( "inFile,i", "Specify the input file name string ", &inputFile );

    string outFile = "out.h5m";
    opts.addOpt< string >( "outFile,o", "Specify the output file name string ", &outFile );

    opts.addOpt< double >( string( "xpos,x" ), string( "x position" ), &x );
    opts.addOpt< double >( string( "ypos,y" ), string( "y position" ), &y );
    opts.addOpt< double >( string( "zpos,z" ), string( "z position" ), &z );

    opts.addOpt< double >( string( "latitude,t" ), string( "north latitude in degrees" ), &lat );
    opts.addOpt< double >( string( "longitude,w" ), string( "east longitude in degrees" ), &lon );

    bool spherical = false;
    opts.addOpt< void >( "spherical,s", "use spherical coords input (default false) ", &spherical );

    double distance = 0.01;
    opts.addOpt< double >( string( "distance,d" ), string( "distance " ), &distance );

    opts.parseCommandLine( argc, argv );

    int         rank = 0;
    int         numProcesses = 1;
    std::string readopts;

#ifdef MOAB_HAVE_MPI
    int fail = MPI_Init( &argc, &argv );
    if( fail ) return 1;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &numProcesses );

    readopts = string( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION" );  // we do not have to
                                                                             // resolve shared ents
#endif
    // Instantiate
    Core mb;

    // Load the file into a new file set
    EntityHandle fileSet;
    ErrorCode    rval = mb.create_meshset( MESHSET_SET, fileSet );MB_CHK_SET_ERR( rval, "Error creating file set" );
    rval = mb.load_file( inputFile.c_str( ), &fileSet, readopts.c_str( ) );MB_CHK_SET_ERR( rval, "Error loading file" );

    if( !rank )
    {
        cout << " reading file " << inputFile << " on " << numProcesses << " task";
        if( numProcesses > 1 ) cout << "s";
        cout << "\n";
    }
    // Get all 2d elements in the file set
    Range elems;
    rval = mb.get_entities_by_dimension( fileSet, 2, elems );MB_CHK_SET_ERR( rval, "Error getting 2d elements" );

    // create a meshset with close elements
    EntityHandle outSet;

    // create meshset
    rval = mb.create_meshset( MESHSET_SET, outSet );MB_CHK_ERR( rval );

    // double sphere radius is 1
    CartVect point( x, y, z );
    if( spherical )
    {
        std::cout << " use spherical coordinates on the sphere of radius 1.\n";
        SphereCoords pos;
        pos.R = 1;  // always on the
        pos.lat = lat * M_PI / 180;  // defined in math.h
        pos.lon = lon * M_PI / 180;
        point = spherical_to_cart( pos );
    }

    Range closeByCells;
    for( Range::iterator it = elems.begin( ); it != elems.end( ); it++ )
    {
        EntityHandle cell = *it;
        CartVect     center;
        rval = mb.get_coords( &cell, 1, &( center[ 0 ] ) );MB_CHK_SET_ERR( rval, "Can't get cell center coords" );
        double dist = ( center - point ).length( );
        if( dist <= distance ) { closeByCells.insert( cell ); }
    }

    rval = mb.add_entities( outSet, closeByCells );MB_CHK_SET_ERR( rval, "Can't add to entity set" );

    int numCells = (int)closeByCells.size( );
#ifdef MOAB_HAVE_MPI
    // Reduce all of the local sums into the global sum
    int globalCells;
    MPI_Reduce( &numCells, &globalCells, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
    numCells = globalCells;
#endif
    if( numCells == 0 )
    {
        if( !rank )
            std::cout << " no close cells to the point " << point << " at distance less than " << distance << "\n";
    }
    else
    {
        if( !rank )
            std::cout << " write file " << outFile << " with cells closer than " << distance << " from " << point
                      << "\n";
        string writeOpts;
        if( numProcesses > 1 ) writeOpts = string( "PARALLEL=WRITE_PART;" );
        rval = mb.write_file( outFile.c_str( ), 0, writeOpts.c_str( ), &outSet, 1 );MB_CHK_SET_ERR( rval, "Can't write file" );
    }
    return 0;
}
