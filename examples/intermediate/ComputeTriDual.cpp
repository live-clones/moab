/** @example ComputeTriDual.cpp \n
 * \brief Compute the polygonal (dual) mesh of the primal simplex grid in 2-D. \n
 * <b>To run</b>: ComputeTriDual  input_file output_file \n
 * An example to demonstrate computing the dual polygonal grid of a triangulation 
 * on a spherical surface.
 */

#include <iostream>
#include <vector>

// Include header for MOAB instance and tag conventions for
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/CartVect.hpp"
#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#endif

///////////////////////////////////////////////////////////////////////////////

moab::ErrorCode compute_dual_mesh( moab::Interface* mb, moab::EntityHandle& dual_set, moab::Range& cells )
{
    moab::ErrorCode rval;

    moab::Range verts;
    rval = mb->get_connectivity( cells, verts );MB_CHK_SET_ERR( rval, "Failed to get connectivity" );

    moab::ReadUtilIface* iface;
    rval = mb->query_interface( iface );MB_CHK_SET_ERR( rval, "Can't get reader interface" );

    moab::EntityHandle startv;
    std::vector< double* > ccenters;
    rval = iface->get_node_coords( 3, cells.size(), 0, startv, ccenters );MB_CHK_SET_ERR( rval, "Can't get node coords" );

    moab::Tag gidTag = mb->globalId_tag();

    std::vector< int > gids( cells.size() );
    rval = mb->get_coords( cells, ccenters[0], ccenters[1], ccenters[2] );MB_CHK_SET_ERR( rval, "Failed to get coordinates" );
    rval = mb->tag_get_data( gidTag, cells, &gids[0] );MB_CHK_SET_ERR( rval, "Can't set global_id tag" );

    moab::Range dualverts( startv, startv + cells.size() - 1 );
    rval = mb->add_entities( dual_set, dualverts );MB_CHK_SET_ERR( rval, "Can't add entities" );
    rval = mb->tag_set_data( gidTag, dualverts, &gids[0] );MB_CHK_SET_ERR( rval, "Can't set global_id tag" );

#define CC( ind )        moab::CartVect( ccenters[0][ind], ccenters[1][ind], ccenters[2][ind] )
#define CCXMY( ind, cv ) moab::CartVect( ccenters[0][ind] - cv[0], ccenters[1][ind] - cv[1], ccenters[2][ind] - cv[2] )

    const moab::EntityHandle svtx = *( cells.begin() );

    // maintain block-wise polygon elements
    moab::Range dualcells;
    std::map< size_t, std::vector< moab::EntityHandle > > polygonConn;
    // Generate new Face array
    for( moab::Range::iterator vit = verts.begin(); vit != verts.end(); vit++ )
    {
        moab::EntityHandle vtx = *vit;

        std::vector< moab::EntityHandle > adjs;
        // generate all required adjacencies
        rval = mb->get_adjacencies( &vtx, 1, 2, true, adjs );MB_CHK_SET_ERR( rval, "Failed to get adjacencies" );
        const size_t nEdges = adjs.size();

        double vxyz[3];
        rval = mb->get_coords( &vtx, 1, vxyz );MB_CHK_SET_ERR( rval, "Failed to get coordinates" );

        // Reorient Faces
        moab::CartVect nodeCentral( vxyz );
        moab::CartVect noderef = CC( adjs[0] - svtx );
        moab::CartVect node0   = noderef - nodeCentral;
        double dNode0Mag       = node0.length();

        moab::CartVect nodeCross = noderef * nodeCentral;

        // Determine the angles about the central Node of each Face Node
        std::vector< double > dAngles;
        dAngles.resize( nEdges );
        dAngles[0] = 0.0;

        for( size_t j = 1; j < nEdges; j++ )
        {
            moab::CartVect nodeDiff = CCXMY( adjs[j] - svtx, nodeCentral );
            double dNodeDiffMag     = nodeDiff.length();

            double dSide    = nodeCross % nodeDiff;
            double dDotNorm = ( node0 % nodeDiff ) / ( dNode0Mag * dNodeDiffMag );

            // double dAngle;
            if( dDotNorm > 1.0 ) { dDotNorm = 1.0; }

            dAngles[j] = acos( dDotNorm );

            if( dSide > 0.0 ) { dAngles[j] = -dAngles[j] + 2.0 * M_PI; }
        }

        std::vector< moab::EntityHandle >& face = polygonConn[nEdges];
        if( face.size() == 0 ) face.reserve( nEdges * verts.size() / 4 );
        face.push_back( dualverts[adjs[0] - svtx] );

        // Orient each Face by putting Nodes in order of increasing angle
        double dCurrentAngle = 0.0;
        for( size_t j = 1; j < nEdges; j++ )
        {
            int ixNextNode    = 1;
            double dNextAngle = 2.0 * M_PI;

            for( size_t k = 1; k < nEdges; k++ )
            {
                if( ( dAngles[k] > dCurrentAngle ) && ( dAngles[k] < dNextAngle ) )
                {
                    ixNextNode = k;
                    dNextAngle = dAngles[k];
                }
            }

            face.push_back( dualverts[adjs[ixNextNode] - svtx] );
            dCurrentAngle = dNextAngle;
        }
    }

    std::map< size_t, std::vector< moab::EntityHandle > >::iterator eit;
    for( eit = polygonConn.begin(); eit != polygonConn.end(); ++eit )
    {
        const size_t nVPerE      = eit->first;
        const size_t nElePerType = static_cast< int >( eit->second.size() / nVPerE );

        moab::EntityHandle starte;  // Connectivity
        moab::EntityHandle* conn;

        rval = iface->get_element_connect( nElePerType, eit->first, moab::MBPOLYGON, 0, starte, conn );MB_CHK_SET_ERR( rval, "Can't get element connectivity" );

        // copy the connectivity that we have accumulated
        std::copy( eit->second.begin(), eit->second.end(), conn );
        eit->second.clear();

        // add this polygon sequence to the aggregate data
        moab::Range mbcells( starte, starte + nElePerType - 1 );
        dualcells.merge( mbcells );
    }

    // add the computed dual cells to mesh
    rval = mb->add_entities( dual_set, dualcells );MB_CHK_SET_ERR( rval, "Can't add polygonal entities" );

    // Assign global IDs to all the dual cells - same as original vertices
    assert( dualcells.size() == verts.size() );
    gids.resize( verts.size() );
    rval = mb->tag_get_data( gidTag, verts, &gids[0] );MB_CHK_SET_ERR( rval, "Can't set global_id tag" );
    if( gids[0] == gids[1] && gids[0] < 0 )
    {
#ifdef MOAB_HAVE_MPI
        moab::ParallelComm pcomm( mb, MPI_COMM_WORLD );

        rval = pcomm.assign_global_ids( dual_set, 2, 1, false, true, true );MB_CHK_SET_ERR( rval, "Can't assign global_ids" );
#else
        // No global ID assigned to input vertices.
        // Can we use std::iota ??
        for( size_t ix = 0; ix < gids.size(); ++ix )
            gids[ix] = ix + 1;

        // set GID tag
        rval = mb->tag_set_data( gidTag, dualcells, &gids[0] );MB_CHK_SET_ERR( rval, "Can't set global_id tag" );
        std::cout << "GIDs: " << gids[0] << " " << gids[1] << " " << gids[2] << " " << gids[3] << " " << gids[4] << " "
                  << gids[5] << "\n";
#endif
    }
    gids.clear();

    // delete the original entities from the mesh
    rval = mb->delete_entities( cells );MB_CHK_SET_ERR( rval, "Can't remove entities" );
    rval = mb->delete_entities( verts );MB_CHK_SET_ERR( rval, "Can't remove entities" );

    return moab::MB_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv )
{
    moab::ErrorCode rval;
    ProgOptions opts;

    // Get MOAB instance
    moab::Interface* mb = new( std::nothrow ) moab::Core;
    if( NULL == mb ) return 1;

    std::string inputFile = std::string( MESH_DIR ) + "/globalmpas_deltri.h5m";
    opts.addOpt< std::string >( "inFile,i", "Specify the input file name string ", &inputFile );
#ifdef MOAB_HAVE_HDF5
    std::string outputFile = "outDualMesh.h5m";
#else
    std::string outputFile = "outDualMesh.vtk";
#endif
    opts.addOpt< std::string >( "outFile,o", "Specify the input file name string ", &outputFile );

    opts.parseCommandLine( argc, argv );

    int rank         = 0;
    int numProcesses = 1;

#ifdef MOAB_HAVE_MPI
    if( MPI_Init( &argc, &argv ) ) return 1;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &numProcesses );

    // we need one ghost layer so that each process can compute the dual independently
    const char* readopts = ( numProcesses == 1 ? ""
                                               : std::string( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;"
                                                              "PARALLEL_RESOLVE_SHARED_ENTS;PARALLEL_GHOSTS=2.0.1" )
                                                     .c_str() );
#else
    const char* readopts = "";
#endif

    moab::EntityHandle triangle_set, dual_set;
    rval = mb->create_meshset( moab::MESHSET_SET, triangle_set );MB_CHK_SET_ERR( rval, "Can't create new set" );
    rval = mb->create_meshset( moab::MESHSET_SET, dual_set );MB_CHK_SET_ERR( rval, "Can't create new set" );

    // This file is in the mesh files directory
    rval = mb->load_file( inputFile.c_str(), &triangle_set, readopts );MB_CHK_SET_ERR( rval, "Failed to read" );

    // get all cells of dimension 2;
    moab::Range cells;
    rval = mb->get_entities_by_dimension( triangle_set, 2, cells );MB_CHK_SET_ERR( rval, "Failed to get cells" );

    std::cout << "Original number of triangular cells : " << cells.size() << "\n";

    // call the routine to compute the dual grid
    rval = compute_dual_mesh( mb, dual_set, cells );MB_CHK_SET_ERR( rval, "Failed to compute dual mesh" );

    std::string writeopts="";
#if defined(MOAB_HAVE_MPI) && defined(MOAB_HAVE_HDF5)
    if( outputFile.substr( outputFile.find_last_of( "." ) + 1 ) == "h5m" )
    {
        writeopts = "PARALLEL=WRITE_PART";
    }
#endif

    // write the mesh to disk
    rval = mb->write_file( outputFile.c_str(), 0, writeopts.c_str() );MB_CHK_SET_ERR( rval, "Failed to write new file" );

    cells.clear();
    rval = mb->get_entities_by_dimension( dual_set, 2, cells );MB_CHK_SET_ERR( rval, "Failed to get cells" );
    std::cout << "Wrote the dual mesh: " << outputFile << " containing " << cells.size() << " polygonal cells\n";

    delete mb;

#ifdef MOAB_HAVE_MPI
    if( MPI_Finalize() ) return 1;
#endif
    return 0;
}

