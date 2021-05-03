//
// Usage:
// tools/mbinferredpart -d 2 --spherical -m mpas/x1.2562.grid.h5m -s mpas/x1.10242.grid.h5m -o mpas_inferred.h5m -e 1e-8 -b 1e-6 
//
#include <iostream>
#include <exception>
#include <cmath>
#include <vector>
#include <string>

#include "moab/ProgOptions.hpp"
#include "moab/Core.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "moab/BVHTree.hpp"

#include "moab/IntxMesh/IntxUtils.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#endif

using namespace moab;

// A function to get the non-default value from a std::map
template < typename K, typename V >
static V get_map_value( const std::map< K, V >& m, const K& key, const V& defval )
{
    typename std::map< K, V >::const_iterator it = m.find( key );
    if( it == m.end() ) { return defval; }
    else
    {
        return it->second;
    }
}

int main( int argc, char* argv[] )
{
    int nprocs = 1, dimension = 3;
    int dimension_source = 3;
#ifdef MOAB_HAVE_MPI
    int proc_id = 0;
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &proc_id );
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
#endif

    int defaultpart  = 0;
    double tolerance = 1e-2, treetolerance = 1e-13, btolerance = 1e-4;
    std::string primaryfile, secondaryfile, outfile( "secondarymesh.h5m" );
    bool keepsparts    = false;
    bool use_spherical = false;
    ProgOptions opts;

    opts.addOpt< std::string >( "primary,m", "Primary mesh filename", &primaryfile );
    opts.addOpt< std::string >( "secondary,s", "Secondary mesh filename", &secondaryfile );
    opts.addOpt< std::string >( "output,o", "Output partitioned mesh filename", &outfile );
    opts.addOpt< int >( "dims,c", "Dimension of entities in source used for partitioning", &dimension_source );
    opts.addOpt< int >( "dim,d", "Dimension of entities in target to use for partitioning", &dimension );
    opts.addOpt< int >( "defaultpart,p", "Default partition number if target element is not found on source grid",
                        &defaultpart );
    opts.addOpt< double >( "eps,e", "Tolerance for the point search", &tolerance );
    opts.addOpt< double >( "beps,b", "Tolerance for the bounding box search", &btolerance );
    opts.addOpt< void >( "keep,K",
                         "Keep the existing partitions in the secondary mesh (use PARALLEL_PARTITION_secondary instead)",
                         &keepsparts );
    opts.addOpt< void >( "spherical", "Hint that the meshes are defined on a spherical surface (Climate problems)",
                         &use_spherical );
    opts.parseCommandLine( argc, argv );

    if( primaryfile.empty() || secondaryfile.empty() )
    {
        opts.printHelp();
#ifdef MOAB_HAVE_MPI
        MPI_Finalize();
#endif
        exit( 1 );
    }

    ErrorCode error;
    Core* mbCore = new Core();

    // Set the read options for parallel file loading
    const std::string partition_set_name = "PARALLEL_PARTITION";
    const std::string global_id_name     = "GLOBAL_ID";
    std::vector< std::string > read_opts, write_opts;
    std::string read_options( "" ), write_options( "" );

    if( nprocs > 1 )
    {
        read_options  = "PARALLEL=READ_PART;PARTITION=" + partition_set_name + ";PARALLEL_RESOLVE_SHARED_ENTS";
        write_options = "PARALLEL=WRITE_PART";
    }

    EntityHandle primaryfileset, secondaryfileset;
    error = mbCore->create_meshset( moab::MESHSET_TRACK_OWNER | moab::MESHSET_SET, primaryfileset );MB_CHK_ERR( error );
    error = mbCore->create_meshset( moab::MESHSET_TRACK_OWNER | moab::MESHSET_SET, secondaryfileset );MB_CHK_ERR( error );

    // Load file
    error = mbCore->load_file( primaryfile.c_str(), &primaryfileset, read_options.c_str() );MB_CHK_ERR( error );
    error = mbCore->load_file( secondaryfile.c_str(), &secondaryfileset, read_options.c_str() );MB_CHK_ERR( error );
    // if (error != MB_SUCCESS && size > 1)
    // {
    //   std::string newread_options = "PARALLEL=BCAST_DELETE;PARALLEL_RESOLVE_SHARED_ENTS";
    //   error = mbCore->load_file(secondaryfile.c_str(), &secondaryfileset, newread_options.c_str());
    // }
    // else MB_CHK_ERR(error);

    Tag gidtag = 0, parttag = 0, sparttag = 0;
    int dum_id = -1;
    error      = mbCore->tag_get_handle( partition_set_name.c_str(), 1, MB_TYPE_INTEGER, parttag,
                                    MB_TAG_SPARSE | MB_TAG_CREAT, &dum_id );MB_CHK_ERR( error );
    gidtag = mbCore->globalId_tag();
    if( keepsparts )
    {
        error = mbCore->tag_get_handle( std::string( partition_set_name + "_secondary" ).c_str(), 1, MB_TYPE_INTEGER,
                                        sparttag, MB_TAG_CREAT | MB_TAG_SPARSE, &dum_id );MB_CHK_ERR( error );
    }

    Range melems, msets, selems, ssets;

    // Get the partition sets on the primary mesh
    std::map< int, int > mpartvals;
    error = mbCore->get_entities_by_type_and_tag( primaryfileset, MBENTITYSET, &parttag, NULL, 1, msets,
                                                  moab::Interface::UNION, true );MB_CHK_ERR( error );
    if( msets.size() == 0 )
    {
        std::cout << "No partition sets found in the primary mesh. Quitting..." << std::endl;
        exit( 1 );
    }

    for( unsigned i = 0; i < msets.size(); ++i )
    {
        EntityHandle mset = msets[i];

        moab::Range msetelems;
        error = mbCore->get_entities_by_dimension( mset, dimension_source, msetelems );MB_CHK_ERR( error );
        melems.merge( msetelems );

        int partID;
        error = mbCore->tag_get_data( parttag, &mset, 1, &partID );MB_CHK_ERR( error );

        // Get the global ID and use that as the indicator
        std::vector< int > gidMelems( msetelems.size() );
        error = mbCore->tag_get_data( gidtag, msetelems, gidMelems.data() );MB_CHK_ERR( error );

        for( unsigned j = 0; j < msetelems.size(); ++j )
            mpartvals[gidMelems[j]] = partID;
            // mpartvals[msetelems[j]]=partID;
#ifdef VERBOSE
        std::cout << "Part " << partID << " has " << msetelems.size() << " elements." << std::endl;
#endif
    }

    // Get information about the secondary file set
    error = mbCore->get_entities_by_type_and_tag( secondaryfileset, MBENTITYSET, &parttag, NULL, 1, ssets,
                                                  moab::Interface::UNION );MB_CHK_ERR( error );
    // TODO: expand and add other dimensional elements
    error = mbCore->get_entities_by_dimension( secondaryfileset, dimension, selems );MB_CHK_ERR( error );

    std::cout << "Primary (elements, parts) : (" << melems.size() << ", " << msets.size()
              << "), Secondary (elements, parts) : (" << selems.size() << ", " << ssets.size() << ")" << std::endl;

    double primary_radius = 1.0, secondary_radius = 1.0;
    std::vector< double > primarycoords;
    Range primaryverts, secondaryverts;
    {
        error = mbCore->get_entities_by_dimension( primaryfileset, 0, primaryverts );MB_CHK_ERR( error );
        error = mbCore->get_entities_by_dimension( secondaryfileset, 0, secondaryverts );MB_CHK_ERR( error );
    }
    if( use_spherical )
    {
        double points[6];
        EntityHandle mfrontback[2] = { primaryverts[0], primaryverts[primaryverts.size() - 1] };
        error                      = mbCore->get_coords( &mfrontback[0], 2, points );MB_CHK_ERR( error );
        primary_radius = 0.5 * ( std::sqrt( points[0] * points[0] + points[1] * points[1] + points[2] * points[2] ) +
                                std::sqrt( points[3] * points[3] + points[4] * points[4] + points[5] * points[5] ) );
        EntityHandle sfrontback[2] = { secondaryverts[0], secondaryverts[secondaryverts.size() - 1] };
        error                      = mbCore->get_coords( &sfrontback[0], 2, points );MB_CHK_ERR( error );
        secondary_radius = 0.5 * ( std::sqrt( points[0] * points[0] + points[1] * points[1] + points[2] * points[2] ) +
                               std::sqrt( points[3] * points[3] + points[4] * points[4] + points[5] * points[5] ) );
        if( fabs( primary_radius - secondary_radius ) > 1e-14 )
        {
            // Let us rescale both primary and secondary meshes to a unit sphere
            error = moab::IntxUtils::ScaleToRadius( mbCore, primaryfileset, 1.0 );MB_CHK_ERR( error );
            error = moab::IntxUtils::ScaleToRadius( mbCore, secondaryfileset, 1.0 );MB_CHK_ERR( error );
        }
    }

    try
    {
        std::map< int, moab::Range > spartvals;
        int npoints_notfound = 0;
        {
            FileOptions fopts( ( use_spherical ? "SPHERICAL" : "" ) );
            moab::EntityHandle tree_root = primaryfileset;

            moab::AdaptiveKDTree tree( mbCore );
            error = tree.build_tree( melems, &tree_root, &fopts );

            // moab::BVHTree tree(mbCore);
            // error = tree.build_tree(melems, &tree_root);MB_CHK_ERR(error);

            for( size_t ie = 0; ie < selems.size(); ie++ )
            {
                moab::EntityHandle selem, leaf;
                double point[3];
                selem = selems[ie];

                // Get the element centroid to be queried
                error = mbCore->get_coords( &selem, 1, point );MB_CHK_ERR( error );

                std::vector< moab::EntityHandle > leaf_elems;

                // Search for the closest source element in the primary mesh corresponding
                // to the target element centroid in the secondary mesh
                error = tree.point_search( point, leaf, treetolerance, btolerance );MB_CHK_ERR( error );

                // We only care about the dimension that the user specified.
                // MOAB partitions are ordered by elements anyway.
                error = mbCore->get_entities_by_dimension( leaf, dimension_source, leaf_elems, true );MB_CHK_ERR( error );

                if( leaf != 0 && leaf_elems.size() )
                {

                    // Now get the primary element centroids so that we can compute
                    // the minimum distance to the target point
                    std::vector< double > centroids( leaf_elems.size() * 3 );
                    error = mbCore->get_coords( &leaf_elems[0], leaf_elems.size(), &centroids[0] );MB_CHK_ERR( error );

                    if( !leaf_elems.size() )
                        std::cout << ie << ": "
                                  << " No leaf elements found." << std::endl;
                    // else std::cout << ie << " found " << leaf_elems.size() << " leaves for
                    // current point " << std::endl;

                    double dist = 1e10;
                    int pinelem = -1;
                    for( size_t il = 0; il < leaf_elems.size(); ++il )
                    {
                        const double* centroid = &centroids[il * 3];
                        const double locdist   = std::pow( point[0] - centroid[0], 2 ) +
                                               std::pow( point[1] - centroid[1], 2 ) +
                                               std::pow( point[2] - centroid[2], 2 );
                        if( locdist < dist && locdist < 1.0E-2 )
                        {
                            dist    = locdist;
                            pinelem = il;

#ifdef VERBOSE
                            int gidMelem;
                            error = mbCore->tag_get_data( gidtag, &leaf_elems[il], 1, &gidMelem );MB_CHK_ERR( error );
                            std::cout << "\t Trial leaf " << il << " set " << gidMelem
                                      << " and part = " << get_map_value( mpartvals, gidMelem, -1 )
                                      << " with distance = " << locdist << std::endl;
#endif
                        }
                    }

                    if( pinelem < 0 )
                    {
#ifdef VERBOSE
                        std::cout << ie
                                  << ": [Error] - Could not find a minimum distance within the "
                                     "leaf nodes."
                                  << std::endl;
#endif
                        npoints_notfound++;
                        spartvals[defaultpart].insert( selems[ie] );
                    }
                    else
                    {
                        int gidMelem;
                        error = mbCore->tag_get_data( gidtag, &leaf_elems[pinelem], 1, &gidMelem );MB_CHK_ERR( error );

                        int mpartid = get_map_value( mpartvals, gidMelem, -1 );
                        if( mpartid < 0 )
                            std::cout << "[WARNING]: Part number for element " << leaf_elems[pinelem]
                                      << " with global ID = " << gidMelem << " not found.\n";

#ifdef VERBOSE
                        std::cout << "Adding element " << ie << " set " << mpartid << " with distance = " << dist
                                  << std::endl;
#endif
                        spartvals[mpartid].insert( selems[ie] );
                    }
                }
                else
                {
#ifdef VERBOSE
                    std::cout << "[WARNING]: Adding element " << ie << " to default (part) set " << defaultpart
                              << std::endl;
#endif
                    spartvals[defaultpart].insert( selems[ie] );
                }
            }
            error = tree.reset_tree();MB_CHK_ERR( error );
        }
        if( npoints_notfound )
            std::cout << "Could not find " << npoints_notfound
                      << " points in the primary mesh. Added to defaultpart set = " << defaultpart << std::endl;

        if( use_spherical )
        {
            if( fabs( primary_radius - secondary_radius ) > 1e-14 )
            {
                error = moab::IntxUtils::ScaleToRadius( mbCore, secondaryfileset, secondary_radius );MB_CHK_ERR( error );
            }
        }

        error = mbCore->delete_entities( &primaryfileset, 1 );MB_CHK_ERR( error );
        // Find parallel partition sets in the secondary mesh - and delete it since we are going to
        // overwrite the sets
        if( !keepsparts )
        {
            std::cout << "Deleting " << ssets.size() << " sets in the secondary mesh" << std::endl;
            error = mbCore->remove_entities( secondaryfileset, ssets );MB_CHK_ERR( error );
            ssets.clear();
        }

        size_t ntotsecondary_elems = 0, ntotsecondary_parts = 0;
        for( std::map< int, moab::Range >::iterator it = spartvals.begin(); it != spartvals.end(); ++it )
        {
            int partID = it->first;
            moab::EntityHandle pset;
            error = mbCore->create_meshset( moab::MESHSET_SET, pset );MB_CHK_ERR( error );
            error = mbCore->add_entities( pset, it->second );MB_CHK_ERR( error );
            error = mbCore->add_parent_child( secondaryfileset, pset );MB_CHK_ERR( error );

#ifdef VERBOSE
            std::cout << "Secondary Part " << partID << " has " << it->second.size() << " elements." << std::endl;
#endif
            ntotsecondary_elems += it->second.size();
            ntotsecondary_parts++;

            if( keepsparts )
            {
                error = mbCore->tag_set_data( sparttag, &pset, 1, &partID );MB_CHK_ERR( error );
            }
            else
            {
                error = mbCore->tag_set_data( parttag, &pset, 1, &partID );MB_CHK_ERR( error );
            }
        }
        std::cout << "Secondary mesh: given " << selems.size() << " elements, and assigned " << ntotsecondary_elems
                  << " elements to " << ntotsecondary_parts << " parts.\n";
        assert( ntotsecondary_elems == selems.size() );

        // mbCore->print_database();

        // Write the re-partitioned secondary mesh to disk
        if( nprocs == 1 )
        {
            error = mbCore->write_file( "secondarymesh.vtk", "VTK", NULL, &secondaryfileset, 1 );MB_CHK_ERR( error );
        }
        error = mbCore->write_file( outfile.c_str(), NULL, write_options.c_str(), &secondaryfileset, 1 );MB_CHK_ERR( error );
        // error = mbCore->write_file(outfile.c_str(), NULL,
        // write_options.c_str());MB_CHK_ERR(error);
    }
    catch( std::exception& e )
    {
        std::cout << " exception caught during tree initialization " << e.what() << std::endl;
    }
    delete mbCore;

#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif
    exit( 0 );
}
