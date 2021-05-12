/** @example kdtree_arborx_test.cpp
 * \brief This example shows how to perform local point-in-element searches with MOAB's
 * kdtree searching functionality, and with newer GPU capable ArborX  searching functionality.
 *
 * MOAB's SpatialLocator functionality performs point-in-element searches over a local or parallel
 * mesh. SpatialLocator is flexible as to what kind of tree is used and what kind of element basis
 * functions are used to localize elements and interpolate local fields.
 *
 * Usage:
 *  Default mesh: ./kdtreeApp -d 3 -r 2
 *  Custom mesh: ./kdtreeApp -d 3 -i mesh_file.h5m
 */

// C++ includes
#include <iostream>
#include <chrono>

#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/CN.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/NestedRefine.hpp"

#include "moab/AdaptiveKDTree.hpp"
#include "TestUtil.hpp"

#include <ArborX.hpp>
#include <ArborX_Version.hpp>
#include <Kokkos_Core.hpp>
#include <random>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;
using MemorySpace = ExecutionSpace::memory_space;

using namespace moab;
using namespace std;

typedef std::chrono::high_resolution_clock Clock;
typedef std::chrono::high_resolution_clock::time_point Timer;
using std::chrono::duration_cast;

Timer start;
std::map< std::string, std::chrono::nanoseconds > timeLog;


#define PUSH_TIMER()          \
    {                         \
        start = Clock::now(); \
    }

#define POP_TIMER( EventName )                                                                                \
    {                                                                                                         \
        std::chrono::nanoseconds elapsed = duration_cast< std::chrono::nanoseconds >( Clock::now() - start ); \
        timeLog[EventName]               = elapsed;                                                           \
    }

#define PRINT_TIMER( EventName )                                                                                       \
    {                                                                                                                  \
        std::cout << "[ " << EventName                                                                                 \
                  << " ]: elapsed = " << static_cast< double >( timeLog[EventName].count() / 1e6 ) << " milli-seconds" \
                  << std::endl;                                                                                        \
    }

int main( int argc, char** argv )
{

    PUSH_TIMER()

    Kokkos::initialize(argc, argv);

    POP_TIMER( "Init k" )
    PRINT_TIMER( "Init k" )
    {
    std::string test_file_name  = TestDir + string( "/3k-tri-sphere.vtk" );
    int num_queries             = 1000000;
    int dimension               = 2;
    int uniformRefinementLevels = 0;
    int max_depth = 30;
    int max_per_leaf = 6;
    std::ostringstream options;
    bool compare_results = false;
    int num_verif = 10;
    bool skip_moab_kdtree = false;
    //
    srand(1234); // set a seed for repeatability


    ProgOptions opts;

    // Need option handling here for input filename
    opts.addOpt< int >( "dim,d", "Dimension of the problem and mesh (default=2)", &dimension );
    opts.addOpt< int >( "queries,n", "Number of queries to perform on the mesh (default=1E4)", &num_queries );
    opts.addOpt< int >( "refine,r", "Number of levels of uniform refinements to perform on the mesh (default=0)",
                        &uniformRefinementLevels );
    opts.addOpt< int > ("max_depth,m", "Maximum depth in kdtree (default 30)", &max_depth);
    opts.addOpt <int >("max_per_leaf,l", "Maximum per leaf (default 6) ", &max_per_leaf);
    opts.addOpt< std::string >( "file,i", "File name to load the mesh)", &test_file_name );
    opts.addOpt <void> ("compare_results,c", "compare search results ", &compare_results);
    opts.addOpt <int> ("verifications,v", "number of verifications", &num_verif);
    opts.addOpt <void> ("skip_moab,s", "skip regular moab kdtree build", &skip_moab_kdtree);

    opts.parseCommandLine( argc, argv );

    if (compare_results)
    {
        if (skip_moab_kdtree)
        {
            cout << " cannot skip moab kdtree if we want to compare results \n";
            skip_moab_kdtree = false;
        }

    }
    // Instantiate
    ErrorCode rval;
    Core mb;
    EntityHandle baseFileset, fileset;

    rval = mb.create_meshset( MESHSET_SET, baseFileset );MB_CHK_ERR( rval );

    PUSH_TIMER()
    // Load the file
    rval = mb.load_file( test_file_name.c_str(), &baseFileset );MB_CHK_SET_ERR( rval, "Error loading file" );

    cout << "load file:  " << test_file_name << "\n";
    if( uniformRefinementLevels )
    {
        moab::NestedRefine uref( &mb, nullptr, baseFileset );
        std::vector< int > uniformRefinementDegree( uniformRefinementLevels, 2 );
        std::vector< EntityHandle > level_sets;
        rval = uref.generate_mesh_hierarchy( uniformRefinementLevels,
                                                uniformRefinementDegree.data(),
                                                level_sets, true );MB_CHK_ERR( rval );
        assert( (int)level_sets.size() == uniformRefinementLevels + 1 );
        fileset = level_sets[uniformRefinementLevels];
    }
    else
        fileset = baseFileset;
    POP_TIMER( "MeshIO-Refine" )
    PRINT_TIMER( "MeshIO-Refine" )

    // Get all elements in the file
    Range elems;
    rval = mb.get_entities_by_dimension( fileset, dimension, elems );MB_CHK_SET_ERR( rval, "Error getting 3d elements" );

    cout << "mesh has " << elems.size()  << " cells \n";



    CartVect min, max, box_extents, pos;
    BoundBox box;
    AdaptiveKDTree kd( &mb ); // constructor does nothing
    if (!skip_moab_kdtree)
    {
        PUSH_TIMER()
        // Create a moab kdtree to use for the location service
        EntityHandle tree_root;
        rval = mb.create_meshset(MESHSET_SET, tree_root); MB_CHK_ERR( rval );
        options << "MAX_DEPTH=" << max_depth << ";";
        options << "MAX_PER_LEAF=" << max_per_leaf << ";";

        FileOptions file_opts( options.str().c_str() );
        rval                   = kd.build_tree( elems, &tree_root, &file_opts );MB_CHK_ERR( rval );
        POP_TIMER( "MOAB KdTree-Setup" )
        PRINT_TIMER( "MOAB KdTree-Setup" )

        // Get the box extents

        unsigned  depth;
        kd.get_info( tree_root, &min[0], &max[0], depth );

        cout << "box: " << min << " " <<  max << "depth: " << depth << "\n";
        // print various info about the tree
        kd.print();
    }
    else
    {
        min = CartVect(-100, -100, -100);
        max = CartVect( 100,  100,  100);
    }
    box = BoundBox (min,max);
    box_extents  = 1.1 * (max - min);
    std::vector<CartVect> queries;
    queries.reserve(num_queries);
    for( int i = 0; i < num_queries; i++ )
    {
        pos  = box.bMin + CartVect( box_extents[0] * .01 * ( rand() % 100 ), box_extents[1] * .01 * ( rand() % 100 ),
                                   box_extents[2] * .01 * ( rand() % 100 ) );
        // project on sphere of radius 100
        pos = pos/pos.length()*100;
        queries.push_back(pos);
    }



    std::vector<EntityHandle> leaves;
    if (compare_results)
        leaves.resize(num_queries);
    if (!skip_moab_kdtree)
    {
        PUSH_TIMER()
        EntityHandle leaf_out;
        for( int i = 0; i < num_queries; i++ )
        {
            pos  = queries[i];

            // project on a sphere if dim 2?
            rval = kd.point_search( &pos[0], leaf_out);  MB_CHK_ERR( rval );
            if (compare_results)
                leaves[i] = leaf_out; // store it for later comparison
        }
        POP_TIMER( "KdTree-Query" )
        PRINT_TIMER( "KdTree-Query" )
    }



    std::cout << "Execution space: " << ExecutionSpace::name() << "\n";
    //Kokkos::ScopeGuard guard(argc, argv);

    std::cout << "ArborX version: " << ArborX::version() << std::endl;
    std::cout << "ArborX hash   : " << ArborX::gitCommitHash() << std::endl;

    // start ArborX build and query
    PUSH_TIMER()
    // Create the View for the bounding boxes, on device
    Kokkos::View<ArborX::Box*, MemorySpace> bounding_boxes("bounding_boxes", elems.size());

    POP_TIMER( "memory space allocate" )
    PRINT_TIMER( "memory space allocate" )

    PUSH_TIMER()
    // with MemorySpace=Kokkos::CudaSpace, BoundingVolume=ArborX::Box, Enable=void
    //Kokkos::View<ArborX::Box*> bounding_boxes("bounding_boxes", elems.size());
    // mirror view on host, will be populated
    auto h_bounding_boxes = Kokkos::create_mirror_view(bounding_boxes);

    std::vector<CartVect>  coords;
    coords.resize(27);// max possible
    for (size_t i=0; i<elems.size(); i++)
    {
        EntityHandle cell=elems[i];
        const EntityHandle *conn ;
        int nnodes;
        rval = mb.get_connectivity(cell, conn, nnodes); MB_CHK_ERR( rval );
        rval = mb.get_coords(conn, nnodes, &coords[0][0] );MB_CHK_SET_ERR( rval, "can't get coordinates" );
        for (int j=0; j<nnodes; j++)
        {
            ArborX::Details::expand(h_bounding_boxes(i),
                    ArborX::Point{ (float)coords[j][0], (float)coords[j][1], (float)coords[j][2]}  );
        }
    }

    POP_TIMER( "Pre ArborX-Build" )
    PRINT_TIMER( "Pre ArborX-Build" )

    PUSH_TIMER()

    Kokkos::deep_copy(bounding_boxes, h_bounding_boxes);
    POP_TIMER( "Deep copy" )
    PRINT_TIMER( "Deep copy" )

    PUSH_TIMER()

    // Create the bounding volume hierarchy
    ArborX::BVH<MemorySpace> bvh(ExecutionSpace{}, bounding_boxes);


    POP_TIMER( "ArborX-Build" )
    PRINT_TIMER( "ArborX-Build" )

    PUSH_TIMER()

    // queries

    // Create the View for the spatial-based queries
    // Kokkos::View<decltype(ArborX::intersects(ArborX::Box{})) *, DeviceType>
    Kokkos::View< decltype(ArborX::intersects(ArborX::Box{})) * , MemorySpace> queries_ar("queries", num_queries);
    // Fill in the queries on host mirror, then copy to device

    auto h_queries_ar = create_mirror_view(queries_ar);

    for( int i = 0; i < num_queries; i++ )
    {
        pos  = queries[i];
        h_queries_ar(i) = ArborX::intersects(
                                                 ArborX::Box (
                                                 { (float)pos[0], (float)pos[1], (float)pos[2] } ,
                                                 { (float)pos[0], (float)pos[1], (float)pos[2] }
                                                 )
                                              );
    }
    // copy from host to device
    Kokkos::deep_copy(queries_ar, h_queries_ar);

    // Perform the search
    Kokkos::View<int*, ExecutionSpace> offsets("offset", 0);
    Kokkos::View<int*, ExecutionSpace> indices("indices", 0);
    ArborX::query(bvh, ExecutionSpace{}, queries_ar, indices, offsets);

    // create mirror and copy at the same time does not work, somehow
    auto host_indices=Kokkos::create_mirror_view( indices);
    auto host_offsets=Kokkos::create_mirror_view( offsets);
    Kokkos::deep_copy(host_indices, indices);
    Kokkos::deep_copy(host_offsets, offsets);

    cout<< "host_indices.size() " << host_indices.size() << "\n";
    cout<< "host_offsets.size() " << host_offsets.size() << "\n";


    POP_TIMER( "ArborX-Query" )
    PRINT_TIMER( "ArborX-Query" )

    cout << bvh.size() << "\n";

    ArborX::Box bb = bvh.bounds();
    ArborX::Point maxcorn = bb.maxCorner();
    cout << "max_corner: " << maxcorn[0] <<" " <<  maxcorn[1] << " "  << maxcorn[2] << "\n";

    if (compare_results)
    {
        // dump the vtk file with the refined mesh, if refinement >=1
        // for few queries, print out the results, and dump the leaves

        if (uniformRefinementLevels)
        {
            std::stringstream ffs;
            ffs << "refinedMesh_" << uniformRefinementLevels << ".vtk";
            rval = mb.write_mesh( ffs.str().c_str(), &fileset, 1 );MB_CHK_ERR( rval );
        }
        for (int i=0; i<num_verif; i++)
        {
            // write an mb file with just the result leaves
            std::stringstream ffs;
            ffs << "leaf_" << i << ".vtk";
            rval = mb.write_mesh( ffs.str().c_str(), &leaves[i], 1 );MB_CHK_ERR( rval );
            // put the results from ArborX in a set to be dumped out
            Range close_by;
            std::cout <<"query point: " << queries[i][0] << " " << queries[i][1] << " " << queries[i][2] << "\n";
            for (int j=host_offsets(i); j<host_offsets(i+1); j++)
            {
                int index = host_indices(j);
                close_by.insert(elems[index]);
                std::cout << "  " << j << ",  " <<  host_indices(j) << "|"  ;
            }
            std::cout << " \n";
            EntityHandle new_set;
            rval = mb.create_meshset( MESHSET_SET, new_set );MB_CHK_ERR( rval );
            rval = mb.add_entities(new_set, close_by);MB_CHK_ERR( rval );
            std::stringstream fff;
            fff << "close_by_" << i << ".vtk";
            rval = mb.write_mesh( fff.str().c_str(), &new_set, 1 );MB_CHK_ERR( rval );
            // write the search point into a new vertex
            rval = mb.remove_entities(new_set, close_by);MB_CHK_ERR( rval );
            EntityHandle new_vertex;
            rval = mb.create_vertex(&queries[i][0], new_vertex); MB_CHK_ERR( rval );
            std::stringstream ffp;
            rval = mb.add_entities(new_set, &new_vertex, 1 ); MB_CHK_ERR( rval );
           ffp << "search_point_" << i << ".vtk";
           rval = mb.write_mesh( ffp.str().c_str(), &new_set, 1 );MB_CHK_ERR( rval );
        }
    }

    //
    }
    Kokkos::finalize();


    return 0;
}

