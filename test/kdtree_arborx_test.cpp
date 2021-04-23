/** @example kdTreeApp.cpp
 * \brief This example shows how to perform local point-in-element searches with MOAB's new tree
 * searching functionality.
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
    std::string test_file_name  = TestDir + string( "/hex_2048.vtk" );
    int num_queries             = 1000000;
    int dimension               = 3;
    int uniformRefinementLevels = 0;

    ProgOptions opts;
    // Need option handling here for input filename
    opts.addOpt< int >( "dim,d", "Dimension of the problem and mesh (default=3)", &dimension );
    opts.addOpt< int >( "queries,n", "Number of queries to perform on the mesh (default=1E4)", &num_queries );
    opts.addOpt< int >( "refine,r", "Number of levels of uniform refinements to perform on the mesh (default=0)",
                        &uniformRefinementLevels );
    opts.addOpt< std::string >( "file,i", "File name to load the mesh)", &test_file_name );

    opts.parseCommandLine( argc, argv );

    // Instantiate
    ErrorCode rval;
    Core mb;
    EntityHandle baseFileset, fileset;

    rval = mb.create_meshset( MESHSET_SET, baseFileset );MB_CHK_ERR( rval );

    PUSH_TIMER()
    // Load the file
    rval = mb.load_file( test_file_name.c_str(), &baseFileset );MB_CHK_SET_ERR( rval, "Error loading file" );

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

    // Get all 3d elements in the file
    Range elems;
    rval = mb.get_entities_by_dimension( fileset, dimension, elems );MB_CHK_SET_ERR( rval, "Error getting 3d elements" );

    PUSH_TIMER()
    // Create a tree to use for the location service
    // Can we accelerate this setup phase on a GPU as well ??
    // Or may be use Kokkos/OpenMP for CPU executor ?
    AdaptiveKDTree kd( &mb );
    EntityHandle tree_root = 0;
    rval                   = kd.build_tree( elems, &tree_root );MB_CHK_ERR( rval );
    //AdaptiveKDTree tree(&mb, elems, &fileset);

    // Build the SpatialLocator
    //SpatialLocator sl( &mb, elems, &tree );
    POP_TIMER( "KdTree-Setup" )
    PRINT_TIMER( "KdTree-Setup" )

    // Get the box extents
    CartVect min, max, box_extents, pos;
    unsigned  depth;
    kd.get_info( tree_root, &min[0], &max[0], depth );

    cout << "box: " << min << " " <<  max << "depth: " << depth << "\n";
    BoundBox box (min,max);
    box_extents  = 1.1 * (max - min);

    // Use Kokkos views to transfer r-o data from CPU to GPU
    // KdTreeView treeview (&tree, KOKKOS_STUFF);

    // Query at random places in the tree
    /*CartVect params;
    int is_inside  = 0;
    int num_inside = 0;
    EntityHandle elem; */

    PUSH_TIMER()

    EntityHandle leaf_out;
    for( int i = 0; i < num_queries; i++ )
    {
        pos  = box.bMin + CartVect( box_extents[0] * .01 * ( rand() % 100 ), box_extents[1] * .01 * ( rand() % 100 ),
                                   box_extents[2] * .01 * ( rand() % 100 ) );
        // Query technically on a GPU datastructure here.
        rval = kd.point_search( &pos[0], leaf_out);  MB_CHK_ERR( rval );

    }
    POP_TIMER( "KdTree-Query" )
    PRINT_TIMER( "KdTree-Query" )



    return 0;
}

