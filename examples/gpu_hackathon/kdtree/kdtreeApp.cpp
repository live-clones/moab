/** \brief This test shows how to perform local point-in-element searches with MOAB's new tree
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

#include <iostream>
#include <cstdlib>
#include <cstdio>

#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/CN.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/NestedRefine.hpp"

#include "moab/AdaptiveKDTree.hpp"
#include "moab/SpatialLocator.hpp"

using namespace moab;
using namespace std;

int main( int argc, char** argv )
{
    std::string test_file_name  = string( MESH_DIR ) + string( "/64bricks_512hex_256part.h5m" );
    int num_queries             = 10000;
    int dimension               = 3;
    int uniformRefinementLevels = 0;

    ProgOptions opts;
    // Need option handling here for input filename
    opts.addOpt< int >( "dim,d", "Dimension of the problem and mesh (default=1)", &dimension );
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
        assert( level_sets.size() == uniformRefinementLevels + 1 );
        fileset = level_sets[uniformRefinementLevels];
    }
    else
        fileset = baseFileset;

    // Get all 3d elements in the file
    Range elems;
    rval = mb.get_entities_by_dimension( fileset, dimension, elems );MB_CHK_SET_ERR( rval, "Error getting 3d elements" );

    // Create a tree to use for the location service
    AdaptiveKDTree tree(&mb, elems, &fileset);

    // Build the SpatialLocator
    SpatialLocator sl( &mb, elems, &tree );

    // Get the box extents
    CartVect box_extents, pos;
    BoundBox box = sl.local_box();
    box_extents  = 1.1 * (box.bMax - box.bMin);

    // Query at random places in the tree
    CartVect params;
    int is_inside  = 0;
    int num_inside = 0;
    EntityHandle elem;
    for( int i = 0; i < num_queries; i++ )
    {
        pos  = box.bMin + CartVect( box_extents[0] * .01 * ( rand() % 100 ), box_extents[1] * .01 * ( rand() % 100 ),
                                   box_extents[2] * .01 * ( rand() % 100 ) );
        rval = sl.locate_point( pos.array(), elem, params.array(), &is_inside, 0.0, 0.0 );MB_CHK_ERR( rval );
        if( is_inside ) num_inside++;
    }

    cout << "Mesh contains " << elems.size() << " elements of type "
         << CN::EntityTypeName( mb.type_from_handle( *elems.begin() ) ) << endl;
    cout << "Bounding box min-max = (" << box.bMin[0] << "," << box.bMin[1] << "," << box.bMin[2] << ")-("
         << box.bMax[0] << "," << box.bMax[1] << "," << box.bMax[2] << ")" << endl;
    cout << "Queries inside box = " << num_inside << "/" << num_queries << " = "
         << 100.0 * ( (double)num_inside ) / num_queries << "%" << endl;

    return 0;
}
