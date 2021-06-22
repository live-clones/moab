/** @example CoalescePartitions.cpp
 * Description: read a mesh that has par partitions and coalesce by id \n
 *
 * To run: ./CoalescePartitions [input] [outfile]  \n
 */

#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"
#include <iostream>

using namespace moab;
using namespace std;

#ifndef MESH_DIR
#define MESH_DIR "."
#endif

// at most 20 edges per polygon
#define MAXEDGES 20
// Note: change the file name below to test a trivial "No such file or directory" error
string test_file_name = string( MESH_DIR ) + string( "/io/out3m.h5m" );
string output         = string( "mesh3p.h5m" );

int main( int argc, char** argv )
{
    // Get MOAB instance
    Interface* mb = new( std::nothrow ) Core;
    if( NULL == mb ) return 1;

    // Need option handling here for input filename
    if( argc > 1 )
    {
        // User has input a mesh file
        test_file_name = argv[1];
    }
    if( argc > 2 ) output = argv[2];
    int dimension = 3;
    if( argc > 3 ) dimension = atoi(argv[3]);

    std::cout << "Run: " << argv[0] << " " << test_file_name << " " << output  << " " << dimension << "\n";
    // Load the mesh with multiple parts
    ErrorCode rval = mb->load_mesh( test_file_name.c_str() );MB_CHK_ERR( rval );

    // Get verts entities, by type
    Range verts;
    rval = mb->get_entities_by_type( 0, MBVERTEX, verts );MB_CHK_ERR( rval );

    // this is our convention
    Tag part_set_tag;
    int dum_id = -1;
    rval     = mb->tag_get_handle( "PARALLEL_PARTITION", 1, MB_TYPE_INTEGER, part_set_tag,
                                     MB_TAG_SPARSE | MB_TAG_CREAT, &dum_id );MB_CHK_ERR( rval );

    // get any sets already with this tag, and clear them
    Range tagged_sets;
    rval =  mb->get_entities_by_type_and_tag( 0, MBENTITYSET, &part_set_tag, NULL, 1, tagged_sets, Interface::UNION );MB_CHK_ERR( rval );

    std::vector<int> vals;
    vals.resize(tagged_sets.size());
    rval = mb->tag_get_data(part_set_tag, tagged_sets, &vals[0]); MB_CHK_ERR( rval );
    // coalesce sets by vals

    std::map<int, EntityHandle> mset;
    for (int i=0; i<(int)tagged_sets.size(); i++)
    {
        int val = vals[i];
        EntityHandle eh = tagged_sets[i];
        auto it = mset.find(val);
        if (it != mset.end())
        {
            // move the eh set to *it
            EntityHandle existingSet = it->second;
            Range elems;
            rval = mb->get_entities_by_handle(eh, elems); MB_CHK_ERR( rval );
            rval = mb->add_entities(existingSet, elems);  MB_CHK_ERR( rval );
            // remove this set
            rval = mb->delete_entities(&eh, 1); MB_CHK_ERR( rval );
        }
    }
    rval = mb->write_file( output.c_str() );MB_CHK_ERR( rval );

    delete mb;

    return 0;
}
