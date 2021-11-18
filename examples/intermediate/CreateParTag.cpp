/** \brief  create parallel partition sets from a tag
 *  create a PARALLEL_PARTITION sets with value the same as a partition tag 
    it is output from a e3sm land run
 */

#include <iostream>
#include <cstdlib>
#include <cstdio>

#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/ProgOptions.hpp"
#include <set>

using namespace moab;
using namespace std;

#ifdef MOAB_HAVE_HDF5
string test_file_name = string( MESH_DIR ) + string( "/64bricks_512hex_256part.h5m" );
#endif
int main( int argc, char** argv )
{
    ProgOptions opts;

    string inputFile = test_file_name;
    opts.addOpt< string >( "inFile,i", "Specify the input file name string ", &inputFile );

    string outFile = "out.h5m";
    opts.addOpt< string >( "outFile,o", "Specify the output file name string ", &outFile );

    opts.parseCommandLine( argc, argv );

    // Instantiate
    Core mb;

    ErrorCode rval = mb.load_file( inputFile.c_str() );MB_CHK_SET_ERR( rval, "Error loading file" );

    cout << " reading file " << inputFile << "\n";

    // Get all 2d elements in the file set
    Range elems;
    rval = mb.get_entities_by_dimension( 0, 2, elems );MB_CHK_SET_ERR( rval, "Error getting 2d elements" );

    cout << "number of cells: " << elems.size() << "\n";

    std::map<int, std::set<EntityHandle> > cellsMap;
    std::vector<int> parti;
    parti.resize(elems.size());
    Tag partitionTag;
    rval = mb.tag_get_handle("partition", partitionTag); MB_CHK_SET_ERR( rval, "can't get partition tag" );
    rval = mb.tag_get_data(partitionTag, elems, &parti[0]);  MB_CHK_SET_ERR( rval, "can't get partition tag values" );

    int i=0;
    for (Range::iterator it=elems.begin(); it!= elems.end(); it++, i++)
    {
        EntityHandle cell = *it;
	    cellsMap[ parti[i] ].insert(cell);
    }

    Tag ptag;
    int dum_id       = -1;
    rval = mb.tag_get_handle( "PARALLEL_PARTITION", 1, MB_TYPE_INTEGER, ptag,
                                                       MB_TAG_SPARSE | MB_TAG_CREAT, &dum_id );
    for (auto it = cellsMap.begin(); it!= cellsMap.end(); it++)
    {
        int val = it->first;
	    std::set<EntityHandle> & setCells = it->second;
	    if (setCells.size() > 0)
	    {
            std::vector<EntityHandle> setRange ;
            std::copy( setCells.begin(), setCells.end(), std::back_inserter( setRange ) );
            EntityHandle pset;
            rval = mb.create_meshset(MESHSET_SET, pset); MB_CHK_ERR( rval);
            rval = mb.add_entities(pset, &setRange[0], setRange.size() ); MB_CHK_ERR( rval);
            rval = mb.tag_set_data(ptag, &pset, 1, &val);  MB_CHK_ERR( rval);
	    }
    }


    rval = mb.write_file( outFile.c_str() );MB_CHK_SET_ERR( rval, "Error writing file" );

    return 0;
}
