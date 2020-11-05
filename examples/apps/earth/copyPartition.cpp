/** @example addPCdata.cpp  Add point cloud data
 * this tool will take an existing h5m fine atm mesh file and add data from an h5m  type file with
 * point cloud mesh will support mainly showing the data with Visit
 *
 * example of usage:
 * ./mbaddpcdata -i wholeFineATM.h5m -s wholeLND_proj01.h5m -o atm_a2l.h5m -v a2lTbot_proj
 *
 * it should work also for pentagon file style data
 * ./mbaddpcdata -i <MOABsource>/MeshFiles/unittest/penta3d.h5m -s wholeLND_proj01.h5m -o
 * atm_a2l.h5m -v a2lTbot_proj -p 1
 *
 * Basically, will output a new h5m file (atm_a2l.h5m), which has an extra tag, corresponding to the
 * variable a2lTbot_proj from the file wholeLND_proj01.h5m; matching is based on the global ids
 * between what we think is the order on the original file (wholeFineATM.h5m) and the order of
 * surfdata_ne11np4_simyr1850_c160614.nc
 *
 *  file  wholeFineATM.h5m is obtained from a coupled run in e3sm, with the ne 11, np 4,
 *
 */
#include "moab/ProgOptions.hpp"
#include "moab/Core.hpp"

#include <math.h>
#include <sstream>

using namespace moab;

int main( int argc, char* argv[] )
{

    ProgOptions opts;

    std::string physfile, pg2file, outfile;

    opts.addOpt< std::string >( "physgridfile,p", "phys grid filename", &physfile );
    opts.addOpt< std::string >( "pg2file,g", "pg2 mesh file", &pg2file );
    opts.addOpt< std::string >( "output,o", "output mesh filename", &outfile );


    opts.parseCommandLine( argc, argv );

    std::cout << "phys grid cloud file: " << physfile << "\n";
    std::cout << "pg2 mesh file: " << pg2file << "\n";
    std::cout << "output file: " << outfile << "\n";

    if( physfile.empty() )
    {
        opts.printHelp();
        return 0;
    }
    ErrorCode rval;
    Core* mb = new Core();

    rval = mb->load_file( physfile.c_str() );MB_CHK_SET_ERR( rval, "can't load phys grid file" );

    Core* mb2 = new Core();
    rval      = mb2->load_file( pg2file.c_str() );MB_CHK_SET_ERR( rval, "can't load pg2 mesh file" );

    Tag globalIDTag1 = mb->globalId_tag();
    Tag parti;
    rval = mb->tag_get_handle("partition", parti);MB_CHK_SET_ERR( rval, "can't get partition tag phys grid mesh " );

    Tag globalIDTag2 = mb2->globalId_tag();

    Range verts1;
    rval = mb->get_entities_by_dimension(0, 0, verts1);MB_CHK_SET_ERR( rval, "can't get vertices " );

    std::vector<int>  partValues;
    partValues.resize(verts1.size());
    rval = mb->tag_get_data(parti, verts1, &partValues[0]);MB_CHK_SET_ERR( rval, "can't get parts values on vertices " );

    Range cells;
    rval = mb2->get_entities_by_dimension(0, 2, cells);
    std::vector<int> globalIdsCells;
    globalIdsCells.resize(cells.size());
    rval = mb2->tag_get_data(globalIDTag2, cells, &globalIdsCells[0]); MB_CHK_SET_ERR( rval, "can't get global ids cells " );

    std::vector<int> globalIdsVerts;
    globalIdsVerts.resize(verts1.size());
    rval = mb->tag_get_data(globalIDTag1, verts1, &globalIdsVerts[0]); MB_CHK_SET_ERR( rval, "can't get global ids cells " );

    Tag partTag;
    rval = mb2->tag_get_handle("PARALLEL_PARTITION", partTag);MB_CHK_SET_ERR( rval, "can't partition tag " );

    Range sets;
    rval = mb2->get_entities_by_type_and_tag( 0, MBENTITYSET, &partTag, NULL, 1, sets );MB_CHK_ERR( rval );

    std::vector<int> setValues;
    setValues.resize(sets.size());
    rval = mb2->tag_get_data(partTag, sets, &setValues[0]);MB_CHK_ERR( rval );

    std::map<int, EntityHandle>  valToSet;
    int i=0;
    for (Range::iterator st=sets.begin(); st!=sets.end(); st++, i++)
    {
        valToSet[setValues[i]] = *st;
    }
    // now, every cell will be put into one set, by looking at the global id of cell

    std::map<int, EntityHandle> gidToCell;
    i=0;
    for (Range::iterator it = cells.begin(); it!=cells.end(); it++, i++)
    {
        gidToCell[globalIdsCells[i]] = *it;
    }
    // empty all sets
    rval = mb2->clear_meshset(sets);MB_CHK_ERR( rval );

    // look now at parti values for vertices, and their global ids
    for (i=0; i<(int)verts1.size(); i++)
    {
        int part = partValues[i];
        int gid = globalIdsVerts[i];
        EntityHandle set1 = valToSet[part];
        EntityHandle cell = gidToCell[gid];
        rval = mb2->add_entities(set1, &cell, 1); MB_CHK_ERR( rval );
    }

    rval = mb2->write_file( outfile.c_str() );MB_CHK_SET_ERR( rval, "can't write file" );

    delete mb;
    delete mb2;

    return 0;
}
