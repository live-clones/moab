#include "moab/ProgOptions.hpp"
#include "moab/Core.hpp"
#include <iostream>

using namespace moab;
using namespace std;

int main( int argc, char* argv[] )
{

    ProgOptions opts;

    std::string inputfile, outfile( "out.h5m" ), physgridfile, variable_name;

    opts.addOpt< std::string >( "input,i", "input mesh filename", &inputfile );
    opts.addOpt< std::string >( "output,o", "output mesh filename", &outfile );
    opts.addOpt< std::string >( "phys,p", "phys grid solution filename", &physgridfile );
    opts.addOpt< std::string >( "var,v", "variable to extract and add to output file", &variable_name );

    opts.parseCommandLine( argc, argv );

    ErrorCode rval;
    Core* mb = new Core();

    EntityHandle fset1, fset2;
    rval = mb->create_meshset(MESHSET_SET, fset1);MB_CHK_SET_ERR( rval, "can't create mesh set" );
    rval = mb->load_file( inputfile.c_str() , &fset1);MB_CHK_SET_ERR( rval, "can't load input file" );


    cout << " opened " << inputfile << " with initial h5m data.\n";

    rval = mb->create_meshset(MESHSET_SET, fset2); MB_CHK_SET_ERR( rval, "can't create mesh set" );
    rval = mb->load_file( physgridfile.c_str() , &fset2);MB_CHK_SET_ERR( rval, "can't load phys grid file" );

    Tag tagv;
    rval = mb->tag_get_handle(variable_name.c_str(), tagv); MB_CHK_SET_ERR( rval, "can't get tag handle" );

    Tag gitag = mb->globalId_tag();

    Range verts; // from phys grid
    rval = mb->get_entities_by_dimension(fset2, 0, verts); MB_CHK_SET_ERR( rval, "can't get vertices" );
    std::vector<int> gids;
    gids.resize (verts.size());
    rval = mb->tag_get_data(gitag, verts, &gids[0]); MB_CHK_SET_ERR( rval, "can't get gi tag values" );
    std::vector<double> valsTag;
    valsTag.resize(verts.size());
    rval = mb->tag_get_data(tagv, verts, &valsTag[0]);MB_CHK_SET_ERR( rval, "can't get tag vals" );
    Range cells;

    rval = mb->get_entities_by_dimension(fset1, 2, cells); MB_CHK_SET_ERR( rval, "can't get cells");

    std::map<int, double> valsByID;
    for (int i=0; i<(int)gids.size(); i++)
        valsByID[gids[i]] = valsTag[i];

    // set now cells values
    std::vector<int> cellsIds;
    cellsIds.resize(cells.size());
    rval = mb->tag_get_data(gitag, cells, &cellsIds[0]); MB_CHK_SET_ERR( rval, "can't get cells ids");
    for (int i=0; i<(int)cells.size(); i++)
    {
        valsTag[i] = valsByID[cellsIds[i]];
    }
    rval = mb->tag_set_data(tagv, cells, &valsTag[0] ); MB_CHK_SET_ERR( rval, "can't set  cells tags");

    rval = mb->write_file(outfile.c_str(), 0, 0, &fset1, 1); MB_CHK_SET_ERR( rval, "can't write file");

    return 0;

}
    // open the netcdf file, and see if it has that variable we are looking for
