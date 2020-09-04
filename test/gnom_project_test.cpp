/*
 * gnom_project_test.cpp
 * will test new global gnomonic projection method, to be used by zoltan for partitioning
 *
 */
#include <iostream>
#include <sstream>

#include "moab/Core.hpp"
#include "moab/Interface.hpp"

#include "moab/IntxMesh/IntxUtils.hpp"
#include "moab/ProgOptions.hpp"
#include "TestUtil.hpp"

using namespace moab;
using namespace std;

int main( int  argc, char*  argv[] )
{
    string filein = STRINGIFY( MESHDIR ) "/mbcslam/eulerHomme.vtk";
    string fileout = "project.vtk";

    ProgOptions opts;
    opts.addOpt< std::string >( "model,m", "input file ", &filein );

    opts.addOpt< std::string >( "output,o", "output filename", &fileout );

    opts.parseCommandLine( argc, argv );



    Core moab;
    Interface* mb = &moab;
    EntityHandle sf;
    ErrorCode rval = mb->create_meshset( MESHSET_SET, sf );MB_CHK_ERR( rval );

    rval = mb->load_file( filein.c_str(), &sf );MB_CHK_ERR( rval );

    double R = 1.;  // should be input
    EntityHandle outSet;
    rval = mb->create_meshset( MESHSET_SET, outSet );MB_CHK_ERR( rval );

    rval = IntxUtils::global_gnomonic_projection( mb, sf,  R, false, outSet );MB_CHK_ERR( rval );

    rval = mb->write_file(fileout.c_str(), 0, 0, &outSet, 1); ;MB_CHK_ERR( rval );

    return 0;

}
