/*
 * This imoab_read_map test will simulate coupling between 2 components
 * 2 meshes will be loaded from 2 files (src, tgt), and one map file
 * after the map is read, in parallel, on coupler pes, with distributed rows, the
 * coupler meshes for source and target will be generated, in a migration step,
 * in which we will migrate from target pes according to row ids, to coupler target mesh,
 *  and from source to coverage mesh mesh on coupler. During this migration, par comm graphs
 *  will be established between source and coupler and target and coupler, which will assist
 *  in field transfer from source to target, through coupler
 *
 */

#include "moab/Core.hpp"
#ifndef MOAB_HAVE_MPI
#error mbtempest tool requires MPI configuration
#endif

// MPI includes
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"

#include "moab/iMOAB.h"
#include "TestUtil.hpp"
#include "moab/CpuTimer.hpp"
#include "moab/ProgOptions.hpp"
#include <iostream>
#include <sstream>

#include "imoab_coupler_utils.hpp"

#ifndef MOAB_HAVE_TEMPESTREMAP
#error The climate coupler test example requires MOAB configuration with TempestRemap
#endif

int main( int argc, char* argv[] )
{
    int ierr;
    int rankInGlobalComm, numProcesses;
    MPI_Group jgroup;
    std::string readopts( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS" );

    // Timer data
    moab::CpuTimer timer;
    double timer_ops;
    std::string opName;

    int repartitioner_scheme = 0;
#ifdef MOAB_HAVE_ZOLTAN
    repartitioner_scheme = 2;  // use the graph partitioner in that caseS
#endif

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rankInGlobalComm );
    MPI_Comm_size( MPI_COMM_WORLD, &numProcesses );

    MPI_Comm_group( MPI_COMM_WORLD, &jgroup );  // all processes in jgroup

    std::string atmFilename = TestDir + "/wholeATM_T.h5m";
    // on a regular case,  5 ATM, 6 CPLATM (ATMX), 17 OCN     , 18 CPLOCN (OCNX)  ;
    // intx atm/ocn is not in e3sm yet, give a number
    //   6 * 100+ 18 = 618 : atmocnid
    // 9 LND, 10 CPLLND
    //   6 * 100 + 10 = 610  atmlndid:
    // cmpatm is for atm on atm pes
    // cmpocn is for ocean, on ocean pe
    // cplatm is for atm on coupler pes
    // cplocn is for ocean on coupelr pes
    // atmocnid is for intx atm / ocn on coupler pes
    //
    int rankInAtmComm = -1;
    int cmpatm        = 5,
        cplatm        = 6;  // component ids are unique over all pes, and established in advance;

    std::string ocnFilename = TestDir + "/recMeshOcn.h5m";
    std::string mapFilename = TestDir + "/atm_ocn_map.nc"; //this is a netcdf file!

    std::string baseline    = TestDir + "/baseline1.txt";
    int rankInOcnComm       = -1;
    int cmpocn = 17, cplocn = 18,
        atmocnid = 618;  // component ids are unique over all pes, and established in advance;

    int rankInCouComm = -1;

    int nghlay = 0;  // number of ghost layers for loading the file
    std::vector< int > groupTasks;
    int startG1 = 0, startG2 = 0, endG1 = numProcesses - 1, endG2 = numProcesses - 1;

    int startG4 = startG1, endG4 = endG1;  // these are for coupler layout
    int context_id = -1;                   // used now for freeing buffers

    // default: load atm on 2 proc, ocean on 2,
    // load map on 2 also, in parallel, distributed by rows (which is very bad actually for ocean mesh, because
    // probably all source cells will be involved in coverage mesh on both tasks

    ProgOptions opts;
    opts.addOpt< std::string >( "atmosphere,t", "atm mesh filename (source)", &atmFilename );
    opts.addOpt< std::string >( "ocean,m", "ocean mesh filename (target)", &ocnFilename );
    opts.addOpt< std::string >( "map_file,w", "map file from source to target", &mapFilename );

    opts.addOpt< int >( "startAtm,a", "start task for atmosphere layout", &startG1 );
    opts.addOpt< int >( "endAtm,b", "end task for atmosphere layout", &endG1 );

    opts.addOpt< int >( "startOcn,c", "start task for ocean layout", &startG2 );
    opts.addOpt< int >( "endOcn,d", "end task for ocean layout", &endG2 );

    opts.addOpt< int >( "startCoupler,g", "start task for coupler layout", &startG4 );
    opts.addOpt< int >( "endCoupler,j", "end task for coupler layout", &endG4 );

    opts.addOpt< int >( "partitioning,p", "partitioning option for migration", &repartitioner_scheme );

    int n = 1;  // number of send/receive / project / send back cycles
    opts.addOpt< int >( "iterations,n", "number of iterations for coupler", &n );

    bool no_regression_test = false;
    opts.addOpt< void >( "no_regression,r", "do not do regression test against baseline 1", &no_regression_test );
    opts.parseCommandLine( argc, argv );

    char fileWriteOptions[] = "PARALLEL=WRITE_PART";

    if( !rankInGlobalComm )
    {
        std::cout << " atm file: " << atmFilename << "\n   on tasks : " << startG1 << ":" << endG1 <<
            "\n ocn file: " << ocnFilename << "\n     on tasks : " << startG2 << ":" << endG2 <<
            "\n map file:" << mapFilename << "\n     on tasks : " << startG4 << ":" << endG4 << "\n";
    }

    return 0;
}

