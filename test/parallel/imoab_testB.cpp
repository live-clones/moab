/**
 * imoab_testB.cpp
 *
 *  This imoab_testB test will simulate coupling between 2 components
 *  meshes will be loaded from 2 files (ocean and atm pg2)
 *  ocn will be migrated to coupler pes and compute intx between atm pg2 and ocn, with atm
 *  as target; atm will be already on coupler PEs; Or should we load on atm pes, and "migrate"
 *  to coupler pes, even if they are the same?
 *  We will compute projection weights too, ocn to atmpg2(FV-FV)
 *   and some data will move from ocn to atm, at every "time" step
 *
 * first, intersect ocn and atm pg2
 * maybe later will identify some equivalent to bilinear maps
 */

#include "moab/Core.hpp"
#ifndef MOAB_HAVE_MPI
#error This test requires MPI configuration
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

using namespace moab;

//#define GRAPH_INFO

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

    int rankInAtmComm = -1;
    int cmpatm = 5, cplatm = 6;  //
    std::string atmFilename =
            "../../sandbox/MeshFiles/e3sm/o_ne11pg2/ne11pg2_inf.h5m";
    // get the ocn from the ol_ne4pg2 projection folder, it has some data
    int rankInOcnComm = -1;
    int cmpocn = 17, cplocn = 18, ocnatmid = 1806;  // 18*100+6
    std::string ocnFilename =
            "../../sandbox/MeshFiles/e3sm/ol_ne4pg2/ocn.h5m";

    // load atm directly on coupler pes! and do intx there, after ocean migration
    // ocean is source !

    int rankInCouComm = -1;

    int nghlay = 0;  // number of ghost layers for loading the file
    std::vector< int > groupTasks;
    // group 1 is atm, 2 is ocn; atm not used ? skip group 3 (land)
    int startG1 = 0, startG2 = 0, endG1 = numProcesses - 1, endG2 = numProcesses - 1;
    int startG4 = startG1, endG4 = endG1;  // these are for coupler layout
    int context_id = -1;                   // used now for freeing buffers

    // default: load atm on 2 proc, ocean on 2, land on 2; migrate to 2 procs, then compute intx
    // later, we need to compute weight matrix with tempestremap

    ProgOptions opts;
    opts.addOpt< std::string >( "atmosphere,t", "atm mesh filename (target)", &atmFilename );

    opts.addOpt< std::string >( "ocean,m", "ocean mesh filename (source)", &ocnFilename );

    opts.addOpt< int >( "startAtm,a", "start task for atmosphere layout", &startG1 );
    opts.addOpt< int >( "endAtm,b", "end task for atmosphere layout", &endG1 );

    opts.addOpt< int >( "startOcn,c", "start task for ocean layout", &startG2 );
    opts.addOpt< int >( "endOcn,d", "end task for ocean layout", &endG2 );



    opts.addOpt< int >( "startCoupler,g", "start task for coupler layout", &startG4 );
    opts.addOpt< int >( "endCoupler,j", "end task for coupler layout", &endG4 );

    opts.addOpt< int >( "partitioning,p", "partitioning option for migration", &repartitioner_scheme );

    int n = 1;  // number of send/receive / project / send back cycles
    opts.addOpt< int >( "iterations,n", "number of iterations for coupler", &n );

    opts.parseCommandLine( argc, argv );

    char fileWriteOptions[] = "PARALLEL=WRITE_PART";

    if( !rankInGlobalComm )
    {
        // do not use group 1!
        std::cout << " atm file: " << atmFilename << "\n   on tasks : " << startG4 << ":" << endG4 <<
           "\n ocn file: " << ocnFilename << "\n     on tasks : " << startG2 << ":" << endG2 <<
           "\n  partitioning (0 trivial, 1 graph, 2 geometry, 3 gnomonic :"
                 << repartitioner_scheme << "\n  " << " number of iterations: " << n << "\n";
    }
    MPI_Group ocnPEGroup;
    MPI_Comm ocnComm;
    ierr = create_group_and_comm( startG2, endG2, jgroup, &ocnPEGroup, &ocnComm );
    CHECKIERR( ierr, "Cannot create ocn MPI group and communicator " )

    // we will always have a coupler
    MPI_Group couPEGroup;
    MPI_Comm couComm;
    ierr = create_group_and_comm( startG4, endG4, jgroup, &couPEGroup, &couComm );
    CHECKIERR( ierr, "Cannot create cpl MPI group and communicator " )

    // ocn_coupler
    MPI_Group joinOcnCouGroup;
    MPI_Comm ocnCouComm;
    ierr = create_joint_comm_group( ocnPEGroup, couPEGroup, &joinOcnCouGroup, &ocnCouComm );
    CHECKIERR( ierr, "Cannot create joint ocn cou communicator" )

    ierr = iMOAB_Initialize( argc, argv );  // not really needed anything from argc, argv, yet; maybe we should
    CHECKIERR( ierr, "Cannot initialize iMOAB" )

    int cmpAtmAppID       = -1;
    iMOAB_AppID cmpAtmPID = &cmpAtmAppID;  // atm
    int cplAtmAppID       = -1;            // -1 means it is not initialized
    iMOAB_AppID cplAtmPID = &cplAtmAppID;  // atm on coupler PEs

    int cmpOcnAppID       = -1;
    iMOAB_AppID cmpOcnPID = &cmpOcnAppID;        // ocn
    int cplOcnAppID = -1, cplOcnAtmAppID = -1;   // -1 means it is not initialized
    iMOAB_AppID cplOcnPID    = &cplOcnAppID;     // ocn on coupler PEs
    iMOAB_AppID cplOcnAtmPID = &cplOcnAtmAppID;  // intx ocn - atm on coupler PEs

    if( couComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( couComm, &rankInCouComm );
        // Register all the applications on the coupler PEs
        ierr = iMOAB_RegisterApplication( "ATMX", &couComm, &cplatm, cplAtmPID );  // atm on coupler pes
        CHECKIERR( ierr, "Cannot register ATM over coupler PEs" )

        ierr = iMOAB_RegisterApplication( "OCNX", &couComm, &cplocn, cplOcnPID );  // ocn on coupler pes
        CHECKIERR( ierr, "Cannot register OCN over coupler PEs" )

    }
    if( ocnComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( ocnComm, &rankInOcnComm );
        ierr = iMOAB_RegisterApplication( "OCN1", &ocnComm, &cmpocn, cmpOcnPID );
        CHECKIERR( ierr, "Cannot register OCN App" )
    }

    // load atm on coupler pes directly
    if( couComm != MPI_COMM_NULL )
    {
        // load atm mesh on copler pes
        ierr = iMOAB_LoadMesh( cplAtmPID, atmFilename.c_str(), readopts.c_str(), &nghlay, atmFilename.length(),
                               readopts.length() );
        CHECKIERR( ierr, "Cannot load atm mesh on coupler pes" )

    }

    // load ocn and migrate to coupler
    repartitioner_scheme = 2; // will use RCB
    // ocean
    ierr = setup_component_coupler_meshes( cmpOcnPID, cmpocn, cplOcnPID, cplocn, &ocnComm, &ocnPEGroup, &couComm,
                                            &couPEGroup, &ocnCouComm, ocnFilename, readopts, nghlay, repartitioner_scheme );
    CHECKIERR( ierr, "Cannot load and migrate ocn mesh" )
    if( couComm != MPI_COMM_NULL )
    {
        char outputFileTgt5[] = "recvOcn5.h5m";
        PUSH_TIMER( "Write migrated OCN mesh on coupler PEs" )
        ierr = iMOAB_WriteMesh( cplOcnPID, outputFileTgt3, fileWriteOptions, strlen( outputFileTgt3 ),
                                strlen( fileWriteOptions ) );
        CHECKIERR( ierr, "cannot write ocn mesh after receiving" )
        POP_TIMER( couComm, rankInCouComm )
    }


    return 0;
}
