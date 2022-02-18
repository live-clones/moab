/**
 * imoab_testE.cpp
 *
 *  two atmosphere meshes:  land RCB with atmosphere inferred from land, and ocean RCB with atmosphere
 *   inferred from ocean); will not require new code, but requires a new driver
 *  atm1 will be migrated to coupler pes with inference from land RCB
 *  atm2 will be migrated to coupler pes with inference from ocn RCB
 *
 *  first part of test E: just infer atm from land RCB test
 *
 *  the second part: atm2 inferred from ocean RCB
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
//#define  ENABLE_OCNATM_COUPLING
//#define VERBOSE

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

    // this test runs only with zoltan anyway; do not need to check zoltan build
    int repartitioner_scheme = 2;

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rankInGlobalComm );
    MPI_Comm_size( MPI_COMM_WORLD, &numProcesses );

    MPI_Comm_group( MPI_COMM_WORLD, &jgroup );  // all processes in jgroup

    int rankInAtmComm = -1;
    int cmpatm = 5, cplatm = 6;  // cplatm for land context is 6
    int cplatm2 = 7; // cpl atm 2 for ocean context
    std::string atmFilename =
            "../../sandbox/MeshFiles/e3sm/o_ne11pg2/ne11pg2_inf.h5m";
    // get the ocn from the ol_ne4pg2 projection folder, it has some data
#ifdef ENABLE_OCNATM_COUPLING
    int rankInOcnComm = -1;
    int cmpocn = 17, cplocn = 18, ocnatmid = 1806;  // 18*100+6
    std::string ocnFilename =
            "../../sandbox/MeshFiles/e3sm/ol_ne4pg2/ocn.h5m";
#endif
    std::string lndFilename =
        "../../sandbox/MeshFiles/e3sm/ol_ne4pg2/land.h5m";  // need to add a land result file from land PEs
    int rankInLndComm = -1;
    int cpllnd = 10, cmplnd = 9, lndatmid = 1006;  // component ids are unique over all pes, and established in advance;

    int rankInCouComm = -1;

    int nghlay = 0;  // number of ghost layers for loading the file
    std::vector< int > groupTasks;
    // group 1 is atm, 2 is ocn;
    int startG1 = 0, endG1 = numProcesses - 1;
#ifdef ENABLE_OCNATM_COUPLING
    int startG2 = 0, endG2 = numProcesses - 1;
#endif
    int startG3 = startG1, endG3 = endG1;  // these are for land layout
    int startG4 = startG1, endG4 = endG1;  // these are for coupler layout
    int context_id = -1;                   // used now for freeing buffers

    // default: load atm on 2 proc, ocean on 2, land on 2; migrate to 2 procs, then compute intx
    // later, we need to compute weight matrix with tempestremap

    ProgOptions opts;
    opts.addOpt< std::string >( "atmosphere,t", "atm mesh filename (target)", &atmFilename );
    opts.addOpt< int >( "startAtm,a", "start task for atmosphere layout", &startG1 );
    opts.addOpt< int >( "endAtm,b", "end task for atmosphere layout", &endG1 );

#ifdef ENABLE_OCNATM_COUPLING
    opts.addOpt< std::string >( "ocean,m", "ocean mesh filename (source)", &ocnFilename );

    opts.addOpt< int >( "startOcn,c", "start task for ocean layout", &startG2 );
    opts.addOpt< int >( "endOcn,d", "end task for ocean layout", &endG2 );
#endif

    opts.addOpt< std::string >( "land,l", "land mesh filename (source)", &lndFilename );

    opts.addOpt< int >( "startLnd,e", "start task for land layout", &startG3 );
    opts.addOpt< int >( "endLnd,f", "end task for land layout", &endG3 );

    opts.addOpt< int >( "startCoupler,g", "start task for coupler layout", &startG4 );
    opts.addOpt< int >( "endCoupler,j", "end task for coupler layout", &endG4 );

    opts.addOpt< int >( "partitioning,p", "partitioning option for migration", &repartitioner_scheme );

    int n = 1;  // number of send/receive / project / send back cycles
    opts.addOpt< int >( "iterations,n", "number of iterations for coupler", &n );

    bool analytic_field = false;

    opts.addOpt< void >( "analytic,q", "analytic field", &analytic_field );

    bool skip_writing_map = false;
    opts.addOpt< void >( "skip_writing,s", "do not write map file ", &skip_writing_map );

    opts.parseCommandLine( argc, argv );

    char fileWriteOptions[] = "PARALLEL=WRITE_PART";

    if( !rankInGlobalComm )
    {
        // do not use group 1!
        std::cout << " atm file: " << atmFilename << "\n   on tasks : " << startG4 << ":" << endG4 <<
#ifdef ENABLE_OCNATM_COUPLING
           "\n ocn file: " << ocnFilename << "\n     on tasks : " << startG2 << ":" << endG2 <<
#endif
           "\n land file: " << lndFilename << "\n     on tasks : " << startG3 << ":" << endG3 <<
           "\n  partitioning (0 trivial, 1 graph, 2 geometry, 3 gnomonic, 4 gnomonic + store  :"
                 << repartitioner_scheme << "\n  " << " number of iterations: " << n << "\n";
    }
    MPI_Group atmPEGroup;
    MPI_Comm atmComm;
    ierr = create_group_and_comm( startG1, endG1, jgroup, &atmPEGroup, &atmComm );
    CHECKIERR( ierr, "Cannot create atm MPI group and communicator " )

#ifdef ENABLE_OCNATM_COUPLING
    MPI_Group ocnPEGroup;
    MPI_Comm ocnComm;
    ierr = create_group_and_comm( startG2, endG2, jgroup, &ocnPEGroup, &ocnComm );
    CHECKIERR( ierr, "Cannot create ocn MPI group and communicator " )
#endif
    MPI_Group lndPEGroup;
    MPI_Comm lndComm;
    ierr = create_group_and_comm( startG3, endG3, jgroup, &lndPEGroup, &lndComm );
    CHECKIERR( ierr, "Cannot create lnd MPI group and communicator " )

    // we will always have a coupler
    MPI_Group couPEGroup;
    MPI_Comm couComm;
    ierr = create_group_and_comm( startG4, endG4, jgroup, &couPEGroup, &couComm );
    CHECKIERR( ierr, "Cannot create cpl MPI group and communicator " )

#ifdef ENABLE_OCNATM_COUPLING
    // ocn_coupler
    MPI_Group joinOcnCouGroup;
    MPI_Comm ocnCouComm;
    ierr = create_joint_comm_group( ocnPEGroup, couPEGroup, &joinOcnCouGroup, &ocnCouComm );
    CHECKIERR( ierr, "Cannot create joint ocn cou communicator" )
#endif
    ierr = iMOAB_Initialize( argc, argv );  // not really needed anything from argc, argv, yet; maybe we should
    CHECKIERR( ierr, "Cannot initialize iMOAB" )

    int cmpAtmAppID       = -1;
    iMOAB_AppID cmpAtmPID = &cmpAtmAppID;  // atm
    // cpl atm for land context
    int cplAtmAppID       = -1;            // -1 means it is not initialized
    iMOAB_AppID cplAtmPID = &cplAtmAppID;  // atm on coupler PEs
    // cpl atm 2 for ocean context
    int cplAtm2AppID       = -1;            // -1 means it is not initialized
    iMOAB_AppID cplAtm2PID = &cplAtm2AppID;  // atm on coupler PEs, ocean context

#ifdef ENABLE_OCNATM_COUPLING
    int cmpOcnAppID       = -1;
    iMOAB_AppID cmpOcnPID = &cmpOcnAppID;        // ocn
    int cplOcnAppID = -1, cplOcnAtm2AppID = -1;   // -1 means it is not initialized
    iMOAB_AppID cplOcnPID    = &cplOcnAppID;     // ocn on coupler PEs
    iMOAB_AppID cplOcnAtm2PID = &cplOcnAtm2AppID;  // intx ocn - atm on coupler PEs
#endif
    int cmpLndAppID       = -1;
    iMOAB_AppID cmpLndPID = &cmpLndAppID;       // lnd
    int cplLndAppID = -1, cplLndAtmAppID = -1;  // -1 means it is not initialized
    iMOAB_AppID cplLndPID = &cplLndAppID;       // land on coupler PEs

    iMOAB_AppID cplLndAtmPID = &cplLndAtmAppID;  // intx lnd - atm on coupler PEs

    if( couComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( couComm, &rankInCouComm );
        // Register all the applications on the coupler PEs
        ierr = iMOAB_RegisterApplication( "ATMX", &couComm, &cplatm, cplAtmPID );  // atm on coupler pes, land context
        CHECKIERR( ierr, "Cannot register ATM over coupler PEs" )
        ierr = iMOAB_RegisterApplication( "LNDX", &couComm, &cpllnd, cplLndPID );  // lnd on coupler pes
        CHECKIERR( ierr, "Cannot register LND over coupler PEs" )

#ifdef ENABLE_OCNATM_COUPLING
        ierr = iMOAB_RegisterApplication( "ATM2X", &couComm, &cplatm2, cplAtm2PID );  // atm on coupler pes, ocean context
        CHECKIERR( ierr, "Cannot register ATM 2 over coupler PEs" )
        ierr = iMOAB_RegisterApplication( "OCNX", &couComm, &cplocn, cplOcnPID );  // ocn on coupler pes
        CHECKIERR( ierr, "Cannot register OCN over coupler PEs" )
#endif

    }
    if( atmComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( atmComm, &rankInAtmComm );
        ierr = iMOAB_RegisterApplication( "ATM1", &atmComm, &cmpatm, cmpAtmPID );
        CHECKIERR( ierr, "Cannot register ATM App" )
    }
#ifdef ENABLE_OCNATM_COUPLING
    if( ocnComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( ocnComm, &rankInOcnComm );
        ierr = iMOAB_RegisterApplication( "OCN1", &ocnComm, &cmpocn, cmpOcnPID );
        CHECKIERR( ierr, "Cannot register OCN App" )
    }
#endif
    // atm_coupler
    MPI_Group joinAtmCouGroup;
    MPI_Comm atmCouComm;
    ierr = create_joint_comm_group( atmPEGroup, couPEGroup, &joinAtmCouGroup, &atmCouComm );
    CHECKIERR( ierr, "Cannot create joint atm cou communicator" )
    // lnd_coupler
    MPI_Group joinLndCouGroup;
    MPI_Comm lndCouComm;
    ierr = create_joint_comm_group( lndPEGroup, couPEGroup, &joinLndCouGroup, &lndCouComm );
    CHECKIERR( ierr, "Cannot create joint ocn cou communicator" )

    // land
    if( lndComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( lndComm, &rankInLndComm );
        ierr = iMOAB_RegisterApplication( "LND1", &lndComm, &cmplnd, cmpLndPID );
        CHECKIERR( ierr, "Cannot register LND App " )
    }
    // we can reuse the lnd, or use repartitioner_scheme = 4; RCB + store to reuse it for atm
    repartitioner_scheme = 4;
    ierr = setup_component_coupler_meshes( cmpLndPID, cmplnd, cplLndPID, cpllnd, &lndComm, &lndPEGroup, &couComm,
                                           &couPEGroup, &lndCouComm, lndFilename, readopts, nghlay,
                                           repartitioner_scheme );

    if( couComm != MPI_COMM_NULL )
    {  // write only for n==1 case
        char outputFileLnd[] = "recvLnd.h5m";
        ierr                 = iMOAB_WriteMesh( cplLndPID, outputFileLnd, fileWriteOptions );
        CHECKIERR( ierr, "cannot write lnd mesh after receiving" )
    }

    // infer atm from land:
    // atm
    repartitioner_scheme = 5; // reuse the partition stored at previous step
    ierr =
        setup_component_coupler_meshes( cmpAtmPID, cmpatm, cplAtmPID, cplatm, &atmComm, &atmPEGroup, &couComm,
                                        &couPEGroup, &atmCouComm, atmFilename, readopts, nghlay, repartitioner_scheme );
    CHECKIERR( ierr, "Cannot load and migrate atm mesh for land context " )
#ifdef VERBOSE
    if( couComm != MPI_COMM_NULL )
    {
        char outputFileAtmInf[] = "recvAtmInf.h5m";
        PUSH_TIMER( couComm, "Write migrated ATM mesh on coupler PEs, inferred from LND" )
        ierr = iMOAB_WriteMesh( cplAtmPID, outputFileAtmInf, fileWriteOptions );
        CHECKIERR( ierr, "cannot write atm mesh after receiving" )
        POP_TIMER( couComm, rankInCouComm )
    }
#endif

#ifdef VERBOSE
    if(  atmComm != MPI_COMM_NULL  )
    {
        char outputFileAtm[] = "AtmAfterMigr.h5m";
        ierr                 = iMOAB_WriteMesh( cmpAtmPID, outputFileAtm, fileWriteOptions );
        CHECKIERR( ierr, "cannot write AtmWithProj3.h5m" )
    }
#endif
    MPI_Barrier( MPI_COMM_WORLD );
    if( couComm != MPI_COMM_NULL )
    {
        // now compute intersection between LNDx and ATMx on coupler PEs
        ierr = iMOAB_RegisterApplication( "LNDATM", &couComm, &lndatmid, cplLndAtmPID );
        CHECKIERR( ierr, "Cannot register lnd_atm intx over coupler pes " )

        PUSH_TIMER( couComm, "Compute LND-ATM mesh intersection" )
        ierr = iMOAB_ComputeMeshIntersectionOnSphere( cplLndPID, cplAtmPID, cplLndAtmPID );
        CHECKIERR( ierr, "failed to compute land - atm intx for mapping" );
        POP_TIMER( couComm, rankInCouComm )
        if (1 == numProcesses)
        {
            char outputFile[] = "IntxLndAtm.h5m";
            ierr = iMOAB_WriteMesh( cplLndAtmPID, outputFile, fileWriteOptions );
        }
    }
    if( lndCouComm != MPI_COMM_NULL )
    {
        // the new graph will be for sending data from lnd comp to coverage mesh for atm mesh;
        // it involves initial lnd app; cmpLndPID; also migrate lnd mesh on coupler pes, cplLndPID
        // results are in cplLndAtmPID, intx mesh; remapper also has some info about coverage mesh
        // after this, the sending of tags from lnd pes to coupler pes will use the new par comm graph, that has more
        // precise info about what to send (specifically for atm cover); every time,
        /// we will use the element global id, which should uniquely identify the element
        PUSH_TIMER( lndCouComm, "Compute ATM coverage graph for LND mesh" )
        ierr = iMOAB_CoverageGraph( &lndCouComm, cmpLndPID, cplLndPID, cplLndAtmPID, &cmplnd, &cpllnd,
                                    &cplatm );  // it happens over joint communicator
        CHECKIERR( ierr, "cannot recompute direct coverage graph for land" )
        POP_TIMER( lndCouComm, rankInLndComm )  // hijack this rank
    }

#ifdef ENABLE_OCNATM_COUPLING
    // load ocn and migrate to coupler
    repartitioner_scheme = 4; // will use RCB, ands store the cuts
    // ocean
    ierr = setup_component_coupler_meshes( cmpOcnPID, cmpocn, cplOcnPID, cplocn, &ocnComm, &ocnPEGroup, &couComm,
                                            &couPEGroup, &ocnCouComm, ocnFilename, readopts, nghlay, repartitioner_scheme );
    CHECKIERR( ierr, "Cannot load and migrate ocn mesh" )
    if( couComm != MPI_COMM_NULL )
    {
        char outputFileTgt5[] = "recvOcn5.h5m";
        PUSH_TIMER( couComm, "Write migrated OCN mesh on coupler PEs" )
        ierr = iMOAB_WriteMesh( cplOcnPID, outputFileTgt5, fileWriteOptions );
        CHECKIERR( ierr, "cannot write ocn mesh after receiving" )
        POP_TIMER( couComm, rankInCouComm )
    }

    // infer atm2 from ocean:
    // atm
    repartitioner_scheme = 5; // reuse the partition stored at previous step
    /*ierr =
        setup_component_coupler_meshes( cmpAtmPID, cmpatm, cplAtm2PID, cplatm2, &atmComm, &atmPEGroup, &couComm,
                                        &couPEGroup, &atmCouComm, atmFilename, readopts, nghlay, repartitioner_scheme );
    CHECKIERR( ierr, "Cannot load and migrate atm mesh for ocean context " )*/
    // we cannot use setup_component_coupler_meshes for ocn - atm2, because atm was already loaded
    if ( atmCouComm != MPI_COMM_NULL &&  5 == repartitioner_scheme ) // yes, we just set it to 5
    {
        // need to send the zoltan buffer from coupler root towards the component root (atm in this case)
        int is_fort = 0;
        ierr = iMOAB_RetrieveZBuffer( &atmPEGroup, &couPEGroup, &atmCouComm, &is_fort);
        CHECKIERR( ierr, "Cannot retrieve Zoltan buffer from coupler root towards atm root " )
    }
    if( atmComm != MPI_COMM_NULL )
    {
        // then send mesh to coupler pes, with this new RCB scheme , from ocean
        ierr = iMOAB_SendMesh( cmpAtmPID, &atmCouComm, &couPEGroup, &cplatm2,
                               &repartitioner_scheme );  // send to  coupler pes
        CHECKIERR( ierr, "cannot send elements for atm again" )
    }
    // now, receive mesh, on coupler communicator; second mesh 2, atm2
    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_ReceiveMesh( cplAtm2PID, &atmCouComm, &atmPEGroup,
                                  &cmpatm );  // receive from component
        CHECKIERR( ierr, "cannot receive elements on coupler app" )
    }
    if( atmComm != MPI_COMM_NULL )
    {
        // then send mesh to coupler pes, with this new RCB scheme , from ocean
        ierr = iMOAB_FreeSenderBuffers( cmpAtmPID, &cplatm2 );
        CHECKIERR( ierr, "cannot free buffers used to send atm second mesh" )
    }


    if( couComm != MPI_COMM_NULL )
    {
        // now compute intersection between OCNx and ATMx on coupler PEs
        ierr = iMOAB_RegisterApplication( "OCNATM", &couComm, &ocnatmid, cplOcnAtm2PID );
        CHECKIERR( ierr, "Cannot register ocen atm2 intx over coupler pes " )
    }

    if( couComm != MPI_COMM_NULL )
    {
        PUSH_TIMER( couComm, "Compute OCN-ATM mesh intersection" )
        ierr = iMOAB_ComputeMeshIntersectionOnSphere(
            cplOcnPID, cplAtm2PID, cplOcnAtm2PID );  // coverage mesh was computed here, for cplAtmPID, atm on coupler pes
        // basically, ocn was redistributed according to target (atm) partition, to "cover" the atm partitions
        // check if intx valid, write some h5m intx file
        CHECKIERR( ierr, "cannot compute intersection" )
        POP_TIMER( couComm, rankInCouComm )
    }

    if( ocnCouComm != MPI_COMM_NULL )
    {
        // now for the  intersection, ocn-atm; will be sending data from ocean to atm
        // the new graph will be for sending data from ocn comp to coverage mesh over atm;
        // it involves initial ocn app; cmpOcnPID; also migrated ocn mesh on coupler pes, cplOcnPID
        // results are in cplOcnAtm2PID, intx mesh; remapper also has some info about coverage mesh
        // after this, the sending of tags from ocn pes to coupler pes will use the new par comm graph, that has more
        // precise info about what to send for atm cover ; every time, we will
        //  use the element global id, which should uniquely identify the element
        PUSH_TIMER( ocnCouComm, "Compute ATM coverage graph for OCN mesh" )
        ierr = iMOAB_CoverageGraph( &ocnCouComm, cmpOcnPID, cplOcnPID, cplOcnAtm2PID, &cmpocn, &cplocn,
                                    &cplatm2 );  // it happens over joint communicator, ocean + coupler
        CHECKIERR( ierr, "cannot recompute direct coverage graph for atm" )
        POP_TIMER( ocnCouComm, rankInOcnComm )  // hijack this rank
    }
#endif


    int fMonotoneTypeID = 0, fVolumetric = 0, fValidate = 1, fNoConserve = 0, fNoBubble = 1, fInverseDistanceMap = 0;
    const char* weights_identifiers[2] = { "lndatm", "ocnatm" };
    int disc_orders[3]                 = { 4, 1, 1 };
    const char* disc_methods[3]        = { "cgll", "fv", "pcloud" };
    const char* dof_tag_names[3]       = { "GLOBAL_DOFS", "GLOBAL_ID", "GLOBAL_ID" };
    if( couComm != MPI_COMM_NULL )
    {
        // weight identifier lndatm
        // Compute the weights to project the solution from LND component to ATM component
        PUSH_TIMER( couComm, "Compute LND-ATM remapping weights" )
        ierr = iMOAB_ComputeScalarProjectionWeights(
            cplLndAtmPID, weights_identifiers[0], disc_methods[1], &disc_orders[1], disc_methods[1], &disc_orders[1],
            &fNoBubble, &fMonotoneTypeID, &fVolumetric, &fInverseDistanceMap, &fNoConserve, &fValidate, dof_tag_names[1],
            dof_tag_names[1] );
        CHECKIERR( ierr, "failed to compute remapping projection weights for ATM-LND scalar non-conservative field" );
        POP_TIMER( couComm, rankInCouComm )

        if (!skip_writing_map)
        {
            PUSH_TIMER( couComm, "write LND-ATM remapping weights" )

#ifdef  ENABLE_OCNATM_COUPLING
            const std::string lnd_atm_map_file_name = "lnd_atm_map_testE2.nc";
#else
            const std::string lnd_atm_map_file_name = "lnd_atm_map_testE1.nc";
#endif
            ierr = iMOAB_WriteMappingWeightsToFile( cplLndAtmPID, weights_identifiers[0],
                    lnd_atm_map_file_name.c_str());
            POP_TIMER( couComm, rankInCouComm )
        }
    }

#ifdef ENABLE_OCNATM_COUPLING
    if( couComm != MPI_COMM_NULL )
    {
        PUSH_TIMER( couComm, "Compute the projection weights with TempestRemap" )
        ierr = iMOAB_ComputeScalarProjectionWeights(
            cplOcnAtm2PID, weights_identifiers[1], disc_methods[1], &disc_orders[1],  // fv
            disc_methods[1], &disc_orders[1],                                        // fv
            &fNoBubble, &fMonotoneTypeID, &fVolumetric, &fInverseDistanceMap, &fNoConserve, &fValidate, dof_tag_names[1],
            dof_tag_names[1] );
        CHECKIERR( ierr, "cannot compute scalar projection weights" )
        POP_TIMER( couComm, rankInCouComm )
    }
#endif

    int tagIndex[2];
    int tagTypes[2]  = { DENSE_DOUBLE, DENSE_DOUBLE };
    int atmCompNDoFs = 1 /* FV disc_orders[0]*disc_orders[0] */, ocnCompNDoFs = 1 /*FV*/;

    const char* bottomTempField          = "T_proj";  // on ocean input files
    const char* bottomTempProjectedField = "T_proj2"; // on projected atm
    const char* bottomUVelField          = "u_proj";
    const char* bottomUVelProjectedField = "u_proj2";
    const char* bottomVVelField          = "v_proj";
    const char* bottomVVelProjectedField = "v_proj2";
    // coming from lnd to atm;  ( PS: will have to merge projection from ocn and land on atm )
    const char* bottomTempProjectedField3 = "T_proj3";
    const char* bottomUVelProjectedField3 = "u_proj3";
    const char* bottomVVelProjectedField3 = "v_proj3";
    // need to define tag storage for land; will use the names T_proj, u_proj, v_proj name, on the coupler side
    // use the same ndof and same size as ocnCompNDoFs (1) // FV actually
    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DefineTagStorage( cplLndPID, bottomTempField, &tagTypes[1], &ocnCompNDoFs, &tagIndex[1] );
        CHECKIERR( ierr, "failed to define the field tag T_proj" );
        ierr = iMOAB_DefineTagStorage( cplLndPID, bottomUVelField, &tagTypes[1], &ocnCompNDoFs, &tagIndex[1] );
        CHECKIERR( ierr, "failed to define the field tag u_proj" );
        ierr = iMOAB_DefineTagStorage( cplLndPID, bottomVVelField, &tagTypes[1], &ocnCompNDoFs, &tagIndex[1] );
        CHECKIERR( ierr, "failed to define the field tag v_proj" );
        // for projection on atm, we need new tags, to T_proj3, u_proj3, v_proj3, all on cplAtmPID
        ierr = iMOAB_DefineTagStorage( cplAtmPID, bottomTempProjectedField3, &tagTypes[1], &atmCompNDoFs, &tagIndex[1] );
        CHECKIERR( ierr, "failed to define the field tag T_proj3" );
        ierr = iMOAB_DefineTagStorage( cplAtmPID, bottomUVelProjectedField3, &tagTypes[1], &atmCompNDoFs, &tagIndex[1] );
        CHECKIERR( ierr, "failed to define the field tag u_proj3" );
        ierr = iMOAB_DefineTagStorage( cplAtmPID, bottomVVelProjectedField3, &tagTypes[1], &atmCompNDoFs, &tagIndex[1] );
        CHECKIERR( ierr, "failed to define the field tag v_proj3" );
    }

#ifdef ENABLE_OCNATM_COUPLING
    // ocn - atm coupling T_proj -> T_proj2
    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DefineTagStorage( cplOcnPID, bottomTempField, &tagTypes[0], &ocnCompNDoFs, &tagIndex[0] );
        CHECKIERR( ierr, "failed to define the field tag T_proj" );
        ierr = iMOAB_DefineTagStorage( cplAtm2PID, bottomTempProjectedField, &tagTypes[1], &atmCompNDoFs, &tagIndex[1] );
        CHECKIERR( ierr, "failed to define the field tag T_proj2" );

        ierr = iMOAB_DefineTagStorage( cplOcnPID, bottomUVelField, &tagTypes[0], &ocnCompNDoFs, &tagIndex[0] );
        CHECKIERR( ierr, "failed to define the field tag u_proj" );
        ierr = iMOAB_DefineTagStorage( cplAtm2PID, bottomUVelProjectedField, &tagTypes[1], &atmCompNDoFs, &tagIndex[1] );
        CHECKIERR( ierr, "failed to define the field tag u_proj2" );

        ierr = iMOAB_DefineTagStorage( cplOcnPID, bottomVVelField, &tagTypes[0], &ocnCompNDoFs, &tagIndex[0] );
        CHECKIERR( ierr, "failed to define the field tag v_proj" );
        ierr = iMOAB_DefineTagStorage( cplAtm2PID, bottomVVelProjectedField, &tagTypes[1], &atmCompNDoFs, &tagIndex[1] );
        CHECKIERR( ierr, "failed to define the field tag v_proj2" );
    }


    // need to make sure that the coverage mesh (created during intx method OCN - ATM) received the tag that need to be
    // projected to target (atm); the coverage mesh has only the ids; need to change the migrate method to
    // accommodate any ids;  now send a tag from original ocn (cmpOcnPID) towards migrated coverage mesh
    // (cplOcnPID), using the new coverage graph communicator

    // make the tag 0, to check we are actually sending needed data
    {
        if( cplOcnAppID >= 0 )
        {
            int nverts[3], nelem[3], nblocks[3], nsbc[3], ndbc[3];
            /*
             * Each process in the communicator will have access to a local mesh instance, which will contain the
             * original cells in the local partition and ghost entities. Number of vertices, primary cells, visible
             * blocks, number of sidesets and nodesets boundary conditions will be returned in numProcesses 3 arrays,
             * for local, ghost and total numbers.
             */
            ierr = iMOAB_GetMeshInfo( cplOcnPID, nverts, nelem, nblocks, nsbc, ndbc );
            CHECKIERR( ierr, "failed to get num primary elems" );
            int numAllElem = nelem[2];
            std::vector< double > vals;
            int storLeng = ocnCompNDoFs * numAllElem;
            vals.resize( storLeng );
            for( int k = 0; k < storLeng; k++ )
                vals[k] = 0.;
            int eetype = 1;
            ierr       = iMOAB_SetDoubleTagStorage( cplOcnPID, bottomTempField, &storLeng, &eetype, &vals[0] );
            CHECKIERR( ierr, "cannot make tag T_proj null" )
            ierr = iMOAB_SetDoubleTagStorage( cplOcnPID, bottomUVelField, &storLeng, &eetype, &vals[0] );
            CHECKIERR( ierr, "cannot make tag u_proj null" )
            ierr = iMOAB_SetDoubleTagStorage( cplOcnPID, bottomVVelField, &storLeng, &eetype, &vals[0] );
            CHECKIERR( ierr, "cannot make tag v_proj null" )
            // set the tag to 0
        }
    }
    if (analytic_field && (ocnComm != MPI_COMM_NULL) ) // we are on ocean pes
    {
        // cmpOcnPID, "T_proj;u_proj;v_proj;"
        ierr = iMOAB_DefineTagStorage( cmpOcnPID, bottomTempField, &tagTypes[0], &ocnCompNDoFs, &tagIndex[0] );
        CHECKIERR( ierr, "failed to define the field tag T_proj" );

        ierr = iMOAB_DefineTagStorage( cmpOcnPID, bottomUVelField, &tagTypes[0], &ocnCompNDoFs, &tagIndex[0] );
        CHECKIERR( ierr, "failed to define the field tag u_proj" );

        ierr = iMOAB_DefineTagStorage( cmpOcnPID, bottomVVelField, &tagTypes[0], &ocnCompNDoFs, &tagIndex[0] );
        CHECKIERR( ierr, "failed to define the field tag v_proj" );

        int nverts[3], nelem[3], nblocks[3], nsbc[3], ndbc[3];
        /*
         * Each process in the communicator will have access to a local mesh instance, which will contain the
         * original cells in the local partition and ghost entities. Number of vertices, primary cells, visible
         * blocks, number of sidesets and nodesets boundary conditions will be returned in numProcesses 3 arrays,
         * for local, ghost and total numbers.
         */
        ierr = iMOAB_GetMeshInfo( cmpOcnPID, nverts, nelem, nblocks, nsbc, ndbc );
        CHECKIERR( ierr, "failed to get num primary elems" );
        int numAllElem = nelem[2];
        std::vector< double > vals;
        int storLeng = ocnCompNDoFs * numAllElem;
        vals.resize( storLeng );
        for( int k = 0; k < storLeng; k++ )
            vals[k] = k;
        int eetype = 1;
        ierr       = iMOAB_SetDoubleTagStorage( cmpOcnPID, bottomTempField, &storLeng, &eetype, &vals[0] );
        CHECKIERR( ierr, "cannot make tag T_proj null" )
        ierr = iMOAB_SetDoubleTagStorage( cmpOcnPID, bottomUVelField, &storLeng, &eetype, &vals[0] );
        CHECKIERR( ierr, "cannot make tag u_proj null" )
        ierr = iMOAB_SetDoubleTagStorage( cmpOcnPID, bottomVVelField, &storLeng, &eetype, &vals[0] );
        CHECKIERR( ierr, "cannot make tag v_proj null" )
                    // set the tag to 0

    }
#endif



    MPI_Barrier( MPI_COMM_WORLD );
    // start a virtual loop for number of iterations
    for( int iters = 0; iters < n; iters++ )
    {
#ifdef VERBOSE
        if( ( atmComm != MPI_COMM_NULL ) && ( 1 == n ) )
        {
            char outputFileAtm[] = "AtmAtStartIt.h5m";
            ierr                 = iMOAB_WriteMesh( cmpAtmPID, outputFileAtm, fileWriteOptions );
            CHECKIERR( ierr, "cannot write AtmWithProj3.h5m" )
        }
#endif
        PUSH_TIMER( MPI_COMM_WORLD, "Send/receive data from land component to coupler in atm context" )
        if( lndComm != MPI_COMM_NULL )
        {
            // as always, use nonblocking sends
            // this is for projection to atm, from ocean:
            ierr = iMOAB_SendElementTag( cmpLndPID, "T_proj;u_proj;v_proj;", &lndCouComm, &cplatm );
            CHECKIERR( ierr, "cannot send tag values" )
        }
        if( couComm != MPI_COMM_NULL )
        {
            // receive on lnd on coupler pes, that was redistributed according to coverage
            ierr = iMOAB_ReceiveElementTag( cplLndPID, "T_proj;u_proj;v_proj;", &lndCouComm, &cplatm );
            CHECKIERR( ierr, "cannot receive tag values" )
        }

        // we can now free the sender buffers
        if( lndComm != MPI_COMM_NULL )
        {
            ierr = iMOAB_FreeSenderBuffers( cmpLndPID, &cplatm );  // context is for atm
            CHECKIERR( ierr, "cannot free buffers used to send lnd tag towards the coverage mesh for atm" )
        }
        POP_TIMER( MPI_COMM_WORLD, rankInGlobalComm )
#ifdef VERBOSE
        if( ( atmComm != MPI_COMM_NULL ) && ( 1 == n ) )
        {
            char outputFileAtm[] = "AtmAfterTagMigrToLand.h5m";
            ierr                 = iMOAB_WriteMesh( cmpAtmPID, outputFileAtm, fileWriteOptions );
            CHECKIERR( ierr, "cannot write AtmWithProj3.h5m" )
        }
#endif
        if( couComm != MPI_COMM_NULL )
        {
            const char* concat_fieldname  = "T_proj;u_proj;v_proj;";
            const char* concat_fieldnameT = "T_proj3;u_proj3;v_proj3;";

            /* We have the remapping weights now. Let us apply the weights onto the tag we defined
               on the source mesh and get the projection on the target mesh */
            PUSH_TIMER( couComm, "Apply Scalar projection weights" )
            ierr = iMOAB_ApplyScalarProjectionWeights( cplLndAtmPID, weights_identifiers[0], concat_fieldname,
                                                       concat_fieldnameT );
            CHECKIERR( ierr, "failed to compute projection weight application" );
            POP_TIMER( couComm, rankInCouComm )
            if( 1 == n )
            {
                char outputFileTgt[] = "fAtmOnCpl4.h5m";
                ierr = iMOAB_WriteMesh( cplAtmPID, outputFileTgt, fileWriteOptions );
                CHECKIERR( ierr, "failed to write fAtmOnCpl4.h5m " );
            }
        }
#ifdef VERBOSE
        if( ( atmComm != MPI_COMM_NULL ) && ( 1 == n ) )
        {
            char outputFileAtm[] = "AtmBeforeMigr.h5m";
            ierr                 = iMOAB_WriteMesh( cmpAtmPID, outputFileAtm, fileWriteOptions );
            CHECKIERR( ierr, "cannot write AtmWithProj3.h5m" )
        }
#endif
        // send the tag to atm pes, from atm mesh on coupler pes
        //   from couComm, using common joint comm atm_coupler
        // as always, use nonblocking sends
        // original graph: context is now the cpl
        PUSH_TIMER( MPI_COMM_WORLD, "Send/receive data from atm coupler to atm component with land proj" )
        if( couComm != MPI_COMM_NULL )
        {
            context_id = cmpatm;
            ierr = iMOAB_SendElementTag( cplAtmPID, "T_proj3;u_proj3;v_proj3;", &atmCouComm, &context_id );
            CHECKIERR( ierr, "cannot send tag values back to atm pes " )
        }

        // receive on component atm

        if( atmComm != MPI_COMM_NULL )
        {
            context_id = cplatm;
            ierr = iMOAB_ReceiveElementTag( cmpAtmPID, "T_proj3;u_proj3;v_proj3;", &atmCouComm, &context_id );
            CHECKIERR( ierr, "cannot receive tag values from atm mesh on coupler pes" )
        }

        if( couComm != MPI_COMM_NULL )
        {
            context_id = cmpatm;
            ierr = iMOAB_FreeSenderBuffers( cplAtmPID, &context_id );
            CHECKIERR( ierr, "cannot free buffers related to send tag" )
        }
        POP_TIMER( MPI_COMM_WORLD, rankInGlobalComm )
#ifdef VERBOSE
        if( ( atmComm != MPI_COMM_NULL ) && ( 1 == n ) )
        {
            char outputFileAtm[] = "AtmWithProj3.h5m";
            ierr                 = iMOAB_WriteMesh( cmpAtmPID, outputFileAtm, fileWriteOptions );
            CHECKIERR( ierr, "cannot write AtmWithProj3.h5m" )
        }
#endif

#ifdef ENABLE_OCNATM_COUPLING
        PUSH_TIMER( MPI_COMM_WORLD, "Send/receive data from ocn component to coupler in atm context" )
        if( ocnComm != MPI_COMM_NULL )
        {
          // as always, use nonblocking sends
          // this is for projection to atm, from ocean:
            ierr = iMOAB_SendElementTag( cmpOcnPID, "T_proj;u_proj;v_proj;", &ocnCouComm, &cplatm2 );
            CHECKIERR( ierr, "cannot send tag values" )
        }
        if( couComm != MPI_COMM_NULL )
        {
          // receive on ocn on coupler pes, that was redistributed according to coverage
            ierr = iMOAB_ReceiveElementTag( cplOcnPID, "T_proj;u_proj;v_proj;", &ocnCouComm, &cplatm2 );
            CHECKIERR( ierr, "cannot receive tag values" )
        }


        // we can now free the sender buffers
        if( ocnComm != MPI_COMM_NULL )
        {
            ierr = iMOAB_FreeSenderBuffers( cmpOcnPID, &cplatm2 );  // context is for atm
            CHECKIERR( ierr, "cannot free buffers used to send ocn tag towards the coverage mesh for atm" )
        }

        POP_TIMER( MPI_COMM_WORLD, rankInGlobalComm )
        if( couComm != MPI_COMM_NULL )
        {
            const char* concat_fieldname  = "T_proj;u_proj;v_proj;";
            const char* concat_fieldnameT = "T_proj2;u_proj2;v_proj2;";

         /* We have the remapping weights computed earlier, and te field. Let us apply the weights onto the tag
          * we defined  on the source mesh and get the projection on the target mesh */
            PUSH_TIMER( couComm, "Apply Scalar projection weights" )
            ierr = iMOAB_ApplyScalarProjectionWeights( cplOcnAtm2PID, weights_identifiers[1], concat_fieldname,
                                                    concat_fieldnameT );
            CHECKIERR( ierr, "failed to compute projection weight application" );
            POP_TIMER( couComm, rankInCouComm )
         // do not write if iters > 0)
           if( 1 == n )
           {
               char outputFileTgt[] = "fAtmOnCpl5.h5m";
               ierr = iMOAB_WriteMesh( cplAtm2PID, outputFileTgt, fileWriteOptions );
               CHECKIERR( ierr, "failed to write fAtmOnCpl3.h5m " );
           }
        }

        PUSH_TIMER( MPI_COMM_WORLD, "Send/receive data from atm2 coupler to atm component with ocean proj" )
        if( couComm != MPI_COMM_NULL )
        {
            context_id = cmpatm;
            ierr = iMOAB_SendElementTag( cplAtm2PID, "T_proj2;u_proj2;v_proj2;", &atmCouComm, &context_id );
            CHECKIERR( ierr, "cannot send tag values back to atm pes from ocean proj " )
        }

        // receive on component atm

        if( atmComm != MPI_COMM_NULL )
        {
            context_id = cplatm2;
            ierr = iMOAB_ReceiveElementTag( cmpAtmPID, "T_proj2;u_proj2;v_proj2;", &atmCouComm, &context_id );
            CHECKIERR( ierr, "cannot receive tag values from atm mesh on coupler pes, ocean context" )
        }

        if( couComm != MPI_COMM_NULL )
        {
            context_id = cmpatm;
            ierr = iMOAB_FreeSenderBuffers( cplAtm2PID, &context_id );
            CHECKIERR( ierr, "cannot free buffers related to send tag" )
        }
        POP_TIMER( MPI_COMM_WORLD, rankInGlobalComm )

#endif
    }
     // free up the MPI objects and finalize
    // order of registering:
 /*   ierr = iMOAB_RegisterApplication( "ATMX", &couComm, &cplatm, cplAtmPID );  // atm on coupler pes, land context
    ierr = iMOAB_RegisterApplication( "LNDX", &couComm, &cpllnd, cplLndPID );  // lnd on coupler pes
    ierr = iMOAB_RegisterApplication( "ATM2X", &couComm, &cplatm2, cplAtm2PID );  // atm on coupler pes, ocean context
    ierr = iMOAB_RegisterApplication( "OCNX", &couComm, &cplocn, cplOcnPID );  // ocn on coupler pes
    ierr = iMOAB_RegisterApplication( "ATM1", &atmComm, &cmpatm, cmpAtmPID );
    ierr = iMOAB_RegisterApplication( "OCN1", &ocnComm, &cmpocn, cmpOcnPID );
    ierr = iMOAB_RegisterApplication( "LND1", &lndComm, &cmplnd, cmpLndPID );
    ierr = iMOAB_RegisterApplication( "LNDATM", &couComm, &lndatmid, cplLndAtmPID );
    ierr = iMOAB_RegisterApplication( "OCNATM", &couComm, &ocnatmid, cplOcnAtm2PID ); */
    // de-registering in reverse order

#ifdef ENABLE_OCNATM_COUPLING
    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplOcnAtm2PID );
        CHECKIERR( ierr, "cannot deregister app ocn-atm2 " )
    }
#endif
    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplLndAtmPID );
        CHECKIERR( ierr, "cannot deregister app LNDATM" )
    }
    if( lndComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cmpLndPID );
        CHECKIERR( ierr, "cannot deregister app LND1" )
    }

#ifdef ENABLE_OCNATM_COUPLING
    if( ocnComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cmpOcnPID );
        CHECKIERR( ierr, "cannot deregister app OCN1" )
    }
#endif

    if( atmComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cmpAtmPID );
        CHECKIERR( ierr, "cannot deregister app ATM1" )
    }
#ifdef ENABLE_OCNATM_COUPLING
    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplOcnPID );
        CHECKIERR( ierr, "cannot deregister app OCNX" )
    }
    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplAtm2PID );
        CHECKIERR( ierr, "cannot deregister app ATM2X" )
    }

#endif

    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplLndPID );
        CHECKIERR( ierr, "cannot deregister app LNDX" )
    }

    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplAtmPID );
        CHECKIERR( ierr, "cannot deregister app ATMX" )
    }

    ierr = iMOAB_Finalize();
    CHECKIERR( ierr, "did not finalize iMOAB" )

    // free atm coupler group and comm
#ifdef ENABLE_OCNATM_COUPLING
    if( MPI_COMM_NULL != ocnCouComm ) MPI_Comm_free( &ocnCouComm );
    MPI_Group_free( &joinOcnCouGroup );
    if( MPI_COMM_NULL != ocnComm ) MPI_Comm_free( &ocnComm );
    MPI_Group_free( &ocnPEGroup );

#endif
    if( MPI_COMM_NULL != atmCouComm ) MPI_Comm_free( &atmCouComm );
    MPI_Group_free( &joinAtmCouGroup );
    if( MPI_COMM_NULL != atmComm ) MPI_Comm_free( &atmComm );


    if( MPI_COMM_NULL != lndComm ) MPI_Comm_free( &lndComm );
    // free land - coupler group and comm
    if( MPI_COMM_NULL != lndCouComm ) MPI_Comm_free( &lndCouComm );
    MPI_Group_free( &joinLndCouGroup );

    if( MPI_COMM_NULL != couComm ) MPI_Comm_free( &couComm );

    MPI_Group_free( &atmPEGroup );


    MPI_Group_free( &lndPEGroup );

    MPI_Group_free( &couPEGroup );
    MPI_Group_free( &jgroup );

    MPI_Finalize();

    return 0;
}
