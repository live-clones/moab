/*
 * This imoab_coupler test will simulate coupling between 3 components
 * 3 meshes will be loaded from 3 files (atm, ocean, lnd), and they will be migrated to
 * all processors (coupler pes); then, intx will be performed between migrated meshes
 * and weights will be generated, such that a field from one component will be transferred to
 * the other component
 * currently, the atm will send some data to be projected to ocean and land components
 *
 * first, intersect atm and ocn, and recompute comm graph 1 between atm and atm_cx, for ocn intx
 * second, intersect atm and lnd, and recompute comm graph 2 between atm and atm_cx for lnd intx

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

using namespace moab;

//#define GRAPH_INFO

#ifndef MOAB_HAVE_TEMPESTREMAP
#error The climate coupler test example requires MOAB configuration with TempestRemap
#endif

#define ENABLE_ATMOCN_COUPLING
#define ENABLE_ATMCPLOCN_COUPLING

#if( !defined( ENABLE_ATMOCN_COUPLING ) && !defined( ENABLE_ATMCPLOCN_COUPLING ) )
#error Enable either OCN (ENABLE_ATMOCN_COUPLING) and/or LND (ENABLE_ATMCPLOCN_COUPLING) for coupling
#endif

int main( int argc, char* argv[] )
{
    int ierr;
    int rankInGlobalComm, numProcesses;
    MPI_Group jgroup;
    std::string readopts( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS" );
    std::string readoptsLnd( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION" );

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

    std::string atmFilename = TestDir + "unittest/wholeATM_T.h5m";
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
#ifdef ENABLE_ATMOCN_COUPLING
    std::string ocnFilename = TestDir + "unittest/recMeshOcn.h5m";
    std::string baseline    = TestDir + "unittest/baseline1.txt";
    int rankInOcnComm       = -1;
    int cmpocn = 17, cplocn = 18,
        atmocnid = 618;  // component ids are unique over all pes, and established in advance;
#endif

#ifdef ENABLE_ATMCPLOCN_COUPLING
    int cplatm2 = 10,
        atm2ocnid = 610;  // component ids are unique over all pes, and established in advance;
#endif

    int rankInCouComm = -1;

    int nghlay = 0;  // number of ghost layers for loading the file
    std::vector< int > groupTasks;
    int startG1 = 0, startG2 = 0,
        endG1 = numProcesses - 1, endG2 = numProcesses - 1; // Support launch of imoab_coupler test on any combo of 2*x processes
    int startG4 = startG1, endG4 = endG1;  // these are for coupler layout
    int context_id = -1;                   // used now for freeing buffers

    // default: load atm on 2 proc, ocean on 2, land on 2; migrate to 2 procs, then compute intx
    // later, we need to compute weight matrix with tempestremap

    ProgOptions opts;
    opts.addOpt< std::string >( "atmosphere,t", "atm mesh filename (source)", &atmFilename );
#ifdef ENABLE_ATMOCN_COUPLING
    opts.addOpt< std::string >( "ocean,m", "ocean mesh filename (target)", &ocnFilename );
#endif
    opts.addOpt< int >( "startAtm,a", "start task for atmosphere layout", &startG1 );
    opts.addOpt< int >( "endAtm,b", "end task for atmosphere layout", &endG1 );
#ifdef ENABLE_ATMOCN_COUPLING
    opts.addOpt< int >( "startOcn,c", "start task for ocean layout", &startG2 );
    opts.addOpt< int >( "endOcn,d", "end task for ocean layout", &endG2 );
#endif

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
#ifdef ENABLE_ATMOCN_COUPLING
            "\n ocn file: " << ocnFilename << "\n     on tasks : " << startG2 << ":" << endG2 <<
#endif
            "\n  partitioning (0 trivial, 1 graph, 2 geometry) " << repartitioner_scheme << "\n  ";
    }

    // load files on 3 different communicators, groups
    // first groups has task 0, second group tasks 0 and 1
    // coupler will be on joint tasks, will be on a third group (0 and 1, again)
    // first groups has task 0, second group tasks 0 and 1
    // coupler will be on joint tasks, will be on a third group (0 and 1, again)
    MPI_Group atmPEGroup;
    MPI_Comm atmComm;
    ierr = create_group_and_comm( startG1, endG1, jgroup, &atmPEGroup, &atmComm );
    CHECKIERR( ierr, "Cannot create atm MPI group and communicator " )

#ifdef ENABLE_ATMOCN_COUPLING
    MPI_Group ocnPEGroup;
    MPI_Comm ocnComm;
    ierr = create_group_and_comm( startG2, endG2, jgroup, &ocnPEGroup, &ocnComm );
    CHECKIERR( ierr, "Cannot create ocn MPI group and communicator " )
#endif

    // we will always have a coupler
    MPI_Group couPEGroup;
    MPI_Comm couComm;
    ierr = create_group_and_comm( startG4, endG4, jgroup, &couPEGroup, &couComm );
    CHECKIERR( ierr, "Cannot create cpl MPI group and communicator " )

    // atm_coupler
    MPI_Group joinAtmCouGroup;
    MPI_Comm atmCouComm;
    ierr = create_joint_comm_group( atmPEGroup, couPEGroup, &joinAtmCouGroup, &atmCouComm );
    CHECKIERR( ierr, "Cannot create joint atm cou communicator" )

#ifdef ENABLE_ATMOCN_COUPLING
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
    int cplAtmAppID       = -1;            // -1 means it is not initialized
    iMOAB_AppID cplAtmPID = &cplAtmAppID;  // atm on coupler PEs
#ifdef ENABLE_ATMOCN_COUPLING
    int cmpOcnAppID       = -1;
    iMOAB_AppID cmpOcnPID = &cmpOcnAppID;        // ocn
    int cplOcnAppID = -1, cplAtmOcnAppID = -1;   // -1 means it is not initialized
    iMOAB_AppID cplOcnPID    = &cplOcnAppID;     // ocn on coupler PEs
    iMOAB_AppID cplAtmOcnPID = &cplAtmOcnAppID;  // intx atm -ocn on coupler PEs
#endif

#ifdef ENABLE_ATMCPLOCN_COUPLING
    int cplAtm2AppID          = -1;                // -1 means it is not initialized
    iMOAB_AppID cplAtm2PID    = &cplAtm2AppID;     // atm on second coupler PEs
    int cplAtm2OcnAppID       = -1;                // -1 means it is not initialized
    iMOAB_AppID cplAtm2OcnPID = &cplAtm2OcnAppID;  // intx atm - lnd on coupler PEs
#endif

    if( couComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( couComm, &rankInCouComm );
        // Register all the applications on the coupler PEs
        ierr = iMOAB_RegisterApplication( "ATMX", &couComm, &cplatm,
                                          cplAtmPID );  // atm on coupler pes
        CHECKIERR( ierr, "Cannot register ATM over coupler PEs" )
#ifdef ENABLE_ATMOCN_COUPLING
        ierr = iMOAB_RegisterApplication( "OCNX", &couComm, &cplocn,
                                          cplOcnPID );  // ocn on coupler pes
        CHECKIERR( ierr, "Cannot register OCN over coupler PEs" )
#endif
#ifdef ENABLE_ATMCPLOCN_COUPLING
        ierr = iMOAB_RegisterApplication( "ATMX2", &couComm, &cplatm2,
                                          cplAtm2PID );  // second atm on coupler pes
        CHECKIERR( ierr, "Cannot register second ATM over coupler PEs" )
#endif
    }

    if( atmComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( atmComm, &rankInAtmComm );
        ierr = iMOAB_RegisterApplication( "ATM1", &atmComm, &cmpatm, cmpAtmPID );
        CHECKIERR( ierr, "Cannot register ATM App" )
    }

#ifdef ENABLE_ATMOCN_COUPLING
    if( ocnComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( ocnComm, &rankInOcnComm );
        ierr = iMOAB_RegisterApplication( "OCN1", &ocnComm, &cmpocn, cmpOcnPID );
        CHECKIERR( ierr, "Cannot register OCN App" )
    }
#endif

    // atm
    ierr =
        setup_component_coupler_meshes( cmpAtmPID, cmpatm, cplAtmPID, cplatm, &atmComm, &atmPEGroup, &couComm,
                                        &couPEGroup, &atmCouComm, atmFilename, readopts, nghlay, repartitioner_scheme );
    CHECKIERR( ierr, "Cannot load and migrate atm mesh" )

#ifdef ENABLE_ATMOCN_COUPLING
    // ocean
    ierr =
        setup_component_coupler_meshes( cmpOcnPID, cmpocn, cplOcnPID, cplocn, &ocnComm, &ocnPEGroup, &couComm,
                                        &couPEGroup, &ocnCouComm, ocnFilename, readopts, nghlay, repartitioner_scheme );
    CHECKIERR( ierr, "Cannot load and migrate ocn mesh" )

#endif  // #ifdef ENABLE_ATMOCN_COUPLING

#ifdef ENABLE_ATMCPLOCN_COUPLING

    if( atmComm != MPI_COMM_NULL )
    {
        // then send mesh to second coupler pes
        ierr = iMOAB_SendMesh( cmpAtmPID, &atmCouComm, &couPEGroup, &cplatm2,
                               &repartitioner_scheme );  // send to  coupler pes
        CHECKIERR( ierr, "cannot send elements to coupler-2" )
    }
    // now, receive mesh, on coupler communicator; first mesh 1, atm
    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_ReceiveMesh( cplAtm2PID, &atmCouComm, &atmPEGroup,
                                  &cmpatm );  // receive from component
        CHECKIERR( ierr, "cannot receive elements on coupler app" )
    }

    // we can now free the sender buffers
    if( atmComm != MPI_COMM_NULL )
    {
        int context_id = cplatm;
        ierr           = iMOAB_FreeSenderBuffers( cmpAtmPID, &context_id );
        CHECKIERR( ierr, "cannot free buffers used to send atm mesh" )
    }

    if( couComm != MPI_COMM_NULL && 1 == n )
    {  // write only for n==1 case
        char outputFileLnd[] = "recvAtm2.h5m";
        ierr                 = iMOAB_WriteMesh( cplAtm2PID, outputFileLnd, fileWriteOptions );
        CHECKIERR( ierr, "cannot write second atm mesh after receiving" )
    }

#endif  // #ifdef ENABLE_ATMCPLOCN_COUPLING

    MPI_Barrier( MPI_COMM_WORLD );

#ifdef ENABLE_ATMOCN_COUPLING
    if( couComm != MPI_COMM_NULL )
    {
        // now compute intersection between OCNx and ATMx on coupler PEs
        ierr = iMOAB_RegisterApplication( "ATMOCN", &couComm, &atmocnid, cplAtmOcnPID );
        CHECKIERR( ierr, "Cannot register atm/ocn intersection over coupler pes " )
    }
#endif
#ifdef ENABLE_ATMCPLOCN_COUPLING
    if( couComm != MPI_COMM_NULL )
    {
        // now compute intersection between LNDx and ATMx on coupler PEs
        ierr = iMOAB_RegisterApplication( "ATM2OCN", &couComm, &atm2ocnid, cplAtm2OcnPID );
        CHECKIERR( ierr, "Cannot register atm2/ocn intersection over coupler pes " )
    }
#endif

    int disc_orders[3]                       = { 4, 1, 1 };
    const std::string weights_identifiers[2] = { "scalar", "scalar-pc" };
    const std::string disc_methods[3]        = { "cgll", "fv", "pcloud" };
    const std::string dof_tag_names[3]       = { "GLOBAL_DOFS", "GLOBAL_ID", "GLOBAL_ID" };
#ifdef ENABLE_ATMOCN_COUPLING
    if( couComm != MPI_COMM_NULL )
    {
        PUSH_TIMER( "Compute ATM-OCN mesh intersection" )
        ierr = iMOAB_ComputeMeshIntersectionOnSphere(
            cplAtmPID, cplOcnPID,
            cplAtmOcnPID );  // coverage mesh was computed here, for cplAtmPID, atm on coupler pes
        // basically, atm was redistributed according to target (ocean) partition, to "cover" the
        // ocean partitions check if intx valid, write some h5m intx file
        CHECKIERR( ierr, "cannot compute intersection for atm/ocn" )
        POP_TIMER( couComm, rankInCouComm )
#ifdef VERBOSE
        char prefix[] = "intx_atmocn";
        ierr          = iMOAB_WriteLocalMesh( cplAtmOcnPID, prefix, strlen( prefix ) );
        CHECKIERR( ierr, "failed to write local intx mesh" );
#endif
    }

    if( atmCouComm != MPI_COMM_NULL )
    {
        // the new graph will be for sending data from atm comp to coverage mesh;
        // it involves initial atm app; cmpAtmPID; also migrate atm mesh on coupler pes, cplAtmPID
        // results are in cplAtmOcnPID, intx mesh; remapper also has some info about coverage mesh
        // after this, the sending of tags from atm pes to coupler pes will use the new par comm
        // graph, that has more precise info about what to send for ocean cover ; every time, we
        // will
        //  use the element global id, which should uniquely identify the element
        PUSH_TIMER( "Compute OCN coverage graph for ATM mesh" )
        ierr = iMOAB_CoverageGraph( &atmCouComm, cmpAtmPID, cplAtmPID, cplAtmOcnPID, &cmpatm, &cplatm,
                                    &cplocn );  // it happens over joint communicator
        CHECKIERR( ierr, "cannot recompute direct coverage graph for ocean" )
        POP_TIMER( atmCouComm, rankInAtmComm )  // hijack this rank
    }
#endif

#ifdef ENABLE_ATMCPLOCN_COUPLING
    if( couComm != MPI_COMM_NULL )
    {
        PUSH_TIMER( "Compute ATM-OCN mesh intersection" )
        ierr = iMOAB_ComputeMeshIntersectionOnSphere(
            cplAtm2PID, cplOcnPID,
            cplAtm2OcnPID );  // coverage mesh was computed here, for cplAtmPID, atm on coupler pes
        // basically, atm was redistributed according to target (ocean) partition, to "cover" the
        // ocean partitions check if intx valid, write some h5m intx file
        CHECKIERR( ierr, "cannot compute intersection for atm2/ocn" )
        POP_TIMER( couComm, rankInCouComm )
    // }
    // if( atmCouComm != MPI_COMM_NULL )
    // {
        // the new graph will be for sending data from atm comp to coverage mesh for land mesh;
        // it involves initial atm app; cmpAtmPID; also migrate atm mesh on coupler pes, cplAtmPID
        // results are in cplAtmLndPID, intx mesh; remapper also has some info about coverage mesh
        // after this, the sending of tags from atm pes to coupler pes will use the new par comm
        // graph, that has more precise info about what to send (specifically for land cover); every
        // time,
        /// we will use the element global id, which should uniquely identify the element
        PUSH_TIMER( "Compute OCN coverage graph for ATM2 mesh" )
        // Context: cplatm2 already holds the comm-graph for communicating between atm component and coupler2
        // We just need to create a comm graph to internally transfer data from coupler atm to coupler ocean
        // ierr = iMOAB_CoverageGraph( &couComm, cplAtm2PID, cplAtm2OcnPID, cplAtm2OcnPID, &cplatm2, &atm2ocnid,
                                    // &cplocn );  // it happens over joint communicator
        int type1 = 1;
        int type2 = 1;
        ierr      = iMOAB_ComputeCommGraph( cplAtm2PID, cplAtm2OcnPID, &couComm, &couPEGroup, &couPEGroup, &type1, &type2,
                                            &cplatm2, &atm2ocnid );
        CHECKIERR( ierr, "cannot recompute direct coverage graph for ocean from atm2" )
        POP_TIMER( couComm, rankInCouComm )  // hijack this rank
    }
#endif

    MPI_Barrier( MPI_COMM_WORLD );

    int fMonotoneTypeID = 0, fVolumetric = 0, fValidate = 0, fNoConserve = 0, fNoBubble = 1, fInverseDistanceMap = 0;

#ifdef ENABLE_ATMOCN_COUPLING
#ifdef VERBOSE
    if( couComm != MPI_COMM_NULL && 1 == n )
    {                                    // write only for n==1 case
        char serialWriteOptions[] = "";  // for writing in serial
        std::stringstream outf;
        outf << "intxAtmOcn_" << rankInCouComm << ".h5m";
        std::string intxfile = outf.str();  // write in serial the intx file, for debugging
        ierr                 = iMOAB_WriteMesh( cplAtmOcnPID, intxfile.c_str(), serialWriteOptions );
        CHECKIERR( ierr, "cannot write intx file result" )
    }
#endif

    if( couComm != MPI_COMM_NULL )
    {
        PUSH_TIMER( "Compute the projection weights with TempestRemap" )
        ierr =
            iMOAB_ComputeScalarProjectionWeights( cplAtmOcnPID, weights_identifiers[0].c_str(), disc_methods[0].c_str(),
                                                  &disc_orders[0], disc_methods[1].c_str(), &disc_orders[1], &fNoBubble,
                                                  &fMonotoneTypeID, &fVolumetric, &fInverseDistanceMap, &fNoConserve,
                                                  &fValidate, dof_tag_names[0].c_str(), dof_tag_names[1].c_str() );
        CHECKIERR( ierr, "cannot compute scalar projection weights" )
        POP_TIMER( couComm, rankInCouComm )

        // Let us now write the map file to disk and then read it back to test the I/O API in iMOAB
#ifdef MOAB_HAVE_NETCDF
        {
            const std::string atmocn_map_file_name = "atm_ocn_map.nc";
            ierr = iMOAB_WriteMappingWeightsToFile( cplAtmOcnPID, weights_identifiers[0].c_str(),
                                                    atmocn_map_file_name.c_str() );
            CHECKIERR( ierr, "failed to write map file to disk" );

            const std::string intx_from_file_identifier = "map-from-file";
            int dummyCpl = -1;
            int dummy_rowcol = -1;
            int dummyType = 0;
            ierr = iMOAB_LoadMappingWeightsFromFile( cplAtmOcnPID, &dummyCpl, &dummy_rowcol, &dummyType,
                 intx_from_file_identifier.c_str(), atmocn_map_file_name.c_str() );
            CHECKIERR( ierr, "failed to load map file from disk" );
        }
#endif
    }

#endif

    MPI_Barrier( MPI_COMM_WORLD );

#ifdef ENABLE_ATMCPLOCN_COUPLING
    if( couComm != MPI_COMM_NULL )
    {
        PUSH_TIMER( "Compute the projection weights with TempestRemap for atm2/ocn" )
        ierr =
            iMOAB_ComputeScalarProjectionWeights( cplAtm2OcnPID, weights_identifiers[0].c_str(), disc_methods[0].c_str(),
                                                  &disc_orders[0], disc_methods[1].c_str(), &disc_orders[1], &fNoBubble,
                                                  &fMonotoneTypeID, &fVolumetric, &fInverseDistanceMap, &fNoConserve,
                                                  &fValidate, dof_tag_names[0].c_str(), dof_tag_names[1].c_str() );
        CHECKIERR( ierr, "cannot compute scalar projection weights for atm2/ocn" )
        POP_TIMER( couComm, rankInCouComm )

        // Let us now write the map file to disk and then read it back to test the I/O API in iMOAB
#ifdef MOAB_HAVE_NETCDF
        {
            const std::string atmocn_map_file_name = "atm2_ocn_map.nc";
            ierr = iMOAB_WriteMappingWeightsToFile( cplAtm2OcnPID, weights_identifiers[0].c_str(),
                                                    atmocn_map_file_name.c_str() );
            CHECKIERR( ierr, "failed to write map file to disk" );

            const std::string intx_from_file_identifier = "map2-from-file";
            int dummyCpl                                = -1;
            int dummy_rowcol                            = -1;
            int dummyType                               = 0;
            ierr = iMOAB_LoadMappingWeightsFromFile( cplAtm2OcnPID, &dummyCpl, &dummy_rowcol, &dummyType,
                                                     intx_from_file_identifier.c_str(), atmocn_map_file_name.c_str() );
            CHECKIERR( ierr, "failed to load map file from disk" );
        }
#endif
    }
#endif

    int tagIndex[2];
    int tagTypes[2]  = { DENSE_DOUBLE, DENSE_DOUBLE };
    int atmCompNDoFs = disc_orders[0] * disc_orders[0], ocnCompNDoFs = 1 /*FV*/;

    const char* bottomFields          = "a2oTbot:a2oUbot:a2oVbot";
    const char* bottomProjectedFields = "a2oTbot_proj:a2oUbot_proj:a2oVbot_proj";
    const char* bottomProjectedFields2 = "a2oT2bot_src:a2oU2bot_src:a2oV2bot_src";
    const char* bottomProjectedFields3 = "a2oT2bot_proj:a2oU2bot_proj:a2oV2bot_proj";

    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DefineTagStorage( cplAtmPID, bottomFields, &tagTypes[0], &atmCompNDoFs, &tagIndex[0] );
        CHECKIERR( ierr, "failed to define the field tag a2oTbot" );
        ierr = iMOAB_DefineTagStorage( cplAtm2PID, bottomFields, &tagTypes[0], &atmCompNDoFs, &tagIndex[0] );
        CHECKIERR( ierr, "failed to define the field tag a2oTbot" );

#ifdef ENABLE_ATMOCN_COUPLING
        ierr = iMOAB_DefineTagStorage( cplOcnPID, bottomProjectedFields, &tagTypes[1], &ocnCompNDoFs, &tagIndex[1] );
        CHECKIERR( ierr, "failed to define the field tag a2oTbot_proj" );
#endif
#ifdef ENABLE_ATMCPLOCN_COUPLING
        ierr = iMOAB_DefineTagStorage( cplAtm2PID, bottomProjectedFields2, &tagTypes[0], &atmCompNDoFs, &tagIndex[0] );
        CHECKIERR( ierr, "failed to define the field tag a2oT2bot_proj" );
        ierr = iMOAB_DefineTagStorage( cplOcnPID, bottomProjectedFields3, &tagTypes[1], &ocnCompNDoFs, &tagIndex[1] );
        CHECKIERR( ierr, "failed to define the field tag a2oT2bot_proj" );
#endif
    }

    // need to make sure that the coverage mesh (created during intx method) received the tag that
    // need to be projected to target so far, the coverage mesh has only the ids and global dofs;
    // need to change the migrate method to accommodate any GLL tag
    // now send a tag from original atmosphere (cmpAtmPID) towards migrated coverage mesh
    // (cplAtmPID), using the new coverage graph communicator

    // make the tag 0, to check we are actually sending needed data
    {
        if( cplAtmAppID >= 0 )
        {
            int nverts[3], nelem[3], nblocks[3], nsbc[3], ndbc[3];
            /*
             * Each process in the communicator will have access to a local mesh instance, which
             * will contain the original cells in the local partition and ghost entities. Number of
             * vertices, primary cells, visible blocks, number of sidesets and nodesets boundary
             * conditions will be returned in numProcesses 3 arrays, for local, ghost and total
             * numbers.
             */
            ierr = iMOAB_GetMeshInfo( cplAtmPID, nverts, nelem, nblocks, nsbc, ndbc );
            CHECKIERR( ierr, "failed to get num primary elems" );
            int numAllElem = nelem[2];
            std::vector< double > vals;
            int storLeng = atmCompNDoFs * numAllElem *3; // 3 tags
            int eetype   = 1;

            vals.resize( storLeng );
            for( int k = 0; k < storLeng; k++ )
                vals[k] = 0.;

            ierr = iMOAB_SetDoubleTagStorage( cplAtmPID, bottomFields, &storLeng, &eetype, &vals[0] );
            CHECKIERR( ierr, "cannot make tag nul" )
            // set the tag to 0
        }
        if( cplAtm2AppID >= 0 )
        {
            int nverts[3], nelem[3], nblocks[3], nsbc[3], ndbc[3];
            /*
             * Each process in the communicator will have access to a local mesh instance, which
             * will contain the original cells in the local partition and ghost entities. Number of
             * vertices, primary cells, visible blocks, number of sidesets and nodesets boundary
             * conditions will be returned in numProcesses 3 arrays, for local, ghost and total
             * numbers.
             */
            ierr = iMOAB_GetMeshInfo( cplAtm2PID, nverts, nelem, nblocks, nsbc, ndbc );
            CHECKIERR( ierr, "failed to get num primary elems" );
            int numAllElem = nelem[2];
            std::vector< double > vals;
            int storLeng = atmCompNDoFs * numAllElem * 3;  // 3 tags
            int eetype   = 1;

            vals.resize( storLeng );
            for( int k = 0; k < storLeng; k++ )
                vals[k] = 0.;

            ierr = iMOAB_SetDoubleTagStorage( cplAtm2PID, bottomProjectedFields2, &storLeng, &eetype, &vals[0] );
            CHECKIERR( ierr, "cannot make tag nul" )
            // set the tag to 0
        }
    }

    // start a virtual loop for number of iterations
    for( int iters = 0; iters < n; iters++ )
    {
#ifdef ENABLE_ATMOCN_COUPLING
        PUSH_TIMER( "Send/receive data from atm component to coupler in ocn context" )
        if( atmComm != MPI_COMM_NULL )
        {
            // as always, use nonblocking sends
            // this is for projection to ocean:
            ierr = iMOAB_SendElementTag( cmpAtmPID, bottomFields, &atmCouComm, &cplocn );
            CHECKIERR( ierr, "cannot send tag values" )
        }
        if( couComm != MPI_COMM_NULL )
        {
            // receive on atm on coupler pes, that was redistributed according to coverage
            ierr = iMOAB_ReceiveElementTag( cplAtmPID, bottomFields, &atmCouComm, &cplocn );
            CHECKIERR( ierr, "cannot receive tag values" )
        }
        POP_TIMER( MPI_COMM_WORLD, rankInGlobalComm )

        // we can now free the sender buffers
        if( atmComm != MPI_COMM_NULL )
        {
            ierr = iMOAB_FreeSenderBuffers( cmpAtmPID, &cplocn );  // context is for ocean
            CHECKIERR( ierr, "cannot free buffers used to resend atm tag towards the coverage mesh" )
        }
#ifdef VERBOSE
        if( couComm != MPI_COMM_NULL && 1 == n )
        {
            // write only for n==1 case
            char outputFileRecvd[] = "recvAtmCoupOcn.h5m";
            ierr                   = iMOAB_WriteMesh( cplAtmPID, outputFileRecvd, fileWriteOptions );
            CHECKIERR( ierr, "could not write recvAtmCoupOcn.h5m to disk" )
        }
#endif

        if( couComm != MPI_COMM_NULL )
        {
            /* We have the remapping weights now. Let us apply the weights onto the tag we defined
               on the source mesh and get the projection on the target mesh */
            PUSH_TIMER( "Apply Scalar projection weights" )
            ierr = iMOAB_ApplyScalarProjectionWeights( cplAtmOcnPID, weights_identifiers[0].c_str(), bottomFields,
                                                       bottomProjectedFields );
            CHECKIERR( ierr, "failed to compute projection weight application" );
            POP_TIMER( couComm, rankInCouComm )
            if( 1 == n )  // write only for n==1 case
            {
                char outputFileTgt[] = "fOcnOnCpl.h5m";
                ierr                 = iMOAB_WriteMesh( cplOcnPID, outputFileTgt, fileWriteOptions );
                CHECKIERR( ierr, "could not write fOcnOnCpl.h5m to disk" )
            }
        }

        // send the projected tag back to ocean pes, with send/receive tag
        if( ocnComm != MPI_COMM_NULL )
        {
            int tagIndexIn2;
            ierr =
                iMOAB_DefineTagStorage( cmpOcnPID, bottomProjectedFields, &tagTypes[1], &ocnCompNDoFs, &tagIndexIn2 );
            CHECKIERR( ierr, "failed to define the field tag for receiving back the tags "
                             "a2oTbot_proj, a2oUbot_proj, a2oVbot_proj on ocn pes" );
        }
        // send the tag to ocean pes, from ocean mesh on coupler pes
        //   from couComm, using common joint comm ocn_coupler
        // as always, use nonblocking sends
        // original graph (context is -1_
        if( couComm != MPI_COMM_NULL )
        {
            // need to use ocean comp id for context
            context_id = cmpocn;  // id for ocean on comp
            ierr =
                iMOAB_SendElementTag( cplOcnPID, bottomProjectedFields, &ocnCouComm, &context_id );
            CHECKIERR( ierr, "cannot send tag values back to ocean pes" )
        }

        // receive on component 2, ocean
        if( ocnComm != MPI_COMM_NULL )
        {
            context_id = cplocn;  // id for ocean on coupler
            ierr       = iMOAB_ReceiveElementTag( cmpOcnPID, bottomProjectedFields, &ocnCouComm,
                                                  &context_id );
            CHECKIERR( ierr, "cannot receive tag values from ocean mesh on coupler pes" )
        }

        MPI_Barrier( MPI_COMM_WORLD );

        if( couComm != MPI_COMM_NULL )
        {
            context_id = cmpocn;
            ierr       = iMOAB_FreeSenderBuffers( cplOcnPID, &context_id );
            CHECKIERR( ierr, "cannot free send/receive buffers for OCN context" )
        }
        if( ocnComm != MPI_COMM_NULL && 1 == n )  // write only for n==1 case
        {
            char outputFileOcn[] = "OcnWithProj.h5m";
            ierr                 = iMOAB_WriteMesh( cmpOcnPID, outputFileOcn, fileWriteOptions );
            CHECKIERR( ierr, "could not write OcnWithProj.h5m to disk" )
            // test results only for n == 1, for bottomTempProjectedField
            if( !no_regression_test )
            {
                // the same as remap test
                // get temp field on ocean, from conservative, the global ids, and dump to the baseline file
                // first get GlobalIds from ocn, and fields:
                int nverts[3], nelem[3];
                ierr = iMOAB_GetMeshInfo( cmpOcnPID, nverts, nelem, 0, 0, 0 );
                CHECKIERR( ierr, "failed to get ocn mesh info" );
                std::vector< int > gidElems;
                gidElems.resize( nelem[2] );
                std::vector< double > tempElems;
                tempElems.resize( nelem[2] );
                // get global id storage
                const std::string GidStr = "GLOBAL_ID";  // hard coded too
                int tag_type = DENSE_INTEGER, ncomp = 1, tagInd = 0;
                ierr = iMOAB_DefineTagStorage( cmpOcnPID, GidStr.c_str(), &tag_type, &ncomp, &tagInd );
                CHECKIERR( ierr, "failed to define global id tag" );

                int ent_type = 1;
                ierr         = iMOAB_GetIntTagStorage( cmpOcnPID, GidStr.c_str(), &nelem[2], &ent_type, &gidElems[0] );
                CHECKIERR( ierr, "failed to get global ids" );
                ierr = iMOAB_GetDoubleTagStorage( cmpOcnPID, "a2oTbot_proj", &nelem[2], &ent_type,
                                                  &tempElems[0] );
                CHECKIERR( ierr, "failed to get temperature field" );
                int err_code = 1;
                check_baseline_file( baseline, gidElems, tempElems, 1.e-9, err_code );
                if( 0 == err_code )
                    std::cout << " passed baseline test atm2ocn on ocean task " << rankInOcnComm << "\n";
            }
        }
#endif

#ifdef ENABLE_ATMCPLOCN_COUPLING
        PUSH_TIMER( "Send/receive data from atm2 component to coupler of all data" )
        if( atmComm != MPI_COMM_NULL )
        {
            // as always, use nonblocking sends
            // this is for projection to ocean:
            ierr = iMOAB_SendElementTag( cmpAtmPID, bottomFields, &atmCouComm, &cplatm2 );
            CHECKIERR( ierr, "cannot send tag values" )
        }
        if( couComm != MPI_COMM_NULL )
        {
            // receive on atm on coupler pes, that was redistributed according to coverage
            ierr = iMOAB_ReceiveElementTag( cplAtm2PID, bottomFields, &atmCouComm, &cmpatm );
            CHECKIERR( ierr, "cannot receive tag values" )
        }
        POP_TIMER( MPI_COMM_WORLD, rankInGlobalComm )

        // we can now free the sender buffers
        if( atmComm != MPI_COMM_NULL )
        {
            ierr = iMOAB_FreeSenderBuffers( cmpAtmPID, &cplatm );  // context is for ocean
            CHECKIERR( ierr, "cannot free buffers used to resend atm tag towards the coverage mesh" )
        }
// #ifdef VERBOSE
        if( couComm != MPI_COMM_NULL && 1 == n )
        {  // write only for n==1 case
            char outputFileRecvd[] = "recvAtm2CoupFull.h5m";
            ierr                   = iMOAB_WriteMesh( cplAtm2PID, outputFileRecvd, fileWriteOptions );
            CHECKIERR( ierr, "could not write recvAtmCoupLnd.h5m to disk" )
        }
// #endif

        PUSH_TIMER( "Send/receive data from atm2 coupler to ocean coupler based on coverage data" )
        if( couComm != MPI_COMM_NULL )
        {
            // as always, use nonblocking sends
            // this is for projection to ocean:
            ierr = iMOAB_SendElementTag( cplAtm2PID, bottomFields, &couComm, &atm2ocnid );
            CHECKIERR( ierr, "cannot send tag values" )

            // receive on atm on coupler pes, that was redistributed according to coverage
            ierr = iMOAB_ReceiveElementTag( cplAtm2PID, bottomProjectedFields2, &couComm, &atm2ocnid );
            CHECKIERR( ierr, "cannot receive tag values" )
        }
        POP_TIMER( MPI_COMM_WORLD, rankInGlobalComm )

        // we can now free the sender buffers
        if( couComm != MPI_COMM_NULL )
        {
            ierr = iMOAB_FreeSenderBuffers( cmpAtmPID, &cplocn );  // context is for ocean
            CHECKIERR( ierr, "cannot free buffers used to resend atm tag towards the coverage mesh" )
        }
// #ifdef VERBOSE
        if( couComm != MPI_COMM_NULL && 1 == n )
        {  // write only for n==1 case
            char outputFileRecvd[] = "recvAtm2CoupOcnCtx.h5m";
            ierr                   = iMOAB_WriteMesh( cplAtm2PID, outputFileRecvd, fileWriteOptions );
            CHECKIERR( ierr, "could not write recvAtmCoupLnd.h5m to disk" )
        }
// #endif

        if( couComm != MPI_COMM_NULL )
        {
            /* We have the remapping weights now. Let us apply the weights onto the tag we defined
               on the source mesh and get the projection on the target mesh */
            PUSH_TIMER( "Apply Scalar projection weights" )
            ierr = iMOAB_ApplyScalarProjectionWeights( cplAtm2OcnPID, weights_identifiers[0].c_str(),
                                                       bottomProjectedFields2, bottomProjectedFields3 );
            CHECKIERR( ierr, "failed to compute projection weight application" );
            POP_TIMER( couComm, rankInCouComm )
            if( 1 == n )  // write only for n==1 case
            {
                char outputFileTgt[] = "fOcnOnCpl2.h5m";
                ierr                 = iMOAB_WriteMesh( cplOcnPID, outputFileTgt, fileWriteOptions );
                CHECKIERR( ierr, "could not write the second fOcnOnCpl.h5m to disk" )
            }
        }

        std::cout << "applied scalar projection\n";
        // send the projected tag back to ocean pes, with send/receive tag
        if( ocnComm != MPI_COMM_NULL )
        {
            int tagIndexIn2;
            ierr = iMOAB_DefineTagStorage( cmpOcnPID, bottomProjectedFields3, &tagTypes[1],
                                           &ocnCompNDoFs, &tagIndexIn2 );
            CHECKIERR( ierr, "failed to define the field tag for receiving back the tags "
                             "a2oTbot_proj, a2oUbot_proj, a2oVbot_proj on ocn pes" );
        }
        std::cout << "defined tag agian on ocn\n";
        // send the tag to ocean pes, from ocean mesh on coupler pes
        //   from couComm, using common joint comm ocn_coupler
        // as always, use nonblocking sends
        // original graph (context is -1_
        if( couComm != MPI_COMM_NULL )
        {
            // need to use ocean comp id for context
            context_id = cmpocn;  // id for ocean on comp
            ierr       = iMOAB_SendElementTag( cplOcnPID, bottomProjectedFields3, &ocnCouComm, &context_id );
            CHECKIERR( ierr, "cannot send tag values back to ocean pes" )
        }
        std::cout << "sent ocn data from coupler to component\n";

        // receive on component 2, ocean
        if( ocnComm != MPI_COMM_NULL )
        {
            context_id = cplocn;  // id for ocean on coupler
            ierr       = iMOAB_ReceiveElementTag( cmpOcnPID, bottomProjectedFields3, &ocnCouComm,
                                                  &context_id );
            CHECKIERR( ierr, "cannot receive tag values from ocean mesh on coupler pes" )
        }
        std::cout << "received ocn data from coupler to component\n";

        MPI_Barrier( MPI_COMM_WORLD );

        if( couComm != MPI_COMM_NULL )
        {
            context_id = cmpocn;
            ierr       = iMOAB_FreeSenderBuffers( cplOcnPID, &context_id );
            CHECKIERR( ierr, "cannot free send/receive buffers for OCN context" )
        }
        std::cout << "freed send/recv ocn data from coupler to component\n";
        if( ocnComm != MPI_COMM_NULL && 1 == n )  // write only for n==1 case
        {
            char outputFileOcn[] = "Ocn2WithProj.h5m";
            ierr                 = iMOAB_WriteMesh( cmpOcnPID, outputFileOcn, fileWriteOptions );
            CHECKIERR( ierr, "could not write Ocn2WithProj.h5m to disk" )
            // test results only for n == 1, for bottomTempProjectedField
            if( !no_regression_test )
            {
                // the same as remap test
                // get temp field on ocean, from conservative, the global ids, and dump to the baseline file
                // first get GlobalIds from ocn, and fields:
                int nverts[3], nelem[3];
                ierr = iMOAB_GetMeshInfo( cmpOcnPID, nverts, nelem, 0, 0, 0 );
                CHECKIERR( ierr, "failed to get ocn mesh info" );
                std::vector< int > gidElems;
                gidElems.resize( nelem[2] );
                std::vector< double > tempElems;
                tempElems.resize( nelem[2] );
                // get global id storage
                const std::string GidStr = "GLOBAL_ID";  // hard coded too
                int tag_type = DENSE_INTEGER, ncomp = 1, tagInd = 0;
                ierr = iMOAB_DefineTagStorage( cmpOcnPID, GidStr.c_str(), &tag_type, &ncomp, &tagInd );
                CHECKIERR( ierr, "failed to define global id tag" );

                int ent_type = 1;
                ierr         = iMOAB_GetIntTagStorage( cmpOcnPID, GidStr.c_str(), &nelem[2], &ent_type, &gidElems[0] );
                CHECKIERR( ierr, "failed to get global ids" );
                ierr = iMOAB_GetDoubleTagStorage( cmpOcnPID, "a2oTbot_proj", &nelem[2], &ent_type, &tempElems[0] );
                CHECKIERR( ierr, "failed to get temperature field" );
                int err_code = 1;
                check_baseline_file( baseline, gidElems, tempElems, 1.e-9, err_code );
                if( 0 == err_code )
                    std::cout << " passed baseline test atm2ocn on ocean task " << rankInOcnComm << "\n";
            }

            std::cout << "wrote ocn data on component to disk\n";
        }
#endif  // ENABLE_ATMCPLOCN_COUPLING

    }  // end loop iterations n
#ifdef ENABLE_ATMCPLOCN_COUPLING
    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplAtm2OcnPID );
        CHECKIERR( ierr, "cannot deregister app intx AO" )
    }
#endif  // ENABLE_ATMCPLOCN_COUPLING

#ifdef ENABLE_ATMOCN_COUPLING
    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplAtmOcnPID );
        CHECKIERR( ierr, "cannot deregister app intx AO" )
    }
    if( ocnComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cmpOcnPID );
        CHECKIERR( ierr, "cannot deregister app OCN1" )
    }
#endif  // ENABLE_ATMOCN_COUPLING

    if( atmComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cmpAtmPID );
        CHECKIERR( ierr, "cannot deregister app ATM1" )
    }

#ifdef ENABLE_ATMOCN_COUPLING
    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplOcnPID );
        CHECKIERR( ierr, "cannot deregister app OCNX" )
    }
#endif  // ENABLE_ATMOCN_COUPLING

    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplAtmPID );
        CHECKIERR( ierr, "cannot deregister app ATMX" )
    }

    //#endif
    ierr = iMOAB_Finalize();
    CHECKIERR( ierr, "did not finalize iMOAB" )

    // free atm coupler group and comm
    if( MPI_COMM_NULL != atmCouComm ) MPI_Comm_free( &atmCouComm );
    MPI_Group_free( &joinAtmCouGroup );
    if( MPI_COMM_NULL != atmComm ) MPI_Comm_free( &atmComm );

#ifdef ENABLE_ATMOCN_COUPLING
    if( MPI_COMM_NULL != ocnComm ) MPI_Comm_free( &ocnComm );
    // free ocn - coupler group and comm
    if( MPI_COMM_NULL != ocnCouComm ) MPI_Comm_free( &ocnCouComm );
    MPI_Group_free( &joinOcnCouGroup );
#endif

    if( MPI_COMM_NULL != couComm ) MPI_Comm_free( &couComm );

    MPI_Group_free( &atmPEGroup );
#ifdef ENABLE_ATMOCN_COUPLING
    MPI_Group_free( &ocnPEGroup );
#endif
    MPI_Group_free( &couPEGroup );
    MPI_Group_free( &jgroup );

    MPI_Finalize();

    return 0;
}
