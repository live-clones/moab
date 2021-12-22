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

    std::string ocnFilename = TestDir + "unittest/recMeshOcn.h5m";
    std::string mapFilename = TestDir + "unittest/atm_ocn_map.nc";  // this is a netcdf file!

    std::string baseline = TestDir + "unittest/baseline1.txt";
    int rankInOcnComm    = -1;
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

    int n = 1;  // number of send/receive / project / send back cycles
    opts.addOpt< int >( "iterations,n", "number of iterations for coupler", &n );

    bool no_regression_test = false;
    opts.addOpt< void >( "no_regression,r", "do not do regression test against baseline 1", &no_regression_test );
    opts.parseCommandLine( argc, argv );

    char fileWriteOptions[] = "PARALLEL=WRITE_PART";

    if( !rankInGlobalComm )
    {
        std::cout << " atm file: " << atmFilename << "\n   on tasks : " << startG1 << ":" << endG1
                  << "\n ocn file: " << ocnFilename << "\n     on tasks : " << startG2 << ":" << endG2
                  << "\n map file:" << mapFilename << "\n     on tasks : " << startG4 << ":" << endG4 << "\n";
        if( !no_regression_test )
        {
            std::cout << " check projection against baseline: " << baseline << "\n";
        }
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

    MPI_Group ocnPEGroup;
    MPI_Comm ocnComm;
    ierr = create_group_and_comm( startG2, endG2, jgroup, &ocnPEGroup, &ocnComm );
    CHECKIERR( ierr, "Cannot create ocn MPI group and communicator " )

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
    int cplOcnAppID = -1, cplAtmOcnAppID = -1;   // -1 means it is not initialized
    iMOAB_AppID cplOcnPID    = &cplOcnAppID;     // ocn on coupler PEs
    iMOAB_AppID cplAtmOcnPID = &cplAtmOcnAppID;  // intx atm -ocn on coupler PEs

    if( couComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( couComm, &rankInCouComm );
        // Register all the applications on the coupler PEs
        ierr = iMOAB_RegisterApplication( "ATMX", &couComm, &cplatm,
                                          cplAtmPID );  // atm on coupler pes
        CHECKIERR( ierr, "Cannot register ATM over coupler PEs" )

        ierr = iMOAB_RegisterApplication( "OCNX", &couComm, &cplocn,
                                          cplOcnPID );  // ocn on coupler pes
        CHECKIERR( ierr, "Cannot register OCN over coupler PEs" )
    }

    if( atmComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( atmComm, &rankInAtmComm );
        ierr = iMOAB_RegisterApplication( "ATM1", &atmComm, &cmpatm, cmpAtmPID );
        CHECKIERR( ierr, "Cannot register ATM App" )
        ierr = iMOAB_LoadMesh( cmpAtmPID, atmFilename.c_str(), readopts.c_str(), &nghlay );
        CHECKIERR( ierr, "Cannot load atm mesh" )
    }

    MPI_Barrier( MPI_COMM_WORLD );
    if( ocnComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( ocnComm, &rankInOcnComm );
        ierr = iMOAB_RegisterApplication( "OCN1", &ocnComm, &cmpocn, cmpOcnPID );
        CHECKIERR( ierr, "Cannot register OCN App" )
        ierr = iMOAB_LoadMesh( cmpOcnPID, ocnFilename.c_str(), readopts.c_str(), &nghlay );
        CHECKIERR( ierr, "Cannot load ocn mesh" )
    }

    MPI_Barrier( MPI_COMM_WORLD );

    if( couComm != MPI_COMM_NULL )
    {
        // now load map between OCNx and ATMx on coupler PEs
        ierr = iMOAB_RegisterApplication( "ATMOCN", &couComm, &atmocnid, cplAtmOcnPID );
        CHECKIERR( ierr, "Cannot register ocn_atm map instance over coupler pes " )
    }

    int disc_orders[3] = { 4, 1, 1 };

    const std::string intx_from_file_identifier = "map-from-file";

    if( couComm != MPI_COMM_NULL )
    {
        int dummyCpl     = -1;
        int dummy_rowcol = -1;
        int dummyType    = 0;
        ierr             = iMOAB_LoadMappingWeightsFromFile( cplAtmOcnPID, &dummyCpl, &dummy_rowcol, &dummyType,
                                                             intx_from_file_identifier.c_str(), mapFilename.c_str() );
        CHECKIERR( ierr, "failed to load map file from disk" );
    }

    if( atmCouComm != MPI_COMM_NULL )
    {
        int type      = 1;  // quads in source set
        int direction = 1;  // from source to coupler; will create a mesh on cplAtmPID
        // because it is like "coverage", context will be cplocn
        ierr = iMOAB_MigrateMapMesh( cmpAtmPID, cplAtmOcnPID, cplAtmPID, &atmCouComm, &atmPEGroup, &couPEGroup, &type,
                                     &cmpatm, &cplocn, &direction );
        CHECKIERR( ierr, "failed to migrate mesh for atm on coupler" );
#ifdef VERBOSE
        if( *cplAtmPID >= 0 )
        {
            char prefix[] = "atmcov";
            ierr          = iMOAB_WriteLocalMesh( cplAtmPID, prefix );
            CHECKIERR( ierr, "failed to write local mesh" );
        }
#endif
    }
    MPI_Barrier( MPI_COMM_WORLD );

    if( ocnCouComm != MPI_COMM_NULL )
    {
        int type      = 3;  // cells with GLOBAL_ID in ocean / target set
        int direction = 2;  // from coupler to target; will create a mesh on cplOcnPID
        // it will be like initial migrate cmpocn <-> cplocn
        ierr = iMOAB_MigrateMapMesh( cmpOcnPID, cplAtmOcnPID, cplOcnPID, &ocnCouComm, &ocnPEGroup, &couPEGroup, &type,
                                     &cmpocn, &cplocn, &direction );
        CHECKIERR( ierr, "failed to migrate mesh for ocn on coupler" );
#ifdef VERBOSE
        if( *cplOcnPID >= 0 )
        {
            char prefix[] = "ocntgt";
            ierr          = iMOAB_WriteLocalMesh( cplOcnPID, prefix );
            CHECKIERR( ierr, "failed to write local ocean mesh" );
            char outputFileRec[] = "CoupOcn.h5m";
            ierr                 = iMOAB_WriteMesh( cplOcnPID, outputFileRec, fileWriteOptions );
            CHECKIERR( ierr, "failed to write ocean global mesh file" );
        }
#endif
    }
    MPI_Barrier( MPI_COMM_WORLD );

    int tagIndex[2];
    int tagTypes[2]  = { DENSE_DOUBLE, DENSE_DOUBLE };
    int atmCompNDoFs = disc_orders[0] * disc_orders[0], ocnCompNDoFs = 1 /*FV*/;

    const char* bottomFields          = "a2oTbot:a2oUbot:a2oVbot";
    const char* bottomProjectedFields = "a2oTbot_proj:a2oUbot_proj:a2oVbot_proj";

    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DefineTagStorage( cplAtmPID, bottomFields, &tagTypes[0], &atmCompNDoFs, &tagIndex[0] );
        CHECKIERR( ierr, "failed to define the field tags a2oTbot:a2oUbot:a2oVbot" );

        ierr = iMOAB_DefineTagStorage( cplOcnPID, bottomProjectedFields, &tagTypes[1], &ocnCompNDoFs, &tagIndex[1] );
        CHECKIERR( ierr, "failed to define the field tags a2oTbot_proj:a2oUbot_proj:a2oVbot_proj " );
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
            int storLeng = atmCompNDoFs * numAllElem * 3;  // 3 tags
            int eetype   = 1;

            vals.resize( storLeng );
            for( int k = 0; k < storLeng; k++ )
                vals[k] = 0.;

            ierr = iMOAB_SetDoubleTagStorage( cplAtmPID, bottomFields, &storLeng, &eetype, &vals[0] );
            CHECKIERR( ierr, "cannot make tag nul" )
        }
    }

    const char* concat_fieldname  = "a2oTbot:a2oUbot:a2oVbot";
    const char* concat_fieldnameT = "a2oTbot_proj:a2oUbot_proj:a2oVbot_proj";

    // start a virtual loop for number of iterations
    for( int iters = 0; iters < n; iters++ )
    {

        PUSH_TIMER( MPI_COMM_WORLD, "Send/receive data from atm component to coupler in ocn context" )
        if( atmComm != MPI_COMM_NULL )
        {
            // as always, use nonblocking sends
            // this is for projection to ocean:
            ierr = iMOAB_SendElementTag( cmpAtmPID, "a2oTbot:a2oUbot:a2oVbot", &atmCouComm, &cplocn );
            CHECKIERR( ierr, "cannot send tag values" )
#ifdef GRAPH_INFO
            int is_sender = 1;
            int context   = cplocn;
            int verbose = 0;
            iMOAB_DumpCommGraph( cmpAtmPID, &context, &is_sender, &verbose, "AtmCovOcnS" );
#endif
        }
        if( couComm != MPI_COMM_NULL )
        {
            // receive on atm on coupler pes, that was redistributed according to coverage
            ierr = iMOAB_ReceiveElementTag( cplAtmPID, "a2oTbot:a2oUbot:a2oVbot", &atmCouComm, &cmpatm );
            CHECKIERR( ierr, "cannot receive tag values" )
#ifdef GRAPH_INFO
            int is_sender = 0;
            int context   = atmocnid;  // the same context
            int verbose = 0;
            iMOAB_DumpCommGraph( cmpAtmPID, &context, &is_sender, &verbose, "AtmCovOcnR" );
#endif
        }

        // we can now free the sender buffers
        if( atmComm != MPI_COMM_NULL )
        {
            ierr = iMOAB_FreeSenderBuffers( cmpAtmPID, &cplocn );  // context is for ocean
            CHECKIERR( ierr, "cannot free buffers used to resend atm tag towards the coverage mesh" )
        }
        POP_TIMER( MPI_COMM_WORLD, rankInGlobalComm )
#ifdef VERBOSE
        if( *cplAtmPID >= 0 && n == 1 )
        {
            char prefix[] = "atmcov_withdata";
            ierr          = iMOAB_WriteLocalMesh( cplAtmPID, prefix );
            CHECKIERR( ierr, "failed to write local atm cov mesh with data" );
        }

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
            PUSH_TIMER( couComm, "Apply Scalar projection weights" )
            ierr = iMOAB_ApplyScalarProjectionWeights( cplAtmOcnPID, intx_from_file_identifier.c_str(),
                                                       concat_fieldname, concat_fieldnameT );
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
                             "a2oTbot_proj:a2oUbot_proj:a2oVbot_proj on ocn pes" );
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
                iMOAB_SendElementTag( cplOcnPID, "a2oTbot_proj:a2oUbot_proj:a2oVbot_proj", &ocnCouComm, &context_id );
            CHECKIERR( ierr, "cannot send tag values back to ocean pes" )
        }

        // receive on component 2, ocean
        if( ocnComm != MPI_COMM_NULL )
        {
            context_id = cplocn;  // id for ocean on coupler
            ierr       = iMOAB_ReceiveElementTag( cmpOcnPID, "a2oTbot_proj:a2oUbot_proj:a2oVbot_proj", &ocnCouComm,
                                                  &context_id );
            CHECKIERR( ierr, "cannot receive tag values from ocean mesh on coupler pes" )
        }

        if( couComm != MPI_COMM_NULL )
        {
            context_id = cmpocn;
            ierr       = iMOAB_FreeSenderBuffers( cplOcnPID, &context_id );
        }
        MPI_Barrier( MPI_COMM_WORLD );

        if( ocnComm != MPI_COMM_NULL && 1 == n )  // write only for n==1 case
        {
#ifdef VERBOSE
            char outputFileOcn[] = "OcnWithProj.h5m";
            ierr                 = iMOAB_WriteMesh( cmpOcnPID, outputFileOcn, fileWriteOptions );
            CHECKIERR( ierr, "could not write OcnWithProj.h5m to disk" )
#endif
            // test results only for n == 1, for bottomTempProjectedField = "a2oTbot_proj"
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
        }

    }  // end loop iterations n

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

    if( atmComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cmpAtmPID );
        CHECKIERR( ierr, "cannot deregister app ATM1" )
    }

    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplOcnPID );
        CHECKIERR( ierr, "cannot deregister app OCNX" )
    }

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

    if( MPI_COMM_NULL != ocnComm ) MPI_Comm_free( &ocnComm );
    // free ocn - coupler group and comm
    if( MPI_COMM_NULL != ocnCouComm ) MPI_Comm_free( &ocnCouComm );
    MPI_Group_free( &joinOcnCouGroup );

    if( MPI_COMM_NULL != couComm ) MPI_Comm_free( &couComm );

    MPI_Group_free( &atmPEGroup );

    MPI_Group_free( &ocnPEGroup );

    MPI_Group_free( &couPEGroup );
    MPI_Group_free( &jgroup );

    MPI_Finalize();

    return 0;
}
