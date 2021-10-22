/*
 * This imoab_coupler test will simulate coupling between 3 components
 * 3 meshes will be loaded from 3 files (src, ocean, lnd), and they will be migrated to
 * all processors (coupler pes); then, intx will be performed between migrated meshes
 * and weights will be generated, such that a field from one component will be transferred to
 * the other component
 * currently, the src will send some data to be projected to ocean and land components
 *
 * first, intersect src and tgt, and recompute comm graph 1 between src and src_cx, for tgt intx
 * second, intersect src and lnd, and recompute comm graph 2 between src and src_cx for lnd intx

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

int main( int argc, char* argv[] )
{
    int ierr;
    int rankInGlobalComm, numProcesses;
    MPI_Group jgroup;
    std::string readopts( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS" );
    std::string readoptsLnd( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION" );
    const std::string srctgt_map_file_name = TestDir + "/src_tgt_map.nc";

    int fMonotoneTypeID = 0, fVolumetric = 0, fValidate = 1, fNoConserve = 0, fNoBubble = 1;

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

    std::string srcFilename = TestDir + "/wholeATM_T.h5m";
    // on a regular case,  5 ATM, 6 CPLATM (ATMX), 17 OCN     , 18 CPLOCN (OCNX)  ;
    // intx src/tgt is not in e3sm yet, give a number
    //   6 * 100+ 18 = 618 : cplsrctgtid
    // 9 LND, 10 CPLLND
    //   6 * 100 + 10 = 610  srclndid:
    // cmpsrc is for src on src pes
    // cmptgt is for ocean, on ocean pe
    // cplsrc is for src on coupler pes
    // cpltgt is for ocean on coupelr pes
    // cplsrctgtid is for intx src / tgt on coupler pes
    //
    int rankInSrcComm       = -1;
    int cmpsrc              = 5;  // component ids are unique over all pes, and established in advance;
    std::string tgtFilename = TestDir + "/recMeshOcn.h5m";
    std::string baseline    = TestDir + "/baseline1.txt";
    int rankInTgtComm       = -1;
    int cmptgt              = 17;
    int cplsrctgtid         = 618;  // component ids are unique over all pes, and established in advance;

    int rankInCouComm = -1, rankInCouCommIO = -1;

    int nghlay = 0;  // number of ghost layers for loading the file
    std::vector< int > groupTasks;
    int startG1 = 0, startG2 = 0, endG1 = numProcesses - 1,
        endG2   = numProcesses - 1;        // Support launch of imoab_coupler test on any combo of 2*x processes
    int startG4 = startG1, endG4 = endG1;  // these are for coupler layout

    // default: load src on 2 proc, ocean on 2, land on 2; migrate to 2 procs, then compute intx
    // later, we need to compute weight matrix with tempestremap

    ProgOptions opts;
    opts.addOpt< std::string >( "srcosphere,t", "src mesh filename (source)", &srcFilename );
    opts.addOpt< std::string >( "ocean,m", "ocean mesh filename (target)", &tgtFilename );

    opts.addOpt< int >( "startSrc,a", "start task for source component layout", &startG1 );
    opts.addOpt< int >( "endSrc,b", "end task for source component layout", &endG1 );
    opts.addOpt< int >( "startTgt,c", "start task for target component layout", &startG2 );
    opts.addOpt< int >( "endTgt,d", "end task for target component layout", &endG2 );

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
        std::cout << " src file: " << srcFilename << "\n   on tasks : " << startG1 << ":" << endG1
                  << "\n tgt file: " << tgtFilename << "\n     on tasks : " << startG2 << ":" << endG2
                  << "\n  partitioning (0 trivial, 1 graph, 2 geometry) " << repartitioner_scheme << "\n  ";
    }

    // load files on 3 different communicators, groups
    // first groups has task 0, second group tasks 0 and 1
    // coupler will be on joint tasks, will be on a third group (0 and 1, again)
    // first groups has task 0, second group tasks 0 and 1
    // coupler will be on joint tasks, will be on a third group (0 and 1, again)
    MPI_Group srcPEGroup;
    MPI_Comm srcComm = MPI_COMM_NULL;
    ierr             = create_group_and_comm( startG1, endG1, jgroup, &srcPEGroup, &srcComm );
    CHECKIERR( ierr, "Cannot create src MPI group and communicator " )

    MPI_Group tgtPEGroup;
    MPI_Comm tgtComm = MPI_COMM_NULL;
    ierr             = create_group_and_comm( startG2, endG2, jgroup, &tgtPEGroup, &tgtComm );
    CHECKIERR( ierr, "Cannot create tgt MPI group and communicator " )

    // we will always have a coupler
    MPI_Group couPEGroup;
    MPI_Comm couComm = MPI_COMM_NULL;
    ierr             = create_group_and_comm( startG4, endG4, jgroup, &couPEGroup, &couComm );
    CHECKIERR( ierr, "Cannot create cpl MPI group and communicator " )

    // src_coupler
    MPI_Group joinSrcCouGroup;
    MPI_Comm srcCouComm = MPI_COMM_NULL;
    ierr                = create_joint_comm_group( srcPEGroup, couPEGroup, &joinSrcCouGroup, &srcCouComm );
    CHECKIERR( ierr, "Cannot create joint src cou communicator" )

    // tgt_coupler
    MPI_Group joinTgtCouGroup;
    MPI_Comm tgtCouComm = MPI_COMM_NULL;
    ierr                = create_joint_comm_group( tgtPEGroup, couPEGroup, &joinTgtCouGroup, &tgtCouComm );
    CHECKIERR( ierr, "Cannot create joint tgt cou communicator" )

    ierr = iMOAB_Initialize( argc, argv );  // not really needed anything from argc, argv, yet; maybe we should
    CHECKIERR( ierr, "Cannot initialize iMOAB" )

    int cmpSrcAppID          = -1;
    iMOAB_AppID cmpSrcPID    = &cmpSrcAppID;  // src
    int cmpTgtAppID          = -1;
    iMOAB_AppID cmpTgtPID    = &cmpTgtAppID;     // tgt
    int cplSrcTgtAppID       = -1;               // -1 means it is not initialized
    iMOAB_AppID cplSrcTgtPID = &cplSrcTgtAppID;  // intx src -tgt on coupler PEs

    if( couComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( couComm, &rankInCouComm );
        // Register all the applications on the coupler PEs
        ierr = iMOAB_RegisterApplication( "CouplerApp", &couComm, &cplsrctgtid,
                                          cplSrcTgtPID );  // src on coupler pes
        CHECKIERR( ierr, "Cannot register application over coupler PEs" )
    }

    if( srcComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( srcComm, &rankInSrcComm );
        ierr = iMOAB_RegisterApplication( "SourceComponent", &srcComm, &cmpsrc, cmpSrcPID );
        CHECKIERR( ierr, "Cannot register ATM App" )
    }

    if( tgtComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( tgtComm, &rankInTgtComm );
        ierr = iMOAB_RegisterApplication( "TargetComponent", &tgtComm, &cmptgt, cmpTgtPID );
        CHECKIERR( ierr, "Cannot register OCN App" )
    }

    // source
    ierr = setup_component_coupler_meshes( cmpSrcPID, cmpsrc, cplSrcTgtPID, cplsrctgtid, &srcComm, &srcPEGroup,
                                           &couComm, &couPEGroup, &srcCouComm, srcFilename, readopts, nghlay,
                                           repartitioner_scheme, false );
    CHECKIERR( ierr, "Cannot load and migrate src mesh" )

    // target
    ierr = setup_component_coupler_meshes( cmpTgtPID, cmptgt, cplSrcTgtPID, cplsrctgtid, &tgtComm, &tgtPEGroup,
                                           &couComm, &couPEGroup, &tgtCouComm, tgtFilename, readopts, nghlay,
                                           repartitioner_scheme, false );
    CHECKIERR( ierr, "Cannot load and migrate tgt mesh" )

    MPI_Barrier( MPI_COMM_WORLD );

    const std::string weights_identifiers = "map-from-file";
    const std::string dof_tag_names[2]    = { "GLOBAL_ID", "GLOBAL_ID" };

    // if( couComm != MPI_COMM_NULL )
    // {
    //     PUSH_TIMER( "Compute Source-Target mesh relations for I/O server" )
    //     ierr = iMOAB_ComputePointDoFIntersection(
    //         cmpSrcPID, cmpTgtPID,
    //         cplSrcTgtPID );  // coverage mesh was computed here, for cplSrcPID, src on coupler pes
    //     CHECKIERR( ierr, "cannot setup source-target relations" )
    //     POP_TIMER( couComm, rankInCouComm )
    // }

    if( couComm != MPI_COMM_NULL || srcCouComm != MPI_COMM_NULL || tgtCouComm != MPI_COMM_NULL )
    {
#ifdef MOAB_HAVE_NETCDF
        ierr = iMOAB_LoadMappingWeightsFromFile( cmpSrcPID, cmpTgtPID, cplSrcTgtPID, weights_identifiers.c_str(),
                                                 srctgt_map_file_name.c_str(), NULL, NULL, NULL,
                                                 weights_identifiers.size(), srctgt_map_file_name.size() );
        CHECKIERR( ierr, "failed to load map file from disk" );
#else
#error "Require NetCDF dependency to run the test"
#endif
    }

    // if( srcCouComm != MPI_COMM_NULL )
    // {
    //     int typeA = 3;  // spectral mesh, with GLOBAL_DOFS tags on cells
    //     int typeB = 3;  // point cloud mesh, with GLOBAL_ID tag on vertices
    //     // the new graph will be for sending data from src comp to coverage mesh;
    //     // it involves initial src app; cmpSrcPID; also migrate src mesh on coupler pes, cplSrcPID
    //     // results are in cplSrcTgtPID, intx mesh; remapper also has some info about coverage mesh
    //     // after this, the sending of tags from src pes to coupler pes will use the new par comm
    //     // graph, that has more precise info about what to send for ocean cover ; every time, we
    //     // will
    //     //  use the element global id, which should uniquely identify the element
    //     PUSH_TIMER( "Compute communication graph for source mesh" )
    //     ierr = iMOAB_ComputeCommGraph( cmpSrcPID, cplIOSrcTgtPID, &srcCouComm, &srcPEGroup, &joinSrcCouGroup,
    //                                    &typeA, &typeB, &cmpsrc, &srctgtidio );  // it happens over joint communicator
    //     CHECKIERR( ierr, "cannot compute communication graph for source mesh" )
    //     POP_TIMER( srcCouComm, rankInSrcComm )  // hijack this rank
    // }

    // if( tgtCouComm != MPI_COMM_NULL )
    // {
    //     int typeA = 3;  // spectral mesh, with GLOBAL_DOFS tags on cells
    //     int typeB = 3;  // point cloud mesh, with GLOBAL_ID tag on vertices
    //     // the new graph will be for sending data from src comp to coverage mesh;
    //     // it involves initial src app; cmpSrcPID; also migrate src mesh on coupler pes, cplSrcPID
    //     // results are in cplSrcTgtPID, intx mesh; remapper also has some info about coverage mesh
    //     // after this, the sending of tags from src pes to coupler pes will use the new par comm
    //     // graph, that has more precise info about what to send for ocean cover ; every time, we
    //     // will
    //     //  use the element global id, which should uniquely identify the element
    //     PUSH_TIMER( "Compute communication graph for target mesh" )
    //     ierr = iMOAB_ComputeCommGraph( cmpTgtPID, cplIOSrcTgtPID, &tgtCouComm, &tgtPEGroup, &joinTgtCouGroup,
    //                                    &typeA, &typeB, &cmptgt, &srctgtidio );  // it happens over joint communicator
    //     CHECKIERR( ierr, "cannot compute communication graph for target mesh" )
    //     POP_TIMER( tgtCouComm, rankInTgtComm )  // hijack this rank
    // }

    int tagIndex[2];
    int srcCompNDoFs = 1;
    int tgtCompNDoFs = 1; /*FV*/
    int tagTypes     = DENSE_DOUBLE;

    const char* bottomTempField          = "srcsoln";
    const char* bottomTempProjectedField = "tgtsoln";

    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DefineTagStorage( cplSrcTgtPID, bottomTempField, &tagTypes, &srcCompNDoFs, &tagIndex[0],
                                       strlen( bottomTempField ) );
        CHECKIERR( ierr, "failed to define the field tag srcsoln" );

        ierr = iMOAB_DefineTagStorage( cplSrcTgtPID, bottomTempProjectedField, &tagTypes, &tgtCompNDoFs, &tagIndex[1],
                                       strlen( bottomTempProjectedField ) );
        CHECKIERR( ierr, "failed to define the field tag tgtsoln" );
    }

    // need to make sure that the coverage mesh (created during intx method) received the tag that
    // need to be projected to target so far, the coverage mesh has only the ids and global dofs;
    // need to change the migrate method to accommodate any GLL tag
    // now send a tag from original srcosphere (cmpSrcPID) towards migrated coverage mesh
    // (cplSrcPID), using the new coverage graph communicator

    // make the tag 0, to check we are actually sending needed data
    {
        if( cplSrcTgtPID >= 0 )
        {
            int nverts[3], nelem[3], nblocks[3], nsbc[3], ndbc[3];
            /*
             * Each process in the communicator will have access to a local mesh instance, which
             * will contain the original cells in the local partition and ghost entities. Number of
             * vertices, primary cells, visible blocks, number of sidesets and nodesets boundary
             * conditions will be returned in numProcesses 3 arrays, for local, ghost and total
             * numbers.
             */
            ierr = iMOAB_GetMeshInfo( cplSrcPID, nverts, nelem, nblocks, nsbc, ndbc );
            CHECKIERR( ierr, "failed to get num primary elems" );
            int numAllElem = nelem[2];
            std::vector< double > vals;
            int storLeng = srcCompNDoFs * numAllElem;
            int eetype   = 1;

            vals.resize( storLeng, -1.0 );

            ierr = iMOAB_SetDoubleTagStorage( cplSrcTgtPID, bottomTempField, &storLeng, &eetype, &vals[0],
                                              strlen( bottomTempField ) );
            CHECKIERR( ierr, "cannot reset double tag storage for source coupler" )

            ierr = iMOAB_SetDoubleTagStorage( cplSrcTgtPID, bottomTempProjectedField, &storLeng, &eetype, &vals[0],
                                              strlen( bottomTempProjectedField ) );
            CHECKIERR( ierr, "cannot reset double tag storage for source coupler" )

            // set the tag to 0
        }

        // send the projected tag back to ocean pes, with send/receive tag
        if( tgtComm != MPI_COMM_NULL )
        {
            int tagIndexIn2;
            ierr = iMOAB_DefineTagStorage( cmpTgtPID, bottomTempProjectedField, &tagTypes[1], &tgtCompNDoFs,
                                           &tagIndexIn2, strlen( bottomTempProjectedField ) );
            CHECKIERR( ierr, "failed to define the field tag for receiving back the tag "
                             "srcsoln_proj on target pes" );
        }

        // send the projected tag back to ocean pes, with send/receive tag
        if( srcComm != MPI_COMM_NULL )
        {
            int tagIndexIn;
            ierr = iMOAB_DefineTagStorage( cmpSrcPID, bottomTempField, &tagTypes[0], &srcCompNDoFs, &tagIndexIn,
                                           strlen( bottomTempField ) );
            CHECKIERR( ierr, "failed to define the field tag for receiving back the tag "
                             "srcsoln on source pes" );

            int nverts[3], nelem[3], nblocks[3], nsbc[3], ndbc[3];
            /*
             * Each process in the communicator will have access to a local mesh instance, which
             * will contain the original cells in the local partition and ghost entities. Number of
             * vertices, primary cells, visible blocks, number of sidesets and nodesets boundary
             * conditions will be returned in numProcesses 3 arrays, for local, ghost and total
             * numbers.
             */
            ierr = iMOAB_GetMeshInfo( cmpSrcPID, nverts, nelem, nblocks, nsbc, ndbc );
            CHECKIERR( ierr, "failed to get num primary elems" );

            int numAllElem = nelem[2];
            int storLeng   = srcCompNDoFs * numAllElem;
            std::vector< double > vals( storLeng, 0.0 );

            int coord_length = numAllElem * 3;
            std::vector< double > coords( coord_length );
            ierr = iMOAB_GetOwnedElementsCoordinates( cmpSrcPID, &coord_length, coords.data() );
            CHECKIERR( ierr, "failed to get element coordinates" );

            for( int index = 0; index < numAllElem; ++index )
            {
                const double* node = &coords[3 * index];
                double dLon        = atan2( node[1], node[0] );
                if( dLon < 0.0 ) { dLon += 2.0 * M_PI; }
                double dLat = asin( node[2] );
                vals[index] = sample_stationary_vortex( dLon, dLat );
                // printf( "%d = %f\n", index, vals[index] );
            }

            int eetype = 1;

            ierr = iMOAB_SetDoubleTagStorage( cmpSrcPID, bottomTempField, &storLeng, &eetype, &vals[0],
                                              strlen( bottomTempField ) );
            CHECKIERR( ierr, "cannot set double tag storage for source component" )

            // ierr = iMOAB_WriteMesh( cmpSrcPID, "fSrcOnCmp.h5m", fileWriteOptions, strlen( "fSrcOnCmp.h5m" ),
            //                         strlen( fileWriteOptions ) );
            // CHECKIERR( ierr, "could not write fTgtOnCpl.h5m to disk" )
        }
    }

    // start a virtual loop for number of iterations
    for( int iters = 0; iters < n; iters++ )
    {
        PUSH_TIMER( "Send/receive data from src component to coupler in tgt context" )
        if( srcComm != MPI_COMM_NULL )
        {
            int is_sender = 1;
            int context   = cpltgt;
            ierr          = iMOAB_DumpCommGraph( cmpSrcPID, &context, &is_sender, "SrcCovTar", strlen( "SrcCovTar" ) );
            CHECKIERR( ierr, "cannot dump communication graph" )

            int nelem[3];
            ierr = iMOAB_GetMeshInfo( cmpSrcPID, nullptr, nelem, nullptr, nullptr, nullptr );
            CHECKIERR( ierr, "failed to get mesh info from source" );

            int storLeng = nelem[0], eetype = 1;
            std::vector< double > vals_test( storLeng, 0.0 );
            ierr = iMOAB_GetDoubleTagStorage( cmpSrcPID, bottomTempField, &storLeng, &eetype, &vals_test[0],
                                              strlen( bottomTempField ) );
            CHECKIERR( ierr, "cannot get component double tag storage for source component" )

            for( int i = 0; i < vals_test.size(); ++i )
            {
                printf( "Sending %d value = %f\n", i, vals_test[i] );
            }

            // as always, use nonblocking sends
            // this is for projection to ocean:
            ierr = iMOAB_SendElementTag( cmpSrcPID, bottomTempField, &srcCouComm[0], &cplSrcTgtPID,
                                         strlen( bottomTempField ) );
            CHECKIERR( ierr, "cannot send tag values to coupler" )
        }
        if( couComm != MPI_COMM_NULL )
        {
            int is_sender = 0;
            int context   = cplSrcTgtPID;
            ierr          = iMOAB_DumpCommGraph( cplSrcTgtPID, &context, &is_sender, "TarSrcR", strlen( "TarSrcR" ) );
            cplSrcTgtPID ierr = iMOAB_GetMeshInfo( cplSrcTgtPID, nullptr, nelem, nullptr, nullptr, nullptr );

            // receive on src on coupler pes, that was redistributed according to coverage
            ierr = iMOAB_ReceiveElementTag( cplSrcTgtPID, bottomTempField, &srcCouComm[0], &cplSrcTgtPID,
                                            strlen( bottomTempField ) );
            CHECKIERR( ierr, "cannot receive tag values from coupler" )

            int nelem[3];
            ierr = iMOAB_GetMeshInfo( cplSrcTgtPID, nullptr, nelem, nullptr, nullptr, nullptr );
            CHECKIERR( ierr, "failed to get mesh info from source mesh on coupler" );

            int storLeng = nelem[0], eetype = 1;
            std::vector< double > vals_test( storLeng, 0.0 );
            ierr = iMOAB_GetDoubleTagStorage( cplSrcTgtPID, bottomTempField, &storLeng, &eetype, &vals_test[0],
                                              strlen( bottomTempField ) );
            CHECKIERR( ierr, "cannot get double tag storage for source coupler" )

            for( int i = 0; i < vals_test.size(); ++i )
            {
                printf( "Received %d value = %f\n", i, vals_test[i] );
            }

            ierr = iMOAB_WriteMesh( cplSrcTgtPID, "fSrcOnCpl.h5m", fileWriteOptions, strlen( "fSrcOnCpl.h5m" ),
                                    strlen( fileWriteOptions ) );
            CHECKIERR( ierr, "could not write fTgtOnCpl.h5m to disk" )
        }

        // we can now free the sender buffers
        if( srcComm != MPI_COMM_NULL )
        {
            ierr = iMOAB_FreeSenderBuffers( cmpSrcPID, &cplSrcTgtPID;  // context is for ocean
            CHECKIERR( ierr, "cannot free buffers used to resend src tag towards the coverage mesh" )
        }

        POP_TIMER( MPI_COMM_WORLD, rankInGlobalComm )

        MPI_Barrier( MPI_COMM_WORLD );

        if( couComm != MPI_COMM_NULL && false )
        {
            /* We have the remapping weights now. Let us apply the weights onto the tag we defined
               on the source mesh and get the projection on the target mesh */
            PUSH_TIMER( "Apply Scalar projection weights" )
            ierr = iMOAB_ApplyScalarProjectionWeights( cplSrcTgtPID, weights_identifiers[0].c_str(), bottomTempField,
                                                       bottomTempProjectedField, weights_identifiers[0].size(),
                                                       strlen( bottomTempField ), strlen( bottomTempProjectedField ) );
            CHECKIERR( ierr, "failed to compute projection weight application" );
            POP_TIMER( couComm, rankInCouComm )
            if( 1 == n )  // write only for n==1 case
            {
                char outputFileTgt[] = "fTgtOnCpl.h5m";
                ierr = iMOAB_WriteMesh( cplTgtPID, outputFileTgt, fileWriteOptions, strlen( outputFileTgt ),
                                        strlen( fileWriteOptions ) );
                CHECKIERR( ierr, "could not write fTgtOnCpl.h5m to disk" )

                ierr = iMOAB_WriteMesh( cplSrcPID, "fSrcOnCpl.h5m", fileWriteOptions, strlen( "fSrcOnCpl.h5m" ),
                                        strlen( fileWriteOptions ) );
                CHECKIERR( ierr, "could not write fTgtOnCpl.h5m to disk" )
            }
        }

        // send the tag to ocean pes, from ocean mesh on coupler pes
        //   from couComm, using common joint comm tgt_coupler
        // as always, use nonblocking sends
        // original graph (context is -1_
        if( couComm != MPI_COMM_NULL )
        {
            ierr = iMOAB_SendElementTag( cplTgtPID, bottomTempProjectedField, &tgtCouComm, &context_id,
                                         strlen( bottomTempProjectedField ) );
            CHECKIERR( ierr, "cannot send tag values back to ocean pes" )
        }

        // receive on component 2, ocean
        if( tgtComm != MPI_COMM_NULL && tgtCouComm[0] != MPI_COMM_NULL )
        {
            ierr = iMOAB_ReceiveElementTag( cmpTgtPID, bottomTempProjectedField, &tgtCouComm, &context_id,
                                            strlen( bottomTempProjectedField ) );
            CHECKIERR( ierr, "cannot receive tag values from ocean mesh on coupler pes" )
        }

        // we can now free the sender buffers
        if( tgtComm != MPI_COMM_NULL && tgtCouComm != MPI_COMM_NULL )
        {
            ierr = iMOAB_FreeSenderBuffers( cmpTgtPID, &context_id );  // context is for ocean
            CHECKIERR( ierr, "cannot free buffers used to resend src tag towards the coverage mesh" )
        }

        MPI_Barrier( MPI_COMM_WORLD );

        // if( couComm != MPI_COMM_NULL ) { ierr = iMOAB_FreeSenderBuffers( cplTgtPID, &context_id ); }
        // if( couComm != MPI_COMM_NULL ) { ierr = iMOAB_FreeSenderBuffers( cplTgtPID, &context_id ); }
        if( tgtComm != MPI_COMM_NULL && 1 == n )  // write only for n==1 case
        {
            char outputFileTgt[] = "TgtWithProj.h5m";
            ierr                 = iMOAB_WriteMesh( cmpTgtPID, outputFileTgt, fileWriteOptions, strlen( outputFileTgt ),
                                    strlen( fileWriteOptions ) );
            CHECKIERR( ierr, "could not write TgtWithProj.h5m to disk" )
            // test results only for n == 1, for bottomTempProjectedField
            if( !no_regression_test )
            {
                // the same as remap test
                // get temp field on ocean, from conservative, the global ids, and dump to the baseline file
                // first get GlobalIds from tgt, and fields:
                int nverts[3], nelem[3];
                ierr = iMOAB_GetMeshInfo( cmpTgtPID, nverts, nelem, 0, 0, 0 );
                CHECKIERR( ierr, "failed to get tgt mesh info" );
                std::vector< int > gidElems;
                gidElems.resize( nelem[2] );
                std::vector< double > tempElems;
                tempElems.resize( nelem[2] );
                // get global id storage
                const std::string GidStr = "GLOBAL_ID";  // hard coded too
                int lenG                 = (int)GidStr.length();
                int tag_type = DENSE_INTEGER, ncomp = 1, tagInd = 0;
                ierr = iMOAB_DefineTagStorage( cmpTgtPID, GidStr.c_str(), &tag_type, &ncomp, &tagInd, lenG );
                CHECKIERR( ierr, "failed to define global id tag" );

                int ent_type = 1;
                ierr = iMOAB_GetIntTagStorage( cmpTgtPID, GidStr.c_str(), &nelem[2], &ent_type, &gidElems[0], lenG );
                CHECKIERR( ierr, "failed to get global ids" );
                ierr = iMOAB_GetDoubleTagStorage( cmpTgtPID, bottomTempProjectedField, &nelem[2], &ent_type,
                                                  &tempElems[0], strlen( bottomTempProjectedField ) );
                CHECKIERR( ierr, "failed to get temperature field" );
                int err_code = 1;
                check_baseline_file( baseline, gidElems, tempElems, 1.e-9, err_code );
                if( 0 == err_code )
                    std::cout << " passed baseline test src2tgt on ocean task " << rankInTgtComm << "\n";
            }
        }

    }  // end loop iterations n

    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplSrcTgtPID );
        CHECKIERR( ierr, "cannot deregister app intx AO" )
    }
    if( tgtComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cmpTgtPID );
        CHECKIERR( ierr, "cannot deregister app OCN1" )
    }

    if( srcComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cmpSrcPID );
        CHECKIERR( ierr, "cannot deregister app ATM1" )
    }

    ierr = iMOAB_Finalize();
    CHECKIERR( ierr, "did not finalize iMOAB" )

    // free src coupler group and comm
    if( MPI_COMM_NULL != srcCouComm )
    {
        MPI_Comm_free( &srcCouComm );
        MPI_Group_free( &joinSrcCouGroup );
    }
    if( MPI_COMM_NULL != srcComm ) MPI_Comm_free( &srcComm );
    if( MPI_COMM_NULL != tgtComm ) MPI_Comm_free( &tgtComm );

    // free tgt - coupler group and comm
    if( MPI_COMM_NULL != tgtCouComm )
    {
        MPI_Comm_free( &tgtCouComm );
        MPI_Group_free( &joinTgtCouGroup );
    }
    if( MPI_COMM_NULL != couComm ) MPI_Comm_free( &couComm );

    MPI_Group_free( &srcPEGroup );
    MPI_Group_free( &tgtPEGroup );
    MPI_Group_free( &couPEGroup );
    MPI_Group_free( &jgroup );

    MPI_Finalize();

    return 0;
}
