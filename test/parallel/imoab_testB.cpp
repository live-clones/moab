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


    int  cplatm = 6;  //
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
    int  startG2 = 0, endG2 = numProcesses - 1;
    int startG4 = startG2, endG4 = endG2;  // these are for coupler layout


    // default: load atm on 2 proc, ocean on 2, land on 2; migrate to 2 procs, then compute intx
    // later, we need to compute weight matrix with tempestremap

    ProgOptions opts;
    opts.addOpt< std::string >( "atmosphere,t", "atm mesh filename (target)", &atmFilename );
    opts.addOpt< std::string >( "ocean,m", "ocean mesh filename (source)", &ocnFilename );
    opts.addOpt< int >( "startOcn,c", "start task for ocean layout", &startG2 );
    opts.addOpt< int >( "endOcn,d", "end task for ocean layout", &endG2 );

    opts.addOpt< int >( "startCoupler,g", "start task for coupler layout", &startG4 );
    opts.addOpt< int >( "endCoupler,j", "end task for coupler layout", &endG4 );

    opts.addOpt< int >( "partitioning,p", "partitioning option for migration", &repartitioner_scheme );
    int n = 1;  // number of send/receive / project / send back cycles
    opts.addOpt< int >( "iterations,n", "number of iterations for coupler", &n );

    bool analytic_field=false;

    opts.addOpt< void >( "analytic,q", "analytic field", &analytic_field );

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
        ierr = iMOAB_LoadMesh( cplAtmPID, atmFilename.c_str(), readopts.c_str(), &nghlay );
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
        PUSH_TIMER( couComm, "Write migrated OCN mesh on coupler PEs" )
        ierr = iMOAB_WriteMesh( cplOcnPID, outputFileTgt5, fileWriteOptions );
        CHECKIERR( ierr, "cannot write ocn mesh after receiving" )
        POP_TIMER( couComm, rankInCouComm )
    }
    if( couComm != MPI_COMM_NULL )
    {
        // now compute intersection between OCNx and ATMx on coupler PEs
        ierr = iMOAB_RegisterApplication( "OCNATM", &couComm, &ocnatmid, cplOcnAtmPID );
        CHECKIERR( ierr, "Cannot register atm_ocn intx over coupler pes " )
    }

    if( couComm != MPI_COMM_NULL )
    {
        PUSH_TIMER( couComm, "Compute OCN-ATM mesh intersection" )
        ierr = iMOAB_ComputeMeshIntersectionOnSphere(
            cplOcnPID, cplAtmPID, cplOcnAtmPID );  // coverage mesh was computed here, for cplAtmPID, atm on coupler pes
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
        // results are in cplOcnAtmPID, intx mesh; remapper also has some info about coverage mesh
        // after this, the sending of tags from ocn pes to coupler pes will use the new par comm graph, that has more
        // precise info about what to send for atm cover ; every time, we will
        //  use the element global id, which should uniquely identify the element
        PUSH_TIMER( ocnCouComm, "Compute ATM coverage graph for OCN mesh" )
        ierr = iMOAB_CoverageGraph( &ocnCouComm, cmpOcnPID, cplOcnPID, cplOcnAtmPID, &cmpocn, &cplocn,
                                    &cplatm );  // it happens over joint communicator, ocean + coupler
        CHECKIERR( ierr, "cannot recompute direct coverage graph for atm" )
        POP_TIMER( ocnCouComm, rankInOcnComm )  // hijack this rank
    }

    MPI_Barrier( MPI_COMM_WORLD );

    int fMonotoneTypeID = 0, fVolumetric = 0, fValidate = 1, fNoConserve = 0, fNoBubble = 1;
    const char* weights_identifiers[2] = { "scalar", "scalar-pc" };
    int disc_orders[3]                 = { 4, 1, 1 };
    const char* disc_methods[3]        = { "cgll", "fv", "pcloud" };
    const char* dof_tag_names[3]       = { "GLOBAL_DOFS", "GLOBAL_ID", "GLOBAL_ID" };
    if( couComm != MPI_COMM_NULL )
    {
        PUSH_TIMER( couComm, "Compute the projection weights with TempestRemap" )
        ierr = iMOAB_ComputeScalarProjectionWeights(
            cplOcnAtmPID, weights_identifiers[0], disc_methods[1], &disc_orders[1],  // fv
            disc_methods[1], &disc_orders[1],                                        // fv
            &fNoBubble, &fMonotoneTypeID, &fVolumetric, &fNoConserve, &fValidate, dof_tag_names[1], dof_tag_names[1] );
        CHECKIERR( ierr, "cannot compute scalar projection weights" )
        POP_TIMER( couComm, rankInCouComm )
    }

    int tagIndex[2];
    int tagTypes[2]  = { DENSE_DOUBLE, DENSE_DOUBLE };
    int atmCompNDoFs = 1 /* FV disc_orders[0]*disc_orders[0] */, ocnCompNDoFs = 1 /*FV*/;

    const char* bottomTempField          = "T_proj";  // on ocean input files
    const char* bottomTempProjectedField = "T_proj2"; // on projected atm
    const char* bottomUVelField          = "u_proj";
    const char* bottomUVelProjectedField = "u_proj2";
    const char* bottomVVelField          = "v_proj";
    const char* bottomVVelProjectedField = "v_proj2";

    // ocn - atm coupling T_proj -> T_proj2
    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DefineTagStorage( cplOcnPID, bottomTempField, &tagTypes[0], &ocnCompNDoFs, &tagIndex[0] );
        CHECKIERR( ierr, "failed to define the field tag T_proj" );
        ierr = iMOAB_DefineTagStorage( cplAtmPID, bottomTempProjectedField, &tagTypes[1], &atmCompNDoFs, &tagIndex[1] );
        CHECKIERR( ierr, "failed to define the field tag T_proj2" );

        ierr = iMOAB_DefineTagStorage( cplOcnPID, bottomUVelField, &tagTypes[0], &ocnCompNDoFs, &tagIndex[0] );
        CHECKIERR( ierr, "failed to define the field tag u_proj" );
        ierr = iMOAB_DefineTagStorage( cplAtmPID, bottomUVelProjectedField, &tagTypes[1], &atmCompNDoFs, &tagIndex[1] );
        CHECKIERR( ierr, "failed to define the field tag u_proj2" );

        ierr = iMOAB_DefineTagStorage( cplOcnPID, bottomVVelField, &tagTypes[0], &ocnCompNDoFs, &tagIndex[0] );
        CHECKIERR( ierr, "failed to define the field tag v_proj" );
        ierr = iMOAB_DefineTagStorage( cplAtmPID, bottomVVelProjectedField, &tagTypes[1], &atmCompNDoFs, &tagIndex[1] );
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
    MPI_Barrier( MPI_COMM_WORLD );
    // start a virtual loop for number of iterations
    for( int iters = 0; iters < n; iters++ )
    {
        PUSH_TIMER( MPI_COMM_WORLD, "Send/receive data from ocn component to coupler in atm context" )
        if( ocnComm != MPI_COMM_NULL )
        {
          // as always, use nonblocking sends
          // this is for projection to atm, from ocean:
            ierr = iMOAB_SendElementTag( cmpOcnPID, "T_proj;u_proj;v_proj;", &ocnCouComm, &cplatm );
            CHECKIERR( ierr, "cannot send tag values" )
        }
        if( couComm != MPI_COMM_NULL )
        {
          // receive on ocn on coupler pes, that was redistributed according to coverage
            ierr = iMOAB_ReceiveElementTag( cplOcnPID, "T_proj;u_proj;v_proj;", &ocnCouComm, &cplatm );
            CHECKIERR( ierr, "cannot receive tag values" )
        }


        // we can now free the sender buffers
        if( ocnComm != MPI_COMM_NULL )
        {
            ierr = iMOAB_FreeSenderBuffers( cmpOcnPID, &cplatm );  // context is for atm
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
            ierr = iMOAB_ApplyScalarProjectionWeights( cplOcnAtmPID, weights_identifiers[0], concat_fieldname,
                                                    concat_fieldnameT );
            CHECKIERR( ierr, "failed to compute projection weight application" );
            POP_TIMER( couComm, rankInCouComm )
         // do not write if iters > 0)
           if( 1 == n )
           {
               char outputFileTgt[] = "fAtmOnCpl5.h5m";
               ierr = iMOAB_WriteMesh( cplAtmPID, outputFileTgt, fileWriteOptions );
               CHECKIERR( ierr, "failed to write fAtmOnCpl3.h5m " );
           }
        }
    }
        // do not need to send the tag to atm pes, from atm mesh on coupler pes; we are already on atm pes
   // free data
    // free up the MPI objects and finalize

    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplOcnAtmPID );
        CHECKIERR( ierr, "cannot deregister app OCNATM" )
    }

    if( ocnComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cmpOcnPID );
        CHECKIERR( ierr, "cannot deregister app OCN1" )
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


    ierr = iMOAB_Finalize();
    CHECKIERR( ierr, "did not finalize iMOAB" )

    // free atm coupler group and comm


    if( MPI_COMM_NULL != ocnComm ) MPI_Comm_free( &ocnComm );
    // free ocn - coupler group and comm
    if( MPI_COMM_NULL != ocnCouComm ) MPI_Comm_free( &ocnCouComm );
    MPI_Group_free( &joinOcnCouGroup );

    if( MPI_COMM_NULL != couComm ) MPI_Comm_free( &couComm );


    MPI_Group_free( &ocnPEGroup );

    MPI_Group_free( &couPEGroup );
    MPI_Group_free( &jgroup );

    MPI_Finalize();

    return 0;
}
