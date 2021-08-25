/**
 * imoab_apg2_ol_coupler.cpp
 *
 *  This imoab_apg2_ol_coupler test will simulate coupling between 3 components
 *  meshes will be loaded from 3 files (atm phys + atm pg2, ocean, and land), and
 *  Atm and ocn will be migrated to coupler pes and compute intx on them
 *  then, intx will be performed between migrated meshes
 * and weights will be generated, such that a field from phys atm will be transferred to ocn
 * currently, the phys atm will send some data to be projected to ocean
 *  similar for land, but land will be intersected too with the atm pg2
 *  Land is defined by the domain mesh
 *
 * first, intersect atm and ocn, and recompute comm graph 1 between atm phys and atm_cx, for ocn intx
 * repeat the same for land; for atm/lnd coupling will use similar Fv -Fv maps; maybe later will identify some
 * equivalent to bilinear maps
 */

#include "moab/Core.hpp"
#ifndef MOAB_HAVE_MPI
#error imoab coupler test requires MPI configuration
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

// #define VERBOSE

#ifndef MOAB_HAVE_TEMPESTREMAP
#error The climate coupler test example requires MOAB configuration with TempestRemap
#endif

#define ENABLE_ATMOCN_COUPLING
// for land coupling with phys atm, we do not have to "compute" intersection
// it is enough to know that the ids are the same, we can project from atm to land that way, using a computed graph
#define ENABLE_ATMLND_COUPLING

int main( int argc, char* argv[] )
{
    int ierr;
    int rankInGlobalComm, numProcesses;
    MPI_Group jgroup;
    std::string readopts( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS" );
    std::string readoptsPhysAtm( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION" );

    /*if (argc < 2)
      return 0; // no test; we will force naming the land domain file*/
    // Timer data
    moab::CpuTimer timer;
    double timer_ops;
    std::string opName;

    int repartitioner_scheme = 0;
#ifdef MOAB_HAVE_ZOLTAN
    repartitioner_scheme = 2;  // use the geometry partitioner in that case
#endif

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rankInGlobalComm );
    MPI_Comm_size( MPI_COMM_WORLD, &numProcesses );

    MPI_Comm_group( MPI_COMM_WORLD, &jgroup );  // all processes in jgroup

    std::string atmFilename =
        "../../sandbox/MeshFiles/e3sm/ne4pg2_o240/ne4pg2_p8.h5m";  // we should use only mesh from here
    std::string atmPhysMesh =
        "../../sandbox/MeshFiles/e3sm/ne4pg2_o240/AtmPhys_pg2.h5m";  // it has some data associated to vertices, T_ph,
                                                                     // u_ph, v_ph
    // we will eventually project that data to ocean mesh, after intx atm/ocn

    // on a regular case,  5 ATM, 6 CPLATM (ATMX), 17 OCN     , 18 CPLOCN (OCNX)  ;
    // intx atm/ocn is not in e3sm yet, give a number
    //   6 * 100 + 18 = 618 : atmocnid
    //  18 * 100 + 6  = 1806 : ocnatmid on coupler pes!
    // 9 LND, 10 CPLLND
    //   6 * 100 + 10 = 610  atmlndid:
    //   10 * 100 + 6 = 1006  lndatmid: on coupler pes
    // cmpatm is for atm on atm pes
    // cmpocn is for ocean, on ocean pe
    // cplatm is for atm on coupler pes
    // cplocn is for ocean on coupler pes
    // atmocnid is for intx atm / ocn on coupler pes
    //
    int rankInAtmComm = -1;
    int cmpatm = 5, cplatm = 6;  // component ids are unique over all pes, and established in advance;
    int cmpPhysAtm = 105;        // different from atm spectral ?
#ifdef ENABLE_ATMOCN_COUPLING
    std::string ocnFilename = TestDir + "/recMeshOcn.h5m";
    int rankInOcnComm       = -1;
    int cmpocn = 17, cplocn = 18, atmocnid = 618,
        ocnatmid = 1806;  // component ids are unique over all pes, and established in advance;
#endif
#ifdef ENABLE_ATMLND_COUPLING
    std::string lndFilename = "../../sandbox/MeshFiles/e3sm/ne4pg2_o240/land_p8.h5m";
    int rankInLndComm       = -1;
    int cpllnd = 10, cmplnd = 9, atmlndid = 610,
        lndatmid = 1006;  // component ids are unique over all pes, and established in advance;
#endif

    int rankInCouComm = -1;

    int nghlay = 0;  // number of ghost layers for loading the file
    std::vector< int > groupTasks;
    int startG1 = 0, startG2 = 0, endG1 = numProcesses - 1, endG2 = numProcesses - 1, startG3 = startG1, endG3 = endG1;
    int startG4 = startG1, endG4 = endG1;  // these are for coupler layout
    int context_id = -1;                   // used now for freeing buffers

    // default: load atm on 2 proc, ocean on 2, land on 2; migrate to 2 procs, then compute intx
    // later, we need to compute weight matrix with tempestremap

    ProgOptions opts;
    opts.addOpt< std::string >( "atmosphere,t", "atm mesh filename (source)", &atmFilename );
#ifdef ENABLE_ATMOCN_COUPLING
    opts.addOpt< std::string >( "ocean,m", "ocean mesh filename (target)", &ocnFilename );
#endif

#ifdef ENABLE_ATMLND_COUPLING
    opts.addOpt< std::string >( "land,l", "land mesh filename (target)", &lndFilename );
#endif

    opts.addOpt< std::string >( "physgrid,q", "physics grid file", &atmPhysMesh );

    opts.addOpt< int >( "startAtm,a", "start task for atmosphere layout", &startG1 );
    opts.addOpt< int >( "endAtm,b", "end task for atmosphere layout", &endG1 );
#ifdef ENABLE_ATMOCN_COUPLING
    opts.addOpt< int >( "startOcn,c", "start task for ocean layout", &startG2 );
    opts.addOpt< int >( "endOcn,d", "end task for ocean layout", &endG2 );
#endif

#ifdef ENABLE_ATMLND_COUPLING
    opts.addOpt< int >( "startLnd,e", "start task for land layout", &startG3 );
    opts.addOpt< int >( "endLnd,f", "end task for land layout", &endG3 );
#endif

    opts.addOpt< int >( "startCoupler,g", "start task for coupler layout", &startG4 );
    opts.addOpt< int >( "endCoupler,j", "end task for coupler layout", &endG4 );

    opts.addOpt< int >( "partitioning,p", "partitioning option for migration", &repartitioner_scheme );

    opts.parseCommandLine( argc, argv );

    char fileWriteOptions[] = "PARALLEL=WRITE_PART";

    if( !rankInGlobalComm )
    {
        std::cout << " atm file: " << atmFilename << "\n   on tasks : " << startG1 << ":" << endG1 <<
#ifdef ENABLE_ATMOCN_COUPLING
            "\n ocn file: " << ocnFilename << "\n     on tasks : " << startG2 << ":" << endG2 <<
#endif
#ifdef ENABLE_ATMLND_COUPLING
            "\n land file: " << lndFilename << "\n     on tasks : " << startG3 << ":" << endG3 <<
#endif

            "\n atm phys file: " << atmPhysMesh << "\n     on tasks : " << startG1 << ":" << endG1 <<

            "\n  partitioning (0 trivial, 1 graph, 2 geometry) " << repartitioner_scheme << "\n  ";
    }
    // load files on 3 different communicators, groups
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

#ifdef ENABLE_ATMLND_COUPLING
    MPI_Group lndPEGroup;
    MPI_Comm lndComm;
    ierr = create_group_and_comm( startG3, endG3, jgroup, &lndPEGroup, &lndComm );
    CHECKIERR( ierr, "Cannot create lnd MPI group and communicator " )
#endif

    // we will always have a coupler
    MPI_Group couPEGroup;
    MPI_Comm couComm;
    ierr = create_group_and_comm( startG4, endG4, jgroup, &couPEGroup, &couComm );
    CHECKIERR( ierr, "Cannot create cpl MPI group and communicator " )

    // now, create the joint communicators atm_coupler, ocn_coupler, lnd_coupler
    // for each, we will have to create the group first, then the communicator

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

#ifdef ENABLE_ATMLND_COUPLING
    // lnd_coupler
    MPI_Group joinLndCouGroup;
    MPI_Comm lndCouComm;
    ierr = create_joint_comm_group( lndPEGroup, couPEGroup, &joinLndCouGroup, &lndCouComm );
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
    iMOAB_AppID cmpOcnPID = &cmpOcnAppID;                            // ocn
    int cplOcnAppID = -1, cplAtmOcnAppID = -1, cplOcnAtmAppID = -1;  // -1 means it is not initialized
    iMOAB_AppID cplOcnPID    = &cplOcnAppID;                         // ocn on coupler PEs
    iMOAB_AppID cplAtmOcnPID = &cplAtmOcnAppID;                      // intx atm - ocn on coupler PEs
    iMOAB_AppID cplOcnAtmPID = &cplOcnAtmAppID;                      // intx ocn - atm on coupler PEs

#endif

#ifdef ENABLE_ATMLND_COUPLING
    int cmpLndAppID       = -1;
    iMOAB_AppID cmpLndPID = &cmpLndAppID;                            // lnd
    int cplLndAppID = -1, cplAtmLndAppID = -1, cplLndAtmAppID = -1;  // -1 means it is not initialized
    iMOAB_AppID cplLndPID    = &cplLndAppID;                         // land on coupler PEs
    iMOAB_AppID cplAtmLndPID = &cplAtmLndAppID;  // intx atm - lnd on coupler PEs will be similar to ocn atm intx

    iMOAB_AppID cplLndAtmPID = &cplLndAtmAppID;
#endif

    int cmpPhysAtmID        = -1;
    iMOAB_AppID cmpPhAtmPID = &cmpPhysAtmID;  // phys atm; we do not need to move it to cpl!

    if( couComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( couComm, &rankInCouComm );
        // Register all the applications on the coupler PEs
        ierr = iMOAB_RegisterApplication( "ATMX", &couComm, &cplatm, cplAtmPID );  // atm on coupler pes
        CHECKIERR( ierr, "Cannot register ATM over coupler PEs" )
#ifdef ENABLE_ATMOCN_COUPLING
        ierr = iMOAB_RegisterApplication( "OCNX", &couComm, &cplocn, cplOcnPID );  // ocn on coupler pes
        CHECKIERR( ierr, "Cannot register OCN over coupler PEs" )
#endif

#ifdef ENABLE_ATMLND_COUPLING
        ierr = iMOAB_RegisterApplication( "LNDX", &couComm, &cpllnd, cplLndPID );  // lnd on coupler pes
        CHECKIERR( ierr, "Cannot register LND over coupler PEs" )
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

    MPI_Barrier( MPI_COMM_WORLD );

#ifdef ENABLE_ATMOCN_COUPLING
    // ocean
    ierr =
        setup_component_coupler_meshes( cmpOcnPID, cmpocn, cplOcnPID, cplocn, &ocnComm, &ocnPEGroup, &couComm,
                                        &couPEGroup, &ocnCouComm, ocnFilename, readopts, nghlay, repartitioner_scheme );
    CHECKIERR( ierr, "Cannot load and migrate ocn mesh" )

    if( couComm != MPI_COMM_NULL )
    {
        char outputFileTgt3[] = "recvOcn3.h5m";
        PUSH_TIMER( "Write migrated OCN mesh on coupler PEs" )
        ierr = iMOAB_WriteMesh( cplOcnPID, outputFileTgt3, fileWriteOptions, strlen( outputFileTgt3 ),
                                strlen( fileWriteOptions ) );
        CHECKIERR( ierr, "cannot write ocn mesh after receiving" )
        POP_TIMER( couComm, rankInCouComm )
    }
#endif  // #ifdef ENABLE_ATMOCN_COUPLING

    MPI_Barrier( MPI_COMM_WORLD );
    // load phys atm mesh, with some data on it already
    if( atmComm != MPI_COMM_NULL )
    {

        ierr = iMOAB_RegisterApplication( "PhysAtm", &atmComm, &cmpPhysAtm, cmpPhAtmPID );
        CHECKIERR( ierr, "Cannot register Phys Atm App " )

        // load the next component mesh
        PUSH_TIMER( "Load Phys Atm  mesh" )
        ierr = iMOAB_LoadMesh( cmpPhAtmPID, atmPhysMesh.c_str(), readoptsPhysAtm.c_str(), &nghlay, atmPhysMesh.length(),
                               readoptsPhysAtm.length() );
        CHECKIERR( ierr, "Cannot load Atm Phys  mesh on atm pes" )
        POP_TIMER( atmComm, rankInAtmComm )

        int nverts[3], nelem[3];
        ierr = iMOAB_GetMeshInfo( cmpPhAtmPID, nverts, nelem, 0, 0, 0 );
        CHECKIERR( ierr, "failed to get mesh info" );
        printf( "Phys Atm Component Mesh: %d vertices and %d elements\n", nverts[0], nelem[0] );
    }

    MPI_Barrier( MPI_COMM_WORLD );

#ifdef ENABLE_ATMLND_COUPLING
    // land
    if( lndComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_RegisterApplication( "LND1", &lndComm, &cmplnd, cmpLndPID );
        CHECKIERR( ierr, "Cannot register LND App " )
    }
    ierr = setup_component_coupler_meshes( cmpLndPID, cmplnd, cplLndPID, cpllnd, &lndComm, &lndPEGroup, &couComm,
                                           &couPEGroup, &lndCouComm, lndFilename, readoptsPhysAtm, nghlay,
                                           repartitioner_scheme );

    if( couComm != MPI_COMM_NULL )
    {  // write only for n==1 case
        char outputFileLnd[] = "recvLnd.h5m";
        ierr                 = iMOAB_WriteMesh( cplLndPID, outputFileLnd, fileWriteOptions, strlen( outputFileLnd ),
                                strlen( fileWriteOptions ) );
        CHECKIERR( ierr, "cannot write lnd mesh after receiving" )
    }

#endif  // #ifdef ENABLE_ATMLND_COUPLING

#ifdef ENABLE_ATMOCN_COUPLING
    if( couComm != MPI_COMM_NULL )
    {
        // now compute intersection between OCNx and ATMx on coupler PEs
        ierr = iMOAB_RegisterApplication( "ATMOCN", &couComm, &atmocnid, cplAtmOcnPID );
        CHECKIERR( ierr, "Cannot register ocn_atm intx over coupler pes " )

        ierr = iMOAB_RegisterApplication( "OCNATM", &couComm, &ocnatmid, cplOcnAtmPID );
        CHECKIERR( ierr, "Cannot register atm_ocn intx over coupler pes " )
    }
#endif

#ifdef ENABLE_ATMLND_COUPLING
    if( couComm != MPI_COMM_NULL )
    {
        // now compute intersection between LNDx and ATMx on coupler PEs
        ierr = iMOAB_RegisterApplication( "ATMLND", &couComm, &atmlndid, cplAtmLndPID );
        CHECKIERR( ierr, "Cannot register atm_lnd intx over coupler pes " )
        // now compute intersection between ATMx and LNDx on coupler PEs
        ierr = iMOAB_RegisterApplication( "LNDATM", &couComm, &lndatmid, cplLndAtmPID );
        CHECKIERR( ierr, "Cannot register lnd_atm intx over coupler pes " )
    }
#endif

    const char* weights_identifiers[2] = { "scalar", "scalar-pc" };
    int disc_orders[3]                 = { 4, 1, 1 };
    const char* disc_methods[3]        = { "cgll", "fv", "pcloud" };
    const char* dof_tag_names[3]       = { "GLOBAL_DOFS", "GLOBAL_ID", "GLOBAL_ID" };
#ifdef ENABLE_ATMOCN_COUPLING
    if( couComm != MPI_COMM_NULL )
    {
        PUSH_TIMER( "Compute ATM-OCN mesh intersection" )
        ierr = iMOAB_ComputeMeshIntersectionOnSphere(
            cplAtmPID, cplOcnPID, cplAtmOcnPID );  // coverage mesh was computed here, for cplAtmPID, atm on coupler pes
        // basically, atm was redistributed according to target (ocean) partition, to "cover" the ocean partitions
        // check if intx valid, write some h5m intx file
        CHECKIERR( ierr, "cannot compute intersection" )
        POP_TIMER( couComm, rankInCouComm )

        PUSH_TIMER( "Compute OCN-ATM mesh intersection" )
        ierr =
            iMOAB_ComputeMeshIntersectionOnSphere( cplOcnPID, cplAtmPID, cplOcnAtmPID );  // coverage mesh was computed
        CHECKIERR( ierr, "cannot compute intersection" )
        POP_TIMER( couComm, rankInCouComm )
    }

    if( atmCouComm != MPI_COMM_NULL )
    {
        // the new graph will be for sending data from atm comp to coverage mesh;
        // it involves initial atm app; cmpAtmPID; also migrate atm mesh on coupler pes, cplAtmPID
        // results are in cplAtmOcnPID, intx mesh; remapper also has some info about coverage mesh
        // after this, the sending of tags from atm pes to coupler pes will use the new par comm graph, that has more
        // precise info about what to send for ocean cover ; every time, we will
        //  use the element global id, which should uniquely identify the element
        PUSH_TIMER( "Compute OCN coverage graph for ATM mesh" )
        ierr = iMOAB_CoverageGraph( &atmCouComm, cmpAtmPID, cplAtmPID, cplAtmOcnPID, &cmpatm, &cplatm,
                                    &cplocn );  // it happens over joint communicator
        CHECKIERR( ierr, "cannot recompute direct coverage graph for ocean" )
        POP_TIMER( atmCouComm, rankInAtmComm )  // hijack this rank
    }
    if( ocnCouComm != MPI_COMM_NULL )
    {
        // now for the second intersection, ocn-atm; will be sending data from ocean to atm
        // Can we reuse the intx atm-ocn? Not sure yet; we will compute everything again :(
        // the new graph will be for sending data from ocn comp to coverage mesh over atm;
        // it involves initial ocn app; cmpOcnPID; also migrated ocn mesh on coupler pes, cplOcnPID
        // results are in cplOcnAtmPID, intx mesh; remapper also has some info about coverage mesh
        // after this, the sending of tags from ocn pes to coupler pes will use the new par comm graph, that has more
        // precise info about what to send for atm cover ; every time, we will
        //  use the element global id, which should uniquely identify the element
        PUSH_TIMER( "Compute ATM coverage graph for OCN mesh" )
        ierr = iMOAB_CoverageGraph( &ocnCouComm, cmpOcnPID, cplOcnPID, cplOcnAtmPID, &cmpocn, &cplocn,
                                    &cplatm );  // it happens over joint communicator, ocean + coupler
        CHECKIERR( ierr, "cannot recompute direct coverage graph for atm" )
        POP_TIMER( ocnCouComm, rankInOcnComm )  // hijack this rank
    }

    // need to compute graph between phys atm and atm/ocn intx coverage
    if( atmCouComm != MPI_COMM_NULL )
    {
        int typeA = 2;  // point cloud, phys mesh
        int typeB = 3;  // cells of atmosphere, dof based; maybe need another type for ParCommGraph graphtype ?
        ierr = iMOAB_ComputeCommGraph( cmpPhAtmPID, cplAtmOcnPID, &atmCouComm, &atmPEGroup, &couPEGroup, &typeA, &typeB,
                                       &cmpatm, &atmocnid );
        CHECKIERR( ierr, "cannot compute graph between phys grid on atm and intx between FV atm and ocn" )
    }

    // also
    // need to compute graph between ocn/atm intx and phys atm mesh
    /*if( atmCouComm != MPI_COMM_NULL )
    {
        int typeA = 3;  // cells of atmosphere, dof based; maybe need another type for ParCommGraph graphtype ?
        int typeB = 2;  // point cloud, phys mesh
        ierr = iMOAB_ComputeCommGraph( cplOcnAtmPID, cmpPhAtmPID, &atmCouComm, &couPEGroup, &atmPEGroup, &typeA, &typeB,
                                       &ocnatmid, &cmpatm );
    }*/
    // need to compute graph between atm on coupler and phys atm mesh on component
    if( atmCouComm != MPI_COMM_NULL )
    {
        int typeA = 2;  // point cloud, phys mesh
        int typeB = 3;  // cells of atmosphere, dof based; need another type for ParCommGraph graphtype ?
        ierr = iMOAB_ComputeCommGraph( cmpPhAtmPID, cplAtmPID, &atmCouComm, &atmPEGroup, &couPEGroup, &typeA, &typeB,
                                       &cmpatm, &cplatm );
        CHECKIERR( ierr, "cannot compute graph between phys grid on atm and FV atm on coupler" )
    }
#endif

#ifdef ENABLE_ATMLND_COUPLING
    if( couComm != MPI_COMM_NULL )
    {
        PUSH_TIMER( "Compute ATM-LND mesh intersection" )
        ierr = iMOAB_ComputeMeshIntersectionOnSphere( cplAtmPID, cplLndPID, cplAtmLndPID );
        CHECKIERR( ierr, "failed to compute atm - land intx for mapping" );
        POP_TIMER( couComm, rankInCouComm )

        PUSH_TIMER( "Compute LND-ATM mesh intersection" )
        ierr =
            iMOAB_ComputeMeshIntersectionOnSphere( cplLndPID, cplAtmPID, cplLndAtmPID );  // coverage mesh was computed
        CHECKIERR( ierr, "cannot compute intersection" )
        POP_TIMER( couComm, rankInCouComm )
    }
    if( atmCouComm != MPI_COMM_NULL )
    {
        // the new graph will be for sending data from atm comp to coverage mesh for land mesh;
        // it involves initial atm app; cmpAtmPID; also migrate atm mesh on coupler pes, cplAtmPID
        // results are in cplAtmLndPID, intx mesh; remapper also has some info about coverage mesh
        // after this, the sending of tags from atm pes to coupler pes will use the new par comm graph, that has more
        // precise info about what to send (specifically for land cover); every time,
        /// we will use the element global id, which should uniquely identify the element
        PUSH_TIMER( "Compute LND coverage graph for ATM mesh" )
        ierr = iMOAB_CoverageGraph( &atmCouComm, cmpAtmPID, cplAtmPID, cplAtmLndPID, &cmpatm, &cplatm,
                                    &cpllnd );  // it happens over joint communicator
        CHECKIERR( ierr, "cannot recompute direct coverage graph for land" )
        POP_TIMER( atmCouComm, rankInAtmComm )  // hijack this rank
    }

    // we will compute comm graph between atm phys and land, directly;
    // on the new TriGrid workflow, we do need intersection between atm and land;
    if( atmCouComm != MPI_COMM_NULL )
    {
        int typeA = 2;  // point cloud
        int typeB = 3;  // type 3 for land on coupler, based on global ids for land cells ?
        ierr = iMOAB_ComputeCommGraph( cmpPhAtmPID, cplAtmLndPID, &atmCouComm, &atmPEGroup, &couPEGroup, &typeA, &typeB,
                                       &cmpatm, &atmlndid );
        CHECKIERR( ierr, "cannot compute comm graph between atm and atm/lnd intersection" )
    }

    // for reverse direction, lnd - atm intx:
    if( lndCouComm != MPI_COMM_NULL )
    {
        // now for the second intersection, lnd-atm; will be sending data from lnd to atm
        // Can we reuse the intx atm-lnd? Not sure yet; we will compute everything again :(
        // the new graph will be for sending data from lnd comp to coverage mesh over atm;
        // it involves initial lnd app; cmpLndPID; also migrated lnd mesh on coupler pes, cplLndPID
        // results are in cplLndAtmPID, intx mesh; remapper also has some info about coverage mesh
        // after this, the sending of tags from lnd pes to coupler pes will use the new par comm graph, that has more
        // precise info about what to send for atm cover ; every time, we will
        //  use the element global id, which should uniquely identify the element
        PUSH_TIMER( "Compute ATM coverage graph for LND mesh" )
        ierr = iMOAB_CoverageGraph( &lndCouComm, cmpLndPID, cplLndPID, cplLndAtmPID, &cmplnd, &cpllnd,
                                    &cplatm );  // it happens over joint communicator, ocean + coupler
        CHECKIERR( ierr, "cannot recompute direct coverage graph for atm for intx with land" )
        POP_TIMER( lndCouComm, rankInLndComm )
    }
#endif
    MPI_Barrier( MPI_COMM_WORLD );

    int fMonotoneTypeID = 0, fVolumetric = 0, fValidate = 1, fNoConserve = 0, fNoBubble = 1;

#ifdef ENABLE_ATMOCN_COUPLING
#ifdef VERBOSE
    if( couComm != MPI_COMM_NULL )
    {
        char serialWriteOptions[] = "";  // for writing in serial
        std::stringstream outf;
        outf << "intxAtmOcn_" << rankInCouComm << ".h5m";
        std::string intxfile = outf.str();  // write in serial the intx file, for debugging
        ierr = iMOAB_WriteMesh( cplAtmOcnPID, (char*)intxfile.c_str(), serialWriteOptions, (int)intxfile.length(),
                                strlen( serialWriteOptions ) );
        CHECKIERR( ierr, "cannot write intx file result" )
    }
#endif

    if( couComm != MPI_COMM_NULL )
    {
        PUSH_TIMER( "Compute the projection weights with TempestRemap" )
        ierr = iMOAB_ComputeScalarProjectionWeights(
            cplAtmOcnPID, weights_identifiers[0], disc_methods[1], &disc_orders[1],  // fv
            disc_methods[1], &disc_orders[1],                                        // fv
            &fNoBubble, &fMonotoneTypeID, &fVolumetric, &fNoConserve, &fValidate, dof_tag_names[1], dof_tag_names[1],
            strlen( weights_identifiers[0] ), strlen( disc_methods[1] ), strlen( disc_methods[1] ),
            strlen( dof_tag_names[1] ), strlen( dof_tag_names[1] ) );
        CHECKIERR( ierr, "cannot compute scalar projection weights" )
        POP_TIMER( couComm, rankInCouComm )
    }

    // now compute weight maps for ocn to atm mapping;
    if( couComm != MPI_COMM_NULL )
    {
        PUSH_TIMER( "Compute the projection weights with TempestRemap for ocn - atm map" )
        ierr = iMOAB_ComputeScalarProjectionWeights(
            cplOcnAtmPID, weights_identifiers[0], disc_methods[1], &disc_orders[1],  // fv
            disc_methods[1], &disc_orders[1],                                        // fv
            &fNoBubble, &fMonotoneTypeID, &fVolumetric, &fNoConserve, &fValidate, dof_tag_names[1], dof_tag_names[1],
            strlen( weights_identifiers[0] ), strlen( disc_methods[1] ), strlen( disc_methods[1] ),
            strlen( dof_tag_names[1] ), strlen( dof_tag_names[1] ) );
        CHECKIERR( ierr, "cannot compute scalar projection weights" )
        POP_TIMER( couComm, rankInCouComm )
    }

#endif

    MPI_Barrier( MPI_COMM_WORLD );

    // we do not need to compute weights ? just send tags between dofs ?

#ifdef ENABLE_ATMLND_COUPLING
    if( couComm != MPI_COMM_NULL )
    {
        // Compute the weights to project the solution from ATM component to LND component
        PUSH_TIMER( "Compute ATM-LND remapping weights" )
        ierr = iMOAB_ComputeScalarProjectionWeights(
            cplAtmLndPID, weights_identifiers[0], disc_methods[1], &disc_orders[1], disc_methods[1], &disc_orders[1],
            &fNoBubble, &fMonotoneTypeID, &fVolumetric, &fNoConserve, &fValidate, dof_tag_names[1], dof_tag_names[1],
            strlen( weights_identifiers[0] ), strlen( disc_methods[1] ), strlen( disc_methods[1] ),
            strlen( dof_tag_names[1] ), strlen( dof_tag_names[1] ) );
        CHECKIERR( ierr, "failed to compute remapping projection weights for ATM-LND scalar non-conservative field" );
        POP_TIMER( couComm, rankInCouComm )

        // Compute the weights to project the solution from LND component to ATM component
        PUSH_TIMER( "Compute LND-ATM remapping weights" )
        ierr = iMOAB_ComputeScalarProjectionWeights(
            cplLndAtmPID, weights_identifiers[0], disc_methods[1], &disc_orders[1], disc_methods[1], &disc_orders[1],
            &fNoBubble, &fMonotoneTypeID, &fVolumetric, &fNoConserve, &fValidate, dof_tag_names[1], dof_tag_names[1],
            strlen( weights_identifiers[0] ), strlen( disc_methods[1] ), strlen( disc_methods[1] ),
            strlen( dof_tag_names[1] ), strlen( dof_tag_names[1] ) );
        CHECKIERR( ierr, "failed to compute remapping projection weights for LND-ATM scalar non-conservative field" );
        POP_TIMER( couComm, rankInCouComm )
    }
#endif

    int tagIndex[2];
    int tagTypes[2]  = { DENSE_DOUBLE, DENSE_DOUBLE };
    int atmCompNDoFs = 1 /* FV disc_orders[0]*disc_orders[0] */, ocnCompNDoFs = 1 /*FV*/;

    const char* bottomTempField          = "T_ph";  // same as on phys atm mesh
    const char* bottomTempProjectedField = "T_proj";
    // Define more fields
    const char* bottomUVelField          = "u_ph";
    const char* bottomUVelProjectedField = "u_proj";
    const char* bottomVVelField          = "v_ph";
    const char* bottomVVelProjectedField = "v_proj";

    // coming from ocn to atm, project back from T_proj
    const char* bottomTempField2 = "T2_ph";  // same as on phys atm mesh
    const char* bottomUVelField2 = "u2_ph";
    const char* bottomVVelField2 = "v2_ph";

    // coming from lnd to atm, project back from T_proj
    const char* bottomTempField3 = "T3_ph";  // same as on phys atm mesh
    const char* bottomUVelField3 = "u3_ph";
    const char* bottomVVelField3 = "v3_ph";

    // tags on phys grid atm mesh
    const char* bottomTempPhProjectedField = "Tph_proj";
    const char* bottomUVelPhProjectedField = "uph_proj";
    const char* bottomVVelPhProjectedField = "vph_proj";

    // tags on phys grid atm mesh
    const char* bottomTempPhLndProjectedField = "TphL_proj";  // L and Lnd signify the fields projected from land
    const char* bottomUVelPhLndProjectedField = "uphL_proj";
    const char* bottomVVelPhLndProjectedField = "vphL_proj";

    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DefineTagStorage( cplAtmPID, bottomTempField, &tagTypes[0], &atmCompNDoFs, &tagIndex[0],
                                       strlen( bottomTempField ) );
        CHECKIERR( ierr, "failed to define the field tag T_ph" );
#ifdef ENABLE_ATMOCN_COUPLING

        ierr = iMOAB_DefineTagStorage( cplOcnPID, bottomTempProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndex[1],
                                       strlen( bottomTempProjectedField ) );
        CHECKIERR( ierr, "failed to define the field tag T_proj" );
#endif

        ierr = iMOAB_DefineTagStorage( cplAtmPID, bottomUVelField, &tagTypes[0], &atmCompNDoFs, &tagIndex[0],
                                       strlen( bottomUVelField ) );
        CHECKIERR( ierr, "failed to define the field tag u_ph" );
#ifdef ENABLE_ATMOCN_COUPLING

        ierr = iMOAB_DefineTagStorage( cplOcnPID, bottomUVelProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndex[1],
                                       strlen( bottomUVelProjectedField ) );
        CHECKIERR( ierr, "failed to define the field tag u_proj" );
#endif

        ierr = iMOAB_DefineTagStorage( cplAtmPID, bottomVVelField, &tagTypes[0], &atmCompNDoFs, &tagIndex[0],
                                       strlen( bottomVVelField ) );
        CHECKIERR( ierr, "failed to define the field tag v_ph" );
#ifdef ENABLE_ATMOCN_COUPLING
        ierr = iMOAB_DefineTagStorage( cplOcnPID, bottomVVelProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndex[1],
                                       strlen( bottomVVelProjectedField ) );
        CHECKIERR( ierr, "failed to define the field tag v_proj" );
#endif

#ifdef ENABLE_ATMLND_COUPLING
        // need to define tag storage for land; will use the same T_proj, u_proj, v_proj name, because it will be
        // used to send between point clouds !
        // use the same ndof and same size as ocnCompNDoFs (1) !!
        ierr = iMOAB_DefineTagStorage( cplLndPID, bottomTempProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndex[1],
                                       strlen( bottomTempProjectedField ) );
        CHECKIERR( ierr, "failed to define the field tag T_proj" );
        ierr = iMOAB_DefineTagStorage( cplLndPID, bottomUVelProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndex[1],
                                       strlen( bottomUVelProjectedField ) );
        CHECKIERR( ierr, "failed to define the field tag u_proj" );
        ierr = iMOAB_DefineTagStorage( cplLndPID, bottomVVelProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndex[1],
                                       strlen( bottomVVelProjectedField ) );
        CHECKIERR( ierr, "failed to define the field tag v_proj" );
#endif
    }
    // need to make sure that the coverage mesh (created during intx method) received the tag that need to be projected
    // to target so far, the coverage mesh has only the ids and global dofs; need to change the migrate method to
    // accommodate any GLL tag now send a tag from original atmosphere (cmpAtmPID) towards migrated coverage mesh
    // (cplAtmPID), using the new coverage graph communicator

    // make the tag 0, to check we are actually sending needed data
    {
        if( cplAtmAppID >= 0 )
        {
            int nverts[3], nelem[3], nblocks[3], nsbc[3], ndbc[3];
            /*
             * Each process in the communicator will have access to a local mesh instance, which will contain the
             * original cells in the local partition and ghost entities. Number of vertices, primary cells, visible
             * blocks, number of sidesets and nodesets boundary conditions will be returned in numProcesses 3 arrays,
             * for local, ghost and total numbers.
             */
            ierr = iMOAB_GetMeshInfo( cplAtmPID, nverts, nelem, nblocks, nsbc, ndbc );
            CHECKIERR( ierr, "failed to get num primary elems" );
            int numAllElem = nelem[2];
            std::vector< double > vals;
            int storLeng = atmCompNDoFs * numAllElem;
            vals.resize( storLeng );
            for( int k = 0; k < storLeng; k++ )
                vals[k] = 0.;
            int eetype = 1;
            ierr       = iMOAB_SetDoubleTagStorage( cplAtmPID, bottomTempField, &storLeng, &eetype, &vals[0],
                                              strlen( bottomTempField ) );
            CHECKIERR( ierr, "cannot make tag nul" )
            ierr = iMOAB_SetDoubleTagStorage( cplAtmPID, bottomUVelField, &storLeng, &eetype, &vals[0],
                                              strlen( bottomUVelField ) );
            CHECKIERR( ierr, "cannot make tag nul" )
            ierr = iMOAB_SetDoubleTagStorage( cplAtmPID, bottomVVelField, &storLeng, &eetype, &vals[0],
                                              strlen( bottomVVelField ) );
            CHECKIERR( ierr, "cannot make tag nul" )
            // set the tag to 0
        }
    }

#ifdef ENABLE_ATMOCN_COUPLING
    PUSH_TIMER( "Send/receive data from atm component to coupler in ocn context" )
    if( atmComm != MPI_COMM_NULL )
    {
        // as always, use nonblocking sends
        // this is for projection to ocean:
        ierr =
            iMOAB_SendElementTag( cmpPhAtmPID, "T_ph;u_ph;v_ph;", &atmCouComm, &atmocnid, strlen( "T_ph;u_ph;v_ph;" ) );
        CHECKIERR( ierr, "cannot send tag values" )
    }
    if( couComm != MPI_COMM_NULL )
    {
        // receive on atm on coupler pes, that was redistributed according to coverage
        ierr = iMOAB_ReceiveElementTag( cplAtmOcnPID, "T_ph;u_ph;v_ph;", &atmCouComm, &cmpatm,
                                        strlen( "T_ph;u_ph;v_ph;" ) );
        CHECKIERR( ierr, "cannot receive tag values" )
    }
    POP_TIMER( MPI_COMM_WORLD, rankInGlobalComm )

    // we can now free the sender buffers
    if( atmComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_FreeSenderBuffers( cmpPhAtmPID, &atmocnid );  // context is for ocean
        CHECKIERR( ierr, "cannot free buffers used to resend atm tag towards the coverage mesh" )
    }
    /*
    #ifdef VERBOSE
      if (couComm != MPI_COMM_NULL) {
        char outputFileRecvd[] = "recvAtmCoupOcn.h5m";
        ierr = iMOAB_WriteMesh(cplAtmPID, outputFileRecvd, fileWriteOptions,
            strlen(outputFileRecvd), strlen(fileWriteOptions) );
      }
    #endif
    */

    if( couComm != MPI_COMM_NULL )
    {
        const char* concat_fieldname  = "T_ph;u_ph;v_ph;";
        const char* concat_fieldnameT = "T_proj;u_proj;v_proj;";

        /* We have the remapping weights now. Let us apply the weights onto the tag we defined
           on the source mesh and get the projection on the target mesh */
        PUSH_TIMER( "Apply Scalar projection weights" )
        ierr = iMOAB_ApplyScalarProjectionWeights( cplAtmOcnPID, weights_identifiers[0], concat_fieldname,
                                                   concat_fieldnameT, strlen( weights_identifiers[0] ),
                                                   strlen( concat_fieldname ), strlen( concat_fieldnameT ) );
        CHECKIERR( ierr, "failed to compute projection weight application" );
        POP_TIMER( couComm, rankInCouComm )

        char outputFileTgt[] = "fOcnOnCpl3.h5m";
        ierr                 = iMOAB_WriteMesh( cplOcnPID, outputFileTgt, fileWriteOptions, strlen( outputFileTgt ),
                                strlen( fileWriteOptions ) );
        CHECKIERR( ierr, "failed to write fOcnOnCpl3.h5m " );
    }
    // send the projected tag back to ocean pes, with send/receive tag
    if( ocnComm != MPI_COMM_NULL )
    {
        int tagIndexIn2;
        ierr = iMOAB_DefineTagStorage( cmpOcnPID, bottomTempProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndexIn2,
                                       strlen( bottomTempProjectedField ) );
        CHECKIERR( ierr, "failed to define the field tag for receiving back the tag T_proj on ocn pes" );
        ierr = iMOAB_DefineTagStorage( cmpOcnPID, bottomUVelProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndexIn2,
                                       strlen( bottomUVelProjectedField ) );
        CHECKIERR( ierr, "failed to define the field tag for receiving back the tag u_proj on ocn pes" );
        ierr = iMOAB_DefineTagStorage( cmpOcnPID, bottomVVelProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndexIn2,
                                       strlen( bottomVVelProjectedField ) );
        CHECKIERR( ierr, "failed to define the field tag for receiving back the tag v_proj on ocn pes" );
    }
    // send the tag to ocean pes, from ocean mesh on coupler pes
    //   from couComm, using common joint comm ocn_coupler
    // as always, use nonblocking sends
    // original graph (context is -1_
    if( couComm != MPI_COMM_NULL )
    {
        context_id = cmpocn;
        ierr       = iMOAB_SendElementTag( cplOcnPID, "T_proj;u_proj;v_proj;", &ocnCouComm, &context_id,
                                     strlen( "T_proj;u_proj;v_proj;" ) );
        CHECKIERR( ierr, "cannot send tag values back to ocean pes" )
    }

    // receive on component 2, ocean
    if( ocnComm != MPI_COMM_NULL )
    {
        context_id = cplocn;
        ierr       = iMOAB_ReceiveElementTag( cmpOcnPID, "T_proj;u_proj;v_proj;", &ocnCouComm, &context_id,
                                        strlen( "T_proj;u_proj;v_proj;" ) );
        CHECKIERR( ierr, "cannot receive tag values from ocean mesh on coupler pes" )
    }

    MPI_Barrier( MPI_COMM_WORLD );

    if( couComm != MPI_COMM_NULL )
    {
        context_id = cmpocn;
        ierr       = iMOAB_FreeSenderBuffers( cplOcnPID, &context_id );
        CHECKIERR( ierr, "cannot free buffers related to send tag" )
    }
    if( ocnComm != MPI_COMM_NULL )
    {
        char outputFileOcn[] = "OcnWithProj3.h5m";
        ierr                 = iMOAB_WriteMesh( cmpOcnPID, outputFileOcn, fileWriteOptions, strlen( outputFileOcn ),
                                strlen( fileWriteOptions ) );
        CHECKIERR( ierr, "cannot write OcnWithProj3.h5m" )
    }
#endif

    MPI_Barrier( MPI_COMM_WORLD );

#ifdef ENABLE_ATMLND_COUPLING
    // start land proj:

    // we used this to compute
    //  ierr = iMOAB_ComputeCommGraph(cmpPhAtmPID, cplLndPID, &atmCouComm, &atmPEGroup, &couPEGroup,
    //     &typeA, &typeB, &cmpatm, &atmlndid);

    // end copy
    PUSH_TIMER( "Send/receive data from phys comp atm to coupler land, using computed graph" )
    if( atmComm != MPI_COMM_NULL )
    {

        // as always, use nonblocking sends
        // this is for projection to land:
        ierr =
            iMOAB_SendElementTag( cmpPhAtmPID, "T_ph;u_ph;v_ph;", &atmCouComm, &atmlndid, strlen( "T_ph;u_ph;v_ph;" ) );
        CHECKIERR( ierr, "cannot send tag values towards cpl on land" )
    }
    if( couComm != MPI_COMM_NULL )
    {
        // receive on lnd on coupler pes
        ierr = iMOAB_ReceiveElementTag( cplAtmLndPID, "T_ph;u_ph;v_ph;", &atmCouComm, &cmpatm,
                                        strlen( "T_ph;u_ph;v_ph;" ) );
        CHECKIERR( ierr, "cannot receive tag values on land on coupler, for atm coupling" )
    }
    POP_TIMER( MPI_COMM_WORLD, rankInGlobalComm )

    // we can now free the sender buffers
    if( atmComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_FreeSenderBuffers( cmpPhAtmPID, &atmlndid );
        CHECKIERR( ierr, "cannot free buffers used to resend atm tag towards the land on coupler" )
    }

    if( couComm != MPI_COMM_NULL )
    {
        const char* concat_fieldname  = "T_ph;u_ph;v_ph;";
        const char* concat_fieldnameT = "T_proj;u_proj;v_proj;";

        /* We have the remapping weights now. Let us apply the weights onto the tag we defined
           on the source mesh and get the projection on the target mesh */
        PUSH_TIMER( "Apply Scalar projection weights" )
        ierr = iMOAB_ApplyScalarProjectionWeights( cplAtmLndPID, weights_identifiers[0], concat_fieldname,
                                                   concat_fieldnameT, strlen( weights_identifiers[0] ),
                                                   strlen( concat_fieldname ), strlen( concat_fieldnameT ) );
        CHECKIERR( ierr, "failed to compute projection weight application" );
        POP_TIMER( couComm, rankInCouComm )

        char outputFileTgt[] = "fLndOnCpl3.h5m";
        ierr                 = iMOAB_WriteMesh( cplLndPID, outputFileTgt, fileWriteOptions, strlen( outputFileTgt ),
                                strlen( fileWriteOptions ) );
        CHECKIERR( ierr, "cannot write land on coupler" )
    }

    // end land proj
    // send the tags back to land pes, from land mesh on coupler pes
    // send from cplLndPID to cmpLndPID, using common joint comm
    // as always, use nonblocking sends
    // original graph
    // int context_id = -1;
    // the land might not have these tags yet; it should be a different name for land
    // in e3sm we do have different names
    if( lndComm != MPI_COMM_NULL )
    {
        int tagIndexIn2;
        ierr = iMOAB_DefineTagStorage( cmpLndPID, bottomTempProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndexIn2,
                                       strlen( bottomTempProjectedField ) );
        CHECKIERR( ierr, "failed to define the field tag for receiving back the tag a2oTbot_proj on lnd pes" );
        ierr = iMOAB_DefineTagStorage( cmpLndPID, bottomUVelProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndexIn2,
                                       strlen( bottomUVelProjectedField ) );
        CHECKIERR( ierr, "failed to define the field tag for receiving back the tag a2oUbot_proj on lnd pes" );
        ierr = iMOAB_DefineTagStorage( cmpLndPID, bottomVVelProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndexIn2,
                                       strlen( bottomVVelProjectedField ) );
        CHECKIERR( ierr, "failed to define the field tag for receiving back the tag a2oVbot_proj on lnd pes" );
    }
    if( couComm != MPI_COMM_NULL )
    {
        context_id = cmplnd;
        ierr       = iMOAB_SendElementTag( cplLndPID, "T_proj;u_proj;v_proj;", &lndCouComm, &context_id,
                                     strlen( "T_proj;u_proj;v_proj;" ) );
        CHECKIERR( ierr, "cannot send tag values back to land pes" )
    }
    // receive on component 3, land
    if( lndComm != MPI_COMM_NULL )
    {
        context_id = cpllnd;
        ierr       = iMOAB_ReceiveElementTag( cmpLndPID, "T_proj;u_proj;v_proj;", &lndCouComm, &context_id,
                                        strlen( "T_proj;u_proj;v_proj;" ) );
        CHECKIERR( ierr, "cannot receive tag values from land mesh on coupler pes" )
    }

    MPI_Barrier( MPI_COMM_WORLD );
    if( couComm != MPI_COMM_NULL )
    {
        context_id = cmplnd;
        ierr       = iMOAB_FreeSenderBuffers( cplLndPID, &context_id );
        CHECKIERR( ierr, "cannot free buffers related to sending tags from coupler to land pes" )
    }
    if( lndComm != MPI_COMM_NULL )
    {
        char outputFileLnd[] = "LndWithProj3.h5m";
        ierr                 = iMOAB_WriteMesh( cmpLndPID, outputFileLnd, fileWriteOptions, strlen( outputFileLnd ),
                                strlen( fileWriteOptions ) );
        CHECKIERR( ierr, "cannot write land file with projection" )
    }

    // we have received on ocean component some data projected from atm, on coupler

    // path was T_ph (atm phys) towards atm on FV mesh, we projected, then we sent back to ocean comp
// start copy ocn atm coupling; go reverse direction now !!
#ifdef ENABLE_ATMOCN_COUPLING
    PUSH_TIMER( "Send/receive data from ocn component to coupler in atm context" )
    if( ocnComm != MPI_COMM_NULL )
    {
        // as always, use nonblocking sends
        // this is for projection to ocean:
        ierr = iMOAB_SendElementTag( cmpOcnPID, "T_proj;u_proj;v_proj;", &ocnCouComm, &cplatm,
                                     strlen( "T_proj;u_proj;v_proj;" ) );
        CHECKIERR( ierr, "cannot send tag values T_proj, etc towards ocn coupler" )
    }
    if( couComm != MPI_COMM_NULL )
    {
        // receive on ocn on coupler pes, that was redistributed according to coverage
        ierr = iMOAB_ReceiveElementTag( cplOcnPID, "T_proj;u_proj;v_proj;", &ocnCouComm, &cplatm,
                                        strlen( "T_proj;u_proj;v_proj;" ) );
        CHECKIERR( ierr, "cannot receive tag values" )
    }
    POP_TIMER( MPI_COMM_WORLD, rankInGlobalComm )

    // we can now free the sender buffers
    if( ocnComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_FreeSenderBuffers( cmpOcnPID, &cplatm );  // context is for atm
        CHECKIERR( ierr, "cannot free buffers used to send ocn tag towards the coverage mesh for atm" )
    }
    //#ifdef VERBOSE
    if( couComm != MPI_COMM_NULL )
    {
        // write only for n==1 case
        char outputFileRecvd[] = "recvOcnCpl.h5m";
        ierr = iMOAB_WriteMesh( cplOcnPID, outputFileRecvd, fileWriteOptions, strlen( outputFileRecvd ),
                                strlen( fileWriteOptions ) );
        CHECKIERR( ierr, "could not write recvOcnCplOcn.h5m to disk" )
    }
    //#endif

    if( couComm != MPI_COMM_NULL )
    {
        /* We have the remapping weights now. Let us apply the weights onto the tag we defined
           on the source mesh and get the projection on the target mesh */
        PUSH_TIMER( "Apply Scalar projection weights" )
        const char* concat_fieldname  = "T_proj;u_proj;v_proj;";  // this is now source tag
        const char* concat_fieldnameT = "T2_ph;u2_ph;v2_ph;";     // projected tag on
        // make sure the new tags exist on atm coupler mesh;
        ierr = iMOAB_DefineTagStorage( cplAtmPID, bottomTempField2, &tagTypes[0], &atmCompNDoFs, &tagIndex[0],
                                       strlen( bottomTempField2 ) );
        CHECKIERR( ierr, "failed to define the field tag T2_ph" );

        ierr = iMOAB_DefineTagStorage( cplAtmPID, bottomUVelField2, &tagTypes[0], &atmCompNDoFs, &tagIndex[0],
                                       strlen( bottomUVelField2 ) );
        CHECKIERR( ierr, "failed to define the field tag u2_ph" );

        ierr = iMOAB_DefineTagStorage( cplAtmPID, bottomVVelField2, &tagTypes[0], &atmCompNDoFs, &tagIndex[0],
                                       strlen( bottomVVelField2 ) );
        CHECKIERR( ierr, "failed to define the field tag v2_ph" );

        ierr = iMOAB_ApplyScalarProjectionWeights( cplOcnAtmPID, weights_identifiers[0], concat_fieldname,
                                                   concat_fieldnameT, strlen( weights_identifiers[0] ),
                                                   strlen( concat_fieldname ), strlen( concat_fieldnameT ) );
        CHECKIERR( ierr, "failed to compute projection weight application from ocn to atm " );
        POP_TIMER( couComm, rankInCouComm )

        char outputFileTgt[] = "fAtm2OnCpl2.h5m";
        ierr                 = iMOAB_WriteMesh( cplAtmPID, outputFileTgt, fileWriteOptions, strlen( outputFileTgt ),
                                strlen( fileWriteOptions ) );
        CHECKIERR( ierr, "could not write fAtm2OnCpl.h5m to disk" )
    }

    // send the projected tag back to atm pes, with send/receive partial par graph computed
    //  from intx ocn/atm towards atm physics mesh !
    if( atmComm != MPI_COMM_NULL )
    {
        int tagIndexIn2;
        ierr = iMOAB_DefineTagStorage( cmpPhAtmPID, bottomTempPhProjectedField, &tagTypes[1], &atmCompNDoFs,
                                       &tagIndexIn2, strlen( bottomTempPhProjectedField ) );
        CHECKIERR( ierr, "failed to define the field tag for receiving back the tag "
                         "Tph_proj on atm pes" );
        ierr = iMOAB_DefineTagStorage( cmpPhAtmPID, bottomUVelPhProjectedField, &tagTypes[1], &atmCompNDoFs,
                                       &tagIndexIn2, strlen( bottomUVelPhProjectedField ) );
        CHECKIERR( ierr, "failed to define the field tag for receiving back the tag "
                         "uph_proj on atm pes" );
        ierr = iMOAB_DefineTagStorage( cmpPhAtmPID, bottomVVelPhProjectedField, &tagTypes[1], &atmCompNDoFs,
                                       &tagIndexIn2, strlen( bottomVVelPhProjectedField ) );
        CHECKIERR( ierr, "failed to define the field tag for receiving back the tag "
                         "vph_proj on atm pes" );
    }
    // send the tag to atm pes, from atm pg2 mesh on coupler pes
    //   from couComm, using common joint comm atm_coupler, and computed graph, phys-grid -> atm FV on cpl, in reverse
    // as always, use nonblocking sends
    // graph context comes from commgraph ?

    //      ierr = iMOAB_ComputeCommGraph(  cmpPhAtmPID, cplAtmPID, &atmCouComm, &atmPEGroup, &couPEGroup, &typeA,
    //      &typeB,
    //  &cmpatm, &cplatm );

    if( couComm != MPI_COMM_NULL )
    {
        context_id = cmpatm;
        ierr       = iMOAB_SendElementTag( cplAtmPID, "T2_ph;u2_ph;v2_ph;", &atmCouComm, &context_id,
                                     strlen( "T2_ph;u2_ph;v2_ph;" ) );
        CHECKIERR( ierr, "cannot send tag values back to atm pes" )
    }

    // receive on component atm phys mesh
    if( atmComm != MPI_COMM_NULL )
    {
        context_id = cplatm;
        ierr       = iMOAB_ReceiveElementTag( cmpPhAtmPID, "Tph_proj;uph_proj;vph_proj;", &atmCouComm, &context_id,
                                        strlen( "Tph_proj;uph_proj;vph_proj;" ) );
        CHECKIERR( ierr, "cannot receive tag values from atm pg2 mesh on coupler pes" )
    }

    MPI_Barrier( MPI_COMM_WORLD );

    if( couComm != MPI_COMM_NULL )
    {
        context_id = cmpatm;
        ierr       = iMOAB_FreeSenderBuffers( cplAtmPID, &context_id );
        CHECKIERR( ierr, "cannot free buffers for sending T2_ph from cpl to phys atm" )
    }
    if( atmComm != MPI_COMM_NULL )  // write only for n==1 case
    {
        char outputFileAtmPh[] = "AtmPhysProj.h5m";
        ierr = iMOAB_WriteMesh( cmpPhAtmPID, outputFileAtmPh, fileWriteOptions, strlen( outputFileAtmPh ),
                                strlen( fileWriteOptions ) );
        CHECKIERR( ierr, "could not write AtmPhysProj.h5m to disk" )
    }
#endif
    // end copy need more work for land
// start copy
// lnd atm coupling; go reverse direction now, from lnd to atm , using lnd - atm intx, weight gen
// send back the T_proj , etc , from lan comp to land coupler
#ifdef ENABLE_ATMLND_COUPLING
    PUSH_TIMER( "Send/receive data from lnd component to coupler in atm context" )
    if( lndComm != MPI_COMM_NULL )
    {
        // as always, use nonblocking sends
        ierr = iMOAB_SendElementTag( cmpLndPID, "T_proj;u_proj;v_proj;", &lndCouComm, &cplatm,
                                     strlen( "T_proj;u_proj;v_proj;" ) );
        CHECKIERR( ierr, "cannot send tag values T_proj, etc towards lnd coupler" )
    }
    if( couComm != MPI_COMM_NULL )
    {
        // receive on ocn on coupler pes, that was redistributed according to coverage
        ierr = iMOAB_ReceiveElementTag( cplLndPID, "T_proj;u_proj;v_proj;", &lndCouComm, &cplatm,
                                        strlen( "T_proj;u_proj;v_proj;" ) );
        CHECKIERR( ierr, "cannot receive tag values" )
    }
    POP_TIMER( MPI_COMM_WORLD, rankInGlobalComm )

    // we can now free the sender buffers
    if( lndComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_FreeSenderBuffers( cmpLndPID, &cplatm );  // context is for atm
        CHECKIERR( ierr, "cannot free buffers used to send lnd tag towards the coverage mesh for atm" )
    }
    //#ifdef VERBOSE
    if( couComm != MPI_COMM_NULL )
    {

        char outputFileRecvd[] = "recvLndCpl.h5m";
        ierr = iMOAB_WriteMesh( cplLndPID, outputFileRecvd, fileWriteOptions, strlen( outputFileRecvd ),
                                strlen( fileWriteOptions ) );
        CHECKIERR( ierr, "could not write recvLndCpl.h5m to disk" )
    }
    //#endif

    if( couComm != MPI_COMM_NULL )
    {
        PUSH_TIMER( "Apply Scalar projection weights for lnd - atm coupling" )
        const char* concat_fieldname  = "T_proj;u_proj;v_proj;";  // this is now source tag, on land cpl
        const char* concat_fieldnameT = "T3_ph;u3_ph;v3_ph;";     // projected tag on
        // make sure the new tags exist on atm coupler mesh;
        ierr = iMOAB_DefineTagStorage( cplAtmPID, bottomTempField3, &tagTypes[0], &atmCompNDoFs, &tagIndex[0],
                                       strlen( bottomTempField3 ) );
        CHECKIERR( ierr, "failed to define the field tag T3_ph" );

        ierr = iMOAB_DefineTagStorage( cplAtmPID, bottomUVelField3, &tagTypes[0], &atmCompNDoFs, &tagIndex[0],
                                       strlen( bottomUVelField3 ) );
        CHECKIERR( ierr, "failed to define the field tag u3_ph" );

        ierr = iMOAB_DefineTagStorage( cplAtmPID, bottomVVelField3, &tagTypes[0], &atmCompNDoFs, &tagIndex[0],
                                       strlen( bottomVVelField3 ) );
        CHECKIERR( ierr, "failed to define the field tag v3_ph" );

        ierr = iMOAB_ApplyScalarProjectionWeights( cplLndAtmPID, weights_identifiers[0], concat_fieldname,
                                                   concat_fieldnameT, strlen( weights_identifiers[0] ),
                                                   strlen( concat_fieldname ), strlen( concat_fieldnameT ) );
        CHECKIERR( ierr, "failed to compute projection weight application from lnd to atm " );
        POP_TIMER( couComm, rankInCouComm )

        char outputFileTgt[] = "fAtm3OnCpl.h5m";  // this is for T3_ph, etc
        ierr                 = iMOAB_WriteMesh( cplAtmPID, outputFileTgt, fileWriteOptions, strlen( outputFileTgt ),
                                strlen( fileWriteOptions ) );
        CHECKIERR( ierr, "could not write fAtm3OnCpl.h5m to disk" )
    }

    // send the projected tag back to atm pes, with send/receive partial par graph computed
    //  from intx lnd/atm towards atm physics mesh !
    if( atmComm != MPI_COMM_NULL )
    {
        int tagIndexIn2;
        ierr = iMOAB_DefineTagStorage( cmpPhAtmPID, bottomTempPhLndProjectedField, &tagTypes[1], &atmCompNDoFs,
                                       &tagIndexIn2, strlen( bottomTempPhLndProjectedField ) );
        CHECKIERR( ierr, "failed to define the field tag for receiving back the tag "
                         "TphL_proj on atm pes" );
        ierr = iMOAB_DefineTagStorage( cmpPhAtmPID, bottomUVelPhLndProjectedField, &tagTypes[1], &atmCompNDoFs,
                                       &tagIndexIn2, strlen( bottomUVelPhLndProjectedField ) );
        CHECKIERR( ierr, "failed to define the field tag for receiving back the tag "
                         "uphL_proj on atm pes" );
        ierr = iMOAB_DefineTagStorage( cmpPhAtmPID, bottomVVelPhLndProjectedField, &tagTypes[1], &atmCompNDoFs,
                                       &tagIndexIn2, strlen( bottomVVelPhLndProjectedField ) );
        CHECKIERR( ierr, "failed to define the field tag for receiving back the tag "
                         "vphL_proj on atm pes" );
    }
    // send the tag to atm pes, from atm pg2 mesh on coupler pes
    //   from couComm, using common joint comm atm_coupler, and computed graph, phys-grid -> atm FV on cpl, in reverse
    // as always, use nonblocking sends
    // graph context comes from commgraph ?

    //      ierr = iMOAB_ComputeCommGraph(  cmpPhAtmPID, cplAtmPID, &atmCouComm, &atmPEGroup, &couPEGroup, &typeA,
    //      &typeB,
    //  &cmpatm, &cplatm );

    if( couComm != MPI_COMM_NULL )
    {
        context_id = cmpatm;
        ierr       = iMOAB_SendElementTag( cplAtmPID, "T3_ph;u3_ph;v3_ph;", &atmCouComm, &context_id,
                                     strlen( "T3_ph;u3_ph;v3_ph;" ) );
        CHECKIERR( ierr, "cannot send tag values back to atm pes" )
    }

    // receive on component atm phys mesh
    if( atmComm != MPI_COMM_NULL )
    {
        context_id = cplatm;
        ierr       = iMOAB_ReceiveElementTag( cmpPhAtmPID, "TphL_proj;uphL_proj;vphL_proj;", &atmCouComm, &context_id,
                                        strlen( "TphL_proj;uphL_proj;vphL_proj;" ) );
        CHECKIERR( ierr, "cannot receive tag values from atm pg2 mesh on coupler pes" )
    }

    MPI_Barrier( MPI_COMM_WORLD );

    if( couComm != MPI_COMM_NULL )
    {
        context_id = cmpatm;
        ierr       = iMOAB_FreeSenderBuffers( cplAtmPID, &context_id );
        CHECKIERR( ierr, "cannot free buffers used for sending back atm tags " )
    }
    if( atmComm != MPI_COMM_NULL )  // write only for n==1 case
    {
        char outputFileAtmPh[] = "AtmPhysProj3.h5m";
        ierr = iMOAB_WriteMesh( cmpPhAtmPID, outputFileAtmPh, fileWriteOptions, strlen( outputFileAtmPh ),
                                strlen( fileWriteOptions ) );
        CHECKIERR( ierr, "could not write AtmPhysProj3.h5m to disk" )
    }
#endif
    // we could deregister cplLndAtmPID
    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplLndAtmPID );
        CHECKIERR( ierr, "cannot deregister app intx LA" )
    }

    // we could deregister cplAtmLndPID
    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplAtmLndPID );
        CHECKIERR( ierr, "cannot deregister app intx AL" )
    }
#endif  // ENABLE_ATMLND_COUPLING

#ifdef ENABLE_ATMOCN_COUPLING
    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplOcnAtmPID );
        CHECKIERR( ierr, "cannot deregister app intx OA" )
    }
    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplAtmOcnPID );
        CHECKIERR( ierr, "cannot deregister app intx AO" )
    }
#endif

#ifdef ENABLE_ATMLND_COUPLING
    if( lndComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cmpLndPID );
        CHECKIERR( ierr, "cannot deregister app LND1" )
    }
#endif  // ENABLE_ATMLND_COUPLING
    if( atmComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cmpPhAtmPID );
        CHECKIERR( ierr, "cannot deregister app PhysAtm " )
    }
#ifdef ENABLE_ATMOCN_COUPLING
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

#ifdef ENABLE_ATMLND_COUPLING
    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplLndPID );
        CHECKIERR( ierr, "cannot deregister app LNDX" )
    }
#endif  // ENABLE_ATMLND_COUPLING

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

#ifdef ENABLE_ATMLND_COUPLING
    if( MPI_COMM_NULL != lndComm ) MPI_Comm_free( &lndComm );
    // free land - coupler group and comm
    if( MPI_COMM_NULL != lndCouComm ) MPI_Comm_free( &lndCouComm );
    MPI_Group_free( &joinLndCouGroup );
#endif

    if( MPI_COMM_NULL != couComm ) MPI_Comm_free( &couComm );

    MPI_Group_free( &atmPEGroup );
#ifdef ENABLE_ATMOCN_COUPLING
    MPI_Group_free( &ocnPEGroup );
#endif
#ifdef ENABLE_ATMLND_COUPLING
    MPI_Group_free( &lndPEGroup );
#endif
    MPI_Group_free( &couPEGroup );
    MPI_Group_free( &jgroup );

    MPI_Finalize();
    // endif #if 0

    return 0;
}
