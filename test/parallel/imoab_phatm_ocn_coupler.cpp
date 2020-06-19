/**
 * imoab_phAtm_ocn_coupler.cpp
 *
 *  This imoab_phAtm_ocn_coupler test will simulate coupling between 3 components
 *  meshes will be loaded from 3 files (atm, ocean, and phys atm), and
 *  Atm and ocn will be migrated to coupler pes and compute intx on them
 *  then, intx will be performed between migrated meshes
 * and weights will be generated, such that a field from phys atm will be transferred to ocn
 * currently, the phys atm will send some data to be projected to ocean
 *
 * first, intersect atm and ocn, and recompute comm graph 1 between atm phys and atm_cx, for ocn intx
 * maybe repeat the same for land
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

#define CHECKIERR(rc, message)  if (0!=rc) { printf ("%s. ErrorCode = %d\n", message, rc); return 1;}
#define PUSH_TIMER(operation)  { timer_ops = timer.time_since_birth(); opName = operation;}
#define POP_TIMER(localcomm,localrank) { \
  double locElapsed=timer.time_since_birth() - timer_ops, minElapsed=0, maxElapsed=0; \
  MPI_Reduce(&locElapsed, &maxElapsed, 1, MPI_DOUBLE, MPI_MAX, 0, localcomm); \
  MPI_Reduce(&locElapsed, &minElapsed, 1, MPI_DOUBLE, MPI_MIN, 0, localcomm); \
  if (!localrank) std::cout << "[LOG] Time taken to " << opName.c_str() << ": max = " << maxElapsed << ", avg = " << (maxElapsed+minElapsed)/2 << "\n"; \
  opName.clear(); \
}

using namespace moab;

//#define VERBOSE

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
  std::string readopts("PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS");
  std::string readoptsPhysAtm("PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION");

  // Timer data
  moab::CpuTimer timer;
  double timer_ops;
  std::string opName;

  int repartitioner_scheme = 0;
#ifdef MOAB_HAVE_ZOLTAN
  repartitioner_scheme = 2; // use the graph partitioner in that case
#endif

  MPI_Init(&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &rankInGlobalComm );
  MPI_Comm_size( MPI_COMM_WORLD, &numProcesses );

  MPI_Comm_group(MPI_COMM_WORLD, &jgroup);// all processes in jgroup

  std::string atmFilename = TestDir + "/wholeATM_T.h5m"; // we should use only mesh from here
  std::string atmPhysMesh = TestDir + "/AtmPhys_01.h5m"; // it has some data associated to vertices, T_ph, u_ph, v_ph
  // we will eventually project that data to ocean mesh, after intx atm/ocn

  // on a regular case,  5 ATM, 6 CPLATM (ATMX), 17 OCN     , 18 CPLOCN (OCNX)  ;
  // intx atm/ocn is not in e3sm yet, give a number
  //   6 * 100+ 18 = 618 : atmocnid
  // 9 LND, 10 CPLLND
  //   6 * 100 + 10 = 610  atmlndid:
  // cmpatm is for atm on atm pes
  // cmpocn is for ocean, on ocean pe
  // cplatm is for atm on coupler pes
  // cplocn is for ocean on coupler pes
  // atmocnid is for intx atm / ocn on coupler pes
  //
  int rankInAtmComm = -1;
  int cmpatm=5, cplatm=6;  // component ids are unique over all pes, and established in advance;
  int cmpPhysAtm = 105; // different from atm spectral ?
#ifdef ENABLE_ATMOCN_COUPLING
  std::string ocnFilename = TestDir + "/recMeshOcn.h5m";
  int rankInOcnComm = -1;
  int cmpocn=17, cplocn=18, atmocnid=618;  // component ids are unique over all pes, and established in advance;
#endif
#ifdef ENABLE_ATMLND_COUPLING
  std::string lndFilename = TestDir + "/wholeLnd.h5m";
  int rankInLndComm = -1;
  int cpllnd=10, cmplnd=9;  // component ids are unique over all pes, and established in advance;
#endif

  int rankInCouComm = -1;

  int nghlay=0; // number of ghost layers for loading the file
  std::vector<int> groupTasks;
  int startG1=0, startG2=0, endG1=numProcesses-1, endG2=numProcesses-1, startG3=startG1, endG3=endG1;
  int startG4 = startG1, endG4 = endG1; // these are for coupler layout
  int context_id = -1; // used now for freeing buffers

  // default: load atm on 2 proc, ocean on 2, land on 2; migrate to 2 procs, then compute intx
  // later, we need to compute weight matrix with tempestremap

  ProgOptions opts;
  opts.addOpt<std::string>("atmosphere,t", "atm mesh filename (source)", &atmFilename);
#ifdef ENABLE_ATMOCN_COUPLING
  opts.addOpt<std::string>("ocean,m", "ocean mesh filename (target)", &ocnFilename);
#endif

#ifdef ENABLE_ATMLND_COUPLING
  opts.addOpt<std::string>("land,l", "land mesh filename (target)", &lndFilename);
#endif

  opts.addOpt<std::string>("physgrid,q", "physics grid file", &atmPhysMesh);

  opts.addOpt<int>("startAtm,a", "start task for atmosphere layout", &startG1);
  opts.addOpt<int>("endAtm,b", "end task for atmosphere layout", &endG1);
#ifdef ENABLE_ATMOCN_COUPLING
  opts.addOpt<int>("startOcn,c", "start task for ocean layout", &startG2);
  opts.addOpt<int>("endOcn,d", "end task for ocean layout", &endG2);
#endif

#ifdef ENABLE_ATMLND_COUPLING
  opts.addOpt<int>("startLnd,e", "start task for land layout", &startG3);
  opts.addOpt<int>("endLnd,f", "end task for land layout", &endG3);
#endif


  opts.addOpt<int>("startCoupler,g", "start task for coupler layout", &startG4);
  opts.addOpt<int>("endCoupler,j", "end task for coupler layout", &endG4);

  opts.addOpt<int>("partitioning,p", "partitioning option for migration", &repartitioner_scheme);

  opts.parseCommandLine(argc, argv);

  char fileWriteOptions[] ="PARALLEL=WRITE_PART";

  if (!rankInGlobalComm)
  {
    std::cout << " atm file: " << atmFilename << "\n   on tasks : " << startG1 << ":"<<endG1 <<
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
  groupTasks.resize(numProcesses, 0);
  for (int i=startG1; i<=endG1; i++)
    groupTasks [i-startG1] = i;

  ierr = MPI_Group_incl(jgroup, endG1-startG1+1, &groupTasks[0], &atmPEGroup);
  CHECKIERR(ierr, "Cannot create atmPEGroup")

#ifdef ENABLE_ATMOCN_COUPLING
  groupTasks.clear();
  groupTasks.resize(numProcesses, 0);
  MPI_Group ocnPEGroup;
  for (int i=startG2; i<=endG2; i++)
    groupTasks [i-startG2] = i;

  ierr = MPI_Group_incl(jgroup, endG2-startG2+1, &groupTasks[0], &ocnPEGroup);
  CHECKIERR(ierr, "Cannot create ocnPEGroup")
#endif


#ifdef ENABLE_ATMLND_COUPLING
  groupTasks.clear();
  groupTasks.resize(numProcesses, 0);
  MPI_Group lndPEGroup;
  for (int i=startG3; i<=endG3; i++)
    groupTasks [i-startG3] = i;

  ierr = MPI_Group_incl(jgroup, endG3-startG3+1, &groupTasks[0], &lndPEGroup);
  CHECKIERR(ierr, "Cannot create lndPEGroup")
#endif


  // we will always have a coupler
  groupTasks.clear();
  groupTasks.resize(numProcesses, 0);
  MPI_Group couPEGroup;
  for (int i=startG4; i<=endG4; i++)
    groupTasks [i-startG4] = i;

  ierr = MPI_Group_incl(jgroup, endG4-startG4+1, &groupTasks[0], &couPEGroup);
  CHECKIERR(ierr, "Cannot create couPEGroup")

  // create 4 communicators, one for each group
  int ATM_COMM_TAG = 1;
  MPI_Comm atmComm;
  // atmComm is for atmosphere app;
  ierr = MPI_Comm_create_group(MPI_COMM_WORLD, atmPEGroup, ATM_COMM_TAG, &atmComm);
  CHECKIERR(ierr, "Cannot create atmComm")

#ifdef ENABLE_ATMOCN_COUPLING
  int OCN_COMM_TAG = 2;
  MPI_Comm ocnComm;
  // ocnComm is for ocean app
  ierr = MPI_Comm_create_group(MPI_COMM_WORLD, ocnPEGroup, OCN_COMM_TAG, &ocnComm);
  CHECKIERR(ierr, "Cannot create ocnComm")
#endif

#ifdef ENABLE_ATMLND_COUPLING
  int LND_COMM_TAG = 3;
  MPI_Comm lndComm;
  // lndComm is for land app
  ierr = MPI_Comm_create_group(MPI_COMM_WORLD, lndPEGroup, LND_COMM_TAG, &lndComm);
  CHECKIERR(ierr, "Cannot create lndComm")
#endif

  int COU_COMM_TAG = 4;
  MPI_Comm couComm;
  // lndComm is for land app
  ierr = MPI_Comm_create_group(MPI_COMM_WORLD, couPEGroup, COU_COMM_TAG, &couComm);
  CHECKIERR(ierr, "Cannot create lndComm")

  // now, create the joint communicators atm_coupler, ocn_coupler, atmphys_coupler
  // for each, we will have to create the group first, then the communicator

  // atm_coupler
  MPI_Group joinAtmCouGroup;
  ierr = MPI_Group_union(atmPEGroup, couPEGroup, &joinAtmCouGroup);
  CHECKIERR(ierr, "Cannot create joint atm cou group")
  int ATM_COU_COMM_TAG = 5;
  MPI_Comm atmCouComm;
  ierr = MPI_Comm_create_group(MPI_COMM_WORLD, joinAtmCouGroup, ATM_COU_COMM_TAG, &atmCouComm);
  CHECKIERR(ierr, "Cannot create joint atm cou communicator")

#ifdef ENABLE_ATMOCN_COUPLING
  // ocn_coupler
  MPI_Group joinOcnCouGroup;
  ierr = MPI_Group_union(ocnPEGroup, couPEGroup, &joinOcnCouGroup);
  CHECKIERR(ierr, "Cannot create joint ocn cou group")
  int OCN_COU_COMM_TAG = 6;
  MPI_Comm ocnCouComm;
  ierr = MPI_Comm_create_group(MPI_COMM_WORLD, joinOcnCouGroup, OCN_COU_COMM_TAG, &ocnCouComm);
  CHECKIERR(ierr, "Cannot create joint ocn cou communicator")
#endif


#ifdef ENABLE_ATMLND_COUPLING
  // lnd_coupler
  MPI_Group joinLndCouGroup;
  ierr = MPI_Group_union(lndPEGroup, couPEGroup, &joinLndCouGroup);
  CHECKIERR(ierr, "Cannot create joint lnd cou group")
  int LND_COU_COMM_TAG = 7;
  MPI_Comm lndCouComm;
  ierr = MPI_Comm_create_group(MPI_COMM_WORLD, joinLndCouGroup, LND_COU_COMM_TAG, &lndCouComm);
  CHECKIERR(ierr, "Cannot create joint lnd cou communicator")
#endif


  ierr = iMOAB_Initialize(argc, argv); // not really needed anything from argc, argv, yet; maybe we should
  CHECKIERR(ierr, "Cannot initialize iMOAB")

  int cmpAtmAppID =-1;
  iMOAB_AppID cmpAtmPID=&cmpAtmAppID; // atm
  int cplAtmAppID=-1;// -1 means it is not initialized
  iMOAB_AppID cplAtmPID=&cplAtmAppID; // atm on coupler PEs
#ifdef ENABLE_ATMOCN_COUPLING
  int cmpOcnAppID =-1;
  iMOAB_AppID cmpOcnPID=&cmpOcnAppID; // ocn
  int cplOcnAppID=-1, cplAtmOcnAppID=-1;// -1 means it is not initialized
  iMOAB_AppID cplOcnPID=&cplOcnAppID; // ocn on coupler PEs
  iMOAB_AppID cplAtmOcnPID=&cplAtmOcnAppID; // intx atm -ocn on coupler PEs
#endif


#ifdef ENABLE_ATMLND_COUPLING
  int cmpLndAppID =-1;
  iMOAB_AppID cmpLndPID=&cmpLndAppID; // lnd
  int cplLndAppID=-1; //, cplAtmLndAppID=-1;// -1 means it is not initialized
  iMOAB_AppID cplLndPID=&cplLndAppID; // land on coupler PEs
  //iMOAB_AppID cplAtmLndPID=&cplAtmLndAppID; // intx atm - lnd on coupler PEs not needed anymore
#endif


  int cmpPhysAtmID =-1;
  iMOAB_AppID cmpPhAtmPID=&cmpPhysAtmID; // phys atm; we do not need to move it to cpl!

  if (couComm != MPI_COMM_NULL) {
    MPI_Comm_rank( couComm, &rankInCouComm );
    // Register all the applications on the coupler PEs
    ierr = iMOAB_RegisterApplication("ATMX", &couComm, &cplatm, cplAtmPID); // atm on coupler pes
    CHECKIERR(ierr, "Cannot register ATM over coupler PEs")
#ifdef ENABLE_ATMOCN_COUPLING
    ierr = iMOAB_RegisterApplication("OCNX", &couComm, &cplocn, cplOcnPID); // ocn on coupler pes
    CHECKIERR(ierr, "Cannot register OCN over coupler PEs")
#endif

#ifdef ENABLE_ATMLND_COUPLING
    ierr = iMOAB_RegisterApplication("LNDX", &couComm, &cpllnd, cplLndPID); // lnd on coupler pes
    CHECKIERR(ierr, "Cannot register LND over coupler PEs")
#endif

  }

  if (atmComm != MPI_COMM_NULL) {
    MPI_Comm_rank( atmComm, &rankInAtmComm );
    ierr = iMOAB_RegisterApplication("ATM1", &atmComm, &cmpatm, cmpAtmPID);
    CHECKIERR(ierr, "Cannot register ATM App")

    PUSH_TIMER("Load ATM mesh")
    // load first mesh
    ierr = iMOAB_LoadMesh(cmpAtmPID, atmFilename.c_str(), readopts.c_str(), &nghlay, atmFilename.length(), readopts.length() );
    CHECKIERR(ierr, "Cannot load ATM mesh")
    POP_TIMER(atmComm, rankInAtmComm)

    PUSH_TIMER("ATM mesh: compute partition and send mesh")
    // then send mesh to coupler pes, on cplAtmPID
    // atmCouComm is a join communicator, so
    ierr = iMOAB_SendMesh(cmpAtmPID, &atmCouComm, &couPEGroup, &cplatm, &repartitioner_scheme); // send to component 3, on coupler pes
    CHECKIERR(ierr, "cannot send elements" )
    POP_TIMER(atmComm, rankInAtmComm)
  }
  // now, receive mesh, on coupler communicator; first mesh 1, atm
  if (couComm != MPI_COMM_NULL) {
    PUSH_TIMER("Receive ATM mesh and resolve shared entities")
    ierr = iMOAB_ReceiveMesh(cplAtmPID, &atmCouComm, &atmPEGroup, &cmpatm); // receive from component 1
    CHECKIERR(ierr, "cannot receive elements on ATMX app")
    POP_TIMER(couComm, rankInCouComm)
  }

  // we can now free the sender buffers
  if (atmComm != MPI_COMM_NULL) {
    ierr = iMOAB_FreeSenderBuffers(cmpAtmPID, &context_id);
    CHECKIERR(ierr, "cannot free buffers used to send atm mesh")
  }

  MPI_Barrier(MPI_COMM_WORLD);

#ifdef ENABLE_ATMOCN_COUPLING
  if (ocnComm != MPI_COMM_NULL) {
    MPI_Comm_rank( ocnComm, &rankInOcnComm );
    ierr = iMOAB_RegisterApplication("OCN1", &ocnComm, &cmpocn, cmpOcnPID);
    CHECKIERR(ierr, "Cannot register OCN App")

    // load second mesh
    PUSH_TIMER("Load OCN mesh")
    ierr = iMOAB_LoadMesh(cmpOcnPID, ocnFilename.c_str(), readopts.c_str(), &nghlay, ocnFilename.length(), readopts.length() );
    CHECKIERR(ierr, "Cannot load OCN mesh on cmpOcnPID")
    POP_TIMER(ocnComm, rankInOcnComm)

    // then send mesh to coupler pes, on cplAtmPID
    PUSH_TIMER("OCN mesh: compute partition and send mesh")
    // we could have an assert(ocnCouComm != MPI_COMM_NULL);
    ierr = iMOAB_SendMesh(cmpOcnPID, &ocnCouComm, &couPEGroup, &cplocn, &repartitioner_scheme); // send to component 3, on coupler pes
    CHECKIERR(ierr, "cannot send elements" )
    POP_TIMER(ocnComm, rankInOcnComm)
  }
  if (couComm != MPI_COMM_NULL) {
    PUSH_TIMER("Receive OCN mesh and resolve shared entities")
    // we could have an assert(ocnCouComm != MPI_COMM_NULL);
    ierr = iMOAB_ReceiveMesh(cplOcnPID, &ocnCouComm, &ocnPEGroup, &cmpocn); // receive from component 2, ocn, on coupler pes
    CHECKIERR(ierr, "cannot receive elements on OCNX app")
    POP_TIMER(couComm, rankInCouComm)
  }
  if (ocnComm != MPI_COMM_NULL) {
    ierr = iMOAB_FreeSenderBuffers(cmpOcnPID, &context_id);
    CHECKIERR(ierr, "cannot free buffers used to send ocn mesh")
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (couComm != MPI_COMM_NULL) {
    char outputFileTgt3[] = "recvOcn2.h5m";
    PUSH_TIMER("Write migrated OCN mesh on coupler PEs")
    ierr = iMOAB_WriteMesh(cplOcnPID, outputFileTgt3, fileWriteOptions,
      strlen(outputFileTgt3), strlen(fileWriteOptions) );
    CHECKIERR(ierr, "cannot write ocn mesh after receiving")
    POP_TIMER(couComm, rankInCouComm)
  }
#endif // #ifdef ENABLE_ATMOCN_COUPLING

// load phys atm mesh, with some data on it already
  if (atmComm != MPI_COMM_NULL) {

    ierr = iMOAB_RegisterApplication("PhysAtm", &atmComm, &cmpPhysAtm, cmpPhAtmPID);
    CHECKIERR(ierr, "Cannot register Phys Atm App ")

    // load the next component mesh
    PUSH_TIMER("Load Phys Atm  mesh")
    ierr = iMOAB_LoadMesh(cmpPhAtmPID, atmPhysMesh.c_str(), readoptsPhysAtm.c_str(), &nghlay, atmPhysMesh.length(), readoptsPhysAtm.length() );
    CHECKIERR(ierr, "Cannot load Atm Phys  mesh on atm pes")
    POP_TIMER(atmComm, rankInAtmComm)

    int nverts[3], nelem[3];
    ierr = iMOAB_GetMeshInfo( cmpPhAtmPID, nverts, nelem, 0, 0, 0);
    CHECKIERR(ierr, "failed to get mesh info");
    printf("Phys Atm Component Mesh: %d vertices and %d elements\n", nverts[0], nelem[0]);
  }

  MPI_Barrier(MPI_COMM_WORLD);

#ifdef ENABLE_ATMLND_COUPLING
  if (lndComm != MPI_COMM_NULL) {
    MPI_Comm_rank( lndComm, &rankInLndComm );
    ierr = iMOAB_RegisterApplication("LND1", &lndComm, &cmplnd, cmpLndPID);
    CHECKIERR(ierr, "Cannot register LND App ")

    // load the next component mesh
    PUSH_TIMER("Load LND mesh")
    ierr = iMOAB_LoadMesh(cmpLndPID, lndFilename.c_str(), readoptsPhysAtm.c_str(), &nghlay, lndFilename.length(), readoptsPhysAtm.length() );
    CHECKIERR(ierr, "Cannot load LND mesh on cplLndPID")
    POP_TIMER(lndComm, rankInLndComm)

    int nverts[3], nelem[3];
    ierr = iMOAB_GetMeshInfo( cmpLndPID, nverts, nelem, 0, 0, 0);
    CHECKIERR(ierr, "failed to get mesh info");
    printf("Land Component Mesh: %d vertices and %d elements\n", nverts[0], nelem[0]);

    // then send mesh to coupler pes, on cplLndPID
    PUSH_TIMER("LND mesh: compute partition and send mesh")
    // use joint land cou comm
    ierr = iMOAB_SendMesh(cmpLndPID, &lndCouComm, &couPEGroup, &cpllnd, &repartitioner_scheme); // send to component 3, on coupler pes
    CHECKIERR(ierr, "cannot send elements" )
    POP_TIMER(lndComm, rankInLndComm)
  }

  if (couComm != MPI_COMM_NULL) {
    PUSH_TIMER("Receive LND mesh and resolve shared entities")
    ierr = iMOAB_ReceiveMesh(cplLndPID, &lndCouComm, &lndPEGroup, &cmplnd); // receive from component 2, lnd, on coupler pes
    CHECKIERR(ierr, "cannot receive elements on LNDX app")
    POP_TIMER(couComm, rankInCouComm)
  }
  if (lndComm != MPI_COMM_NULL) {
    ierr = iMOAB_FreeSenderBuffers(cmpLndPID, &context_id);
    CHECKIERR(ierr, "cannot free buffers used to send lnd mesh")
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (couComm != MPI_COMM_NULL) {
    char outputFileLnd[] = "recvLnd2.h5m";
    PUSH_TIMER("Write migrated LND mesh on coupler PEs")
    ierr = iMOAB_WriteMesh(cplLndPID, outputFileLnd, fileWriteOptions,
      strlen(outputFileLnd), strlen(fileWriteOptions) );
    CHECKIERR(ierr, "cannot write lnd mesh after receiving")
    POP_TIMER(couComm, rankInCouComm)
  }

#endif // #ifdef ENABLE_ATMLND_COUPLING

#ifdef ENABLE_ATMOCN_COUPLING
  if (couComm != MPI_COMM_NULL) {
    // now compute intersection between OCNx and ATMx on coupler PEs
    ierr = iMOAB_RegisterApplication("ATMOCN", &couComm, &atmocnid, cplAtmOcnPID);
    CHECKIERR(ierr, "Cannot register ocn_atm intx over coupler pes ")
  }
#endif
/*
#ifdef ENABLE_ATMLND_COUPLING ! maybe we do not have to compute intx for PhysAtm - Lnd coupling to work
  if (couComm != MPI_COMM_NULL) {
    // now compute intersection between LNDx and ATMx on coupler PEs
    ierr = iMOAB_RegisterApplication("ATMLND", &couComm, &atmlndid, cplAtmLndPID);
    CHECKIERR(ierr, "Cannot register ocn_atm intx over coupler pes ")
  }
#endif
*/

  const char* weights_identifiers[2] = {"scalar", "scalar-pc"};
  int disc_orders[3] = {4, 1, 1};
  const char* disc_methods[3] = {"cgll", "fv", "pcloud"};
  const char* dof_tag_names[3] = {"GLOBAL_DOFS", "GLOBAL_ID", "GLOBAL_ID"};
#ifdef ENABLE_ATMOCN_COUPLING
  if (couComm != MPI_COMM_NULL) {
    PUSH_TIMER("Compute ATM-OCN mesh intersection")
    ierr = iMOAB_ComputeMeshIntersectionOnSphere(cplAtmPID, cplOcnPID, cplAtmOcnPID); // coverage mesh was computed here, for cplAtmPID, atm on coupler pes
    // basically, atm was redistributed according to target (ocean) partition, to "cover" the ocean partitions
    // check if intx valid, write some h5m intx file
    CHECKIERR(ierr, "cannot compute intersection" )
    POP_TIMER(couComm, rankInCouComm)
  }

  if (atmCouComm != MPI_COMM_NULL)
  {
    // the new graph will be for sending data from atm comp to coverage mesh;
    // it involves initial atm app; cmpAtmPID; also migrate atm mesh on coupler pes, cplAtmPID
    // results are in cplAtmOcnPID, intx mesh; remapper also has some info about coverage mesh
    // after this, the sending of tags from atm pes to coupler pes will use the new par comm graph, that has more precise info about
    // what to send for ocean cover ; every time, we will
    //  use the element global id, which should uniquely identify the element
    PUSH_TIMER("Compute OCN coverage graph for ATM mesh")
    ierr = iMOAB_CoverageGraph(&atmCouComm, cmpAtmPID,  cplAtmPID,  cplAtmOcnPID, &cplocn); // it happens over joint communicator
    CHECKIERR(ierr, "cannot recompute direct coverage graph for ocean" )
    POP_TIMER(atmCouComm, rankInAtmComm) // hijack this rank
  }
#endif
  // need to compute graph between phys atm and atm/ocn intx coverage
  if (atmCouComm != MPI_COMM_NULL)
  {
    int typeA = 2; // point cloud
    int typeB = 1; // quads in coverage set
    ierr = iMOAB_ComputeCommGraph(cmpPhAtmPID, cplAtmOcnPID, &atmCouComm, &atmPEGroup, &couPEGroup,
        &typeA, &typeB, &cmpatm, &atmocnid);
  }

/*

#ifdef ENABLE_ATMLND_COUPLING
  if (couComm != MPI_COMM_NULL) {
    PUSH_TIMER("Compute ATM-LND mesh intersection")
    ierr = iMOAB_ComputePointDoFIntersection(cplAtmPID, cplLndPID, cplAtmLndPID);
    CHECKIERR(ierr, "failed to compute point-cloud mapping");
    POP_TIMER(couComm, rankInCouComm)
  }
  if (atmCouComm != MPI_COMM_NULL)
  {
    // the new graph will be for sending data from atm comp to coverage mesh for land mesh;
    // it involves initial atm app; cmpAtmPID; also migrate atm mesh on coupler pes, cplAtmPID
    // results are in cplAtmLndPID, intx mesh; remapper also has some info about coverage mesh
    // after this, the sending of tags from atm pes to coupler pes will use the new par comm graph, that has more precise info about
    // what to send (specifically for land cover); every time,
    /// we will use the element global id, which should uniquely identify the element
    PUSH_TIMER("Compute LND coverage graph for ATM mesh")
    ierr = iMOAB_CoverageGraph(&atmCouComm, cmpAtmPID,  cplAtmPID,  cplAtmLndPID, &cpllnd); // it happens over joint communicator
    CHECKIERR(ierr, "cannot recompute direct coverage graph for land" )
    POP_TIMER(atmCouComm, rankInAtmComm) // hijack this rank
  }
#endif
*/

#ifdef ENABLE_ATMLND_COUPLING
  // we will compute comm graph between atm phys and land, directly; we do not need any intersection
  // later, we will send data from atm phys to land on coupler; then back to land comp;
  if (atmCouComm != MPI_COMM_NULL)
  {
    int typeA = 2; // point cloud
    int typeB = 2; // point cloud for land on coupler, too
    ierr = iMOAB_ComputeCommGraph(cmpPhAtmPID, cplLndPID, &atmCouComm, &atmPEGroup, &couPEGroup,
        &typeA, &typeB, &cmpatm, &cpllnd);
  }
#endif
  MPI_Barrier(MPI_COMM_WORLD);

  int fMonotoneTypeID=0, fVolumetric=0, fValidate=1, fNoConserve=0;

#ifdef ENABLE_ATMOCN_COUPLING
#ifdef VERBOSE
  if (couComm != MPI_COMM_NULL) {
    char serialWriteOptions[] =""; // for writing in serial
    std::stringstream outf;
    outf<<"intxAtmOcn_" << rankInCouComm<<".h5m";
    std::string intxfile=outf.str(); // write in serial the intx file, for debugging
    ierr = iMOAB_WriteMesh(cplAtmOcnPID, (char*)intxfile.c_str(), serialWriteOptions, (int)intxfile.length(), strlen(serialWriteOptions));
    CHECKIERR(ierr, "cannot write intx file result" )
  }
#endif

  if (couComm != MPI_COMM_NULL) {
    PUSH_TIMER("Compute the projection weights with TempestRemap")
    ierr = iMOAB_ComputeScalarProjectionWeights ( cplAtmOcnPID, weights_identifiers[0],
                                                  disc_methods[0], &disc_orders[0],
                                                  disc_methods[1], &disc_orders[1],
                                                  &fMonotoneTypeID, &fVolumetric, &fNoConserve, &fValidate,
                                                  dof_tag_names[0], dof_tag_names[1],
                                                  strlen(weights_identifiers[0]),
                                                  strlen(disc_methods[0]), strlen(disc_methods[1]),
                                                  strlen(dof_tag_names[0]), strlen(dof_tag_names[1]) );
    CHECKIERR(ierr, "cannot compute scalar projection weights" )
    POP_TIMER(couComm, rankInCouComm)
  }

#endif

  MPI_Barrier(MPI_COMM_WORLD);

  // we do not need to compute weights ? just send tags between dofs ?
/*
#ifdef ENABLE_ATMLND_COUPLING
  if (couComm != MPI_COMM_NULL) {
     Compute the weights to preoject the solution from ATM component to LND compoenent
    PUSH_TIMER("Compute ATM-LND remapping weights")
    ierr = iMOAB_ComputeScalarProjectionWeights ( cplAtmLndPID,
                                                weights_identifiers[1],
                                                disc_methods[0], &disc_orders[0],
                                                disc_methods[2], &disc_orders[2],
                                                &fMonotoneTypeID, &fVolumetric, &fNoConserve, &fValidate,
                                                dof_tag_names[0], dof_tag_names[2],
                                                strlen(weights_identifiers[1]),
                                                strlen(disc_methods[0]), strlen(disc_methods[2]),
                                                strlen(dof_tag_names[0]), strlen(dof_tag_names[2])
                                              );
    CHECKIERR(ierr, "failed to compute remapping projection weights for ATM-LND scalar non-conservative field");
    POP_TIMER(couComm, rankInCouComm)
  }
#endif
*/

  int tagIndex[2];
  int tagTypes[2] = { DENSE_DOUBLE, DENSE_DOUBLE } ;
  int atmCompNDoFs = disc_orders[0]*disc_orders[0], ocnCompNDoFs = 1/*FV*/;

  const char* bottomTempField = "T16_ph";
  const char* bottomTempProjectedField = "T_proj";
  // Define more fields
  const char* bottomUVelField = "u16_ph";
  const char* bottomUVelProjectedField = "u_proj";
  const char* bottomVVelField = "v16_ph";
  const char* bottomVVelProjectedField = "v_proj";

  if (couComm != MPI_COMM_NULL) {
    ierr = iMOAB_DefineTagStorage(cplAtmPID, bottomTempField, &tagTypes[0], &atmCompNDoFs, &tagIndex[0],  strlen(bottomTempField) );
    CHECKIERR(ierr, "failed to define the field tag T16_ph");
#ifdef ENABLE_ATMOCN_COUPLING

    ierr = iMOAB_DefineTagStorage(cplOcnPID, bottomTempProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndex[1],  strlen(bottomTempProjectedField) );
    CHECKIERR(ierr, "failed to define the field tag T_proj");
#endif


    ierr = iMOAB_DefineTagStorage(cplAtmPID, bottomUVelField, &tagTypes[0], &atmCompNDoFs, &tagIndex[0],  strlen(bottomUVelField) );
    CHECKIERR(ierr, "failed to define the field tag u16_ph");
#ifdef ENABLE_ATMOCN_COUPLING

    ierr = iMOAB_DefineTagStorage(cplOcnPID, bottomUVelProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndex[1],  strlen(bottomUVelProjectedField) );
    CHECKIERR(ierr, "failed to define the field tag u_proj");
#endif


    ierr = iMOAB_DefineTagStorage(cplAtmPID, bottomVVelField, &tagTypes[0], &atmCompNDoFs, &tagIndex[0],  strlen(bottomVVelField) );
    CHECKIERR(ierr, "failed to define the field tag v16_ph");
#ifdef ENABLE_ATMOCN_COUPLING
    ierr = iMOAB_DefineTagStorage(cplOcnPID, bottomVVelProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndex[1],  strlen(bottomVVelProjectedField) );
    CHECKIERR(ierr, "failed to define the field tag v_proj");
#endif

#ifdef ENABLE_ATMLND_COUPLING
    // need to define tag storage for land; will use the same T_proj, u_proj, v_proj name, because it will be
    // used to send between point clouds !
    // use the same ndof and same size as ocnCompNDoFs (1) !!
    ierr = iMOAB_DefineTagStorage(cplLndPID, bottomTempProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndex[1],  strlen(bottomTempProjectedField) );
    CHECKIERR(ierr, "failed to define the field tag a2oTbot_proj");
    ierr = iMOAB_DefineTagStorage(cplLndPID, bottomUVelProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndex[1],  strlen(bottomUVelProjectedField) );
    CHECKIERR(ierr, "failed to define the field tag a2oUbot_proj");
    ierr = iMOAB_DefineTagStorage(cplLndPID, bottomVVelProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndex[1],  strlen(bottomVVelProjectedField) );
    CHECKIERR(ierr, "failed to define the field tag a2oUbot_proj");
#endif
  }
  // need to make sure that the coverage mesh (created during intx method) received the tag that need to be projected to target
  // so far, the coverage mesh has only the ids and global dofs;
  // need to change the migrate method to accommodate any GLL tag
  // now send a tag from original atmosphere (cmpAtmPID) towards migrated coverage mesh (cplAtmPID), using the new coverage graph communicator

  // make the tag 0, to check we are actually sending needed data
  {
    if (cplAtmAppID >= 0)
    {
      int nverts[3], nelem[3], nblocks[3], nsbc[3], ndbc[3];
        /*
         * Each process in the communicator will have access to a local mesh instance, which will contain the
         * original cells in the local partition and ghost entities. Number of vertices, primary cells, visible blocks,
         * number of sidesets and nodesets boundary conditions will be returned in numProcesses 3 arrays, for local, ghost and total
         * numbers.
         */
        ierr = iMOAB_GetMeshInfo(  cplAtmPID, nverts, nelem, nblocks, nsbc, ndbc);
        CHECKIERR(ierr, "failed to get num primary elems");
        int numAllElem = nelem[2];
        std::vector<double> vals;
        int storLeng = atmCompNDoFs*numAllElem;
        vals.resize(storLeng);
        for (int k=0; k<storLeng; k++)
          vals[k] = 0.;
        int eetype = 1;
        ierr = iMOAB_SetDoubleTagStorage ( cplAtmPID, bottomTempField, &storLeng, &eetype, &vals[0], strlen(bottomTempField));
        CHECKIERR(ierr, "cannot make tag nul")
        ierr = iMOAB_SetDoubleTagStorage ( cplAtmPID, bottomUVelField, &storLeng, &eetype, &vals[0], strlen(bottomUVelField));
                CHECKIERR(ierr, "cannot make tag nul")
        ierr = iMOAB_SetDoubleTagStorage ( cplAtmPID, bottomVVelField, &storLeng, &eetype, &vals[0], strlen(bottomVVelField));
                CHECKIERR(ierr, "cannot make tag nul")
        // set the tag to 0
    }
  }


#ifdef ENABLE_ATMOCN_COUPLING
  PUSH_TIMER("Send/receive data from atm component to coupler in ocn context")
  if (atmComm != MPI_COMM_NULL ){
     // as always, use nonblocking sends
     // this is for projection to ocean:
     ierr = iMOAB_SendElementTag(cmpPhAtmPID, "T_ph;u_ph;v_ph;", &atmCouComm, &atmocnid, strlen("T_ph;u_ph;v_ph;"));
     CHECKIERR(ierr, "cannot send tag values")
#ifdef VERBOSE
    int is_sender = 1;
    int context = atmocnid; // used to identity the parcommgraph in charge
    iMOAB_DumpCommGraph(cmpPhAtmPID,  &context, &is_sender, "PhysAtmA2OS", strlen("PhysAtmA2OS"));
#endif
  }
  if (couComm != MPI_COMM_NULL) {
    // receive on atm on coupler pes, that was redistributed according to coverage
    ierr = iMOAB_ReceiveElementTag(cplAtmOcnPID, "T16_ph;u16_ph;v16_ph;", &atmCouComm, &cmpatm, strlen("T16_ph;u16_ph;v16_ph;"));
    CHECKIERR(ierr, "cannot receive tag values")
#ifdef VERBOSE
    int is_sender = 0;
    int context = cmpatm; // used to identity the parcommgraph in charge
    iMOAB_DumpCommGraph(cplAtmOcnPID,  &context, &is_sender, "PhysAtmA2OR", strlen("PhysAtmA2OR"));
#endif
  }
  POP_TIMER(MPI_COMM_WORLD, rankInGlobalComm)

  // we can now free the sender buffers
  if (atmComm != MPI_COMM_NULL) {
    ierr = iMOAB_FreeSenderBuffers(cmpPhAtmPID,  &atmocnid); // context is for ocean
    CHECKIERR(ierr, "cannot free buffers used to resend atm tag towards the coverage mesh")
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

  const char * concat_fieldname = "T16_ph;u16_ph;v16_ph;";
  const char * concat_fieldnameT = "T_proj;u_proj;v_proj;";

  if (couComm != MPI_COMM_NULL) {
    /* We have the remapping weights now. Let us apply the weights onto the tag we defined
       on the source mesh and get the projection on the target mesh */
    PUSH_TIMER("Apply Scalar projection weights")
    ierr = iMOAB_ApplyScalarProjectionWeights ( cplAtmOcnPID, weights_identifiers[0],
                                              concat_fieldname,
                                              concat_fieldnameT,
                                              strlen(weights_identifiers[0]),
                                              strlen(concat_fieldname),
                                              strlen(concat_fieldnameT)
                                              );
    CHECKIERR(ierr, "failed to compute projection weight application");
    POP_TIMER(couComm, rankInCouComm)

    char outputFileTgt[] = "fOcnOnCpl2.h5m";
    ierr = iMOAB_WriteMesh(cplOcnPID, outputFileTgt, fileWriteOptions,
      strlen(outputFileTgt), strlen(fileWriteOptions) );
  }
  // send the projected tag back to ocean pes, with send/receive tag
  if (ocnComm != MPI_COMM_NULL) {
    int tagIndexIn2;
    ierr = iMOAB_DefineTagStorage(cmpOcnPID, bottomTempProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndexIn2,  strlen(bottomTempProjectedField) );
    CHECKIERR(ierr, "failed to define the field tag for receiving back the tag T_proj on ocn pes");
    ierr = iMOAB_DefineTagStorage(cmpOcnPID, bottomUVelProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndexIn2,  strlen(bottomUVelProjectedField) );
    CHECKIERR(ierr, "failed to define the field tag for receiving back the tag u_proj on ocn pes");
    ierr = iMOAB_DefineTagStorage(cmpOcnPID, bottomVVelProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndexIn2,  strlen(bottomVVelProjectedField) );
    CHECKIERR(ierr, "failed to define the field tag for receiving back the tag v_proj on ocn pes");
  }
  // send the tag to ocean pes, from ocean mesh on coupler pes
  //   from couComm, using common joint comm ocn_coupler
      // as always, use nonblocking sends
  // original graph (context is -1_
  if (couComm != MPI_COMM_NULL){
    ierr = iMOAB_SendElementTag(cplOcnPID, "T_proj;u_proj;v_proj;", &ocnCouComm, &context_id,
        strlen("T_proj;u_proj;v_proj;"));
    CHECKIERR(ierr, "cannot send tag values back to ocean pes")
  }

  // receive on component 2, ocean
  if (ocnComm != MPI_COMM_NULL)
  {
    ierr = iMOAB_ReceiveElementTag(cmpOcnPID, "T_proj;u_proj;v_proj;",
        &ocnCouComm, &context_id, strlen("T_proj;u_proj;v_proj;"));
    CHECKIERR(ierr, "cannot receive tag values from ocean mesh on coupler pes")
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (couComm != MPI_COMM_NULL){
    ierr = iMOAB_FreeSenderBuffers(cplOcnPID, &context_id);
  }
  if (ocnComm != MPI_COMM_NULL)
  {
    char outputFileOcn[] = "OcnWithProj2.h5m";
    ierr = iMOAB_WriteMesh(cmpOcnPID, outputFileOcn, fileWriteOptions,
        strlen(outputFileOcn), strlen(fileWriteOptions) );
  }
#endif

  MPI_Barrier(MPI_COMM_WORLD);

#ifdef ENABLE_ATMLND_COUPLING
  // start land proj:

// we used this to compute
   //  ierr = iMOAB_ComputeCommGraph(cmpPhAtmPID, cplLndPID, &atmCouComm, &atmPEGroup, &couPEGroup,
    //     &typeA, &typeB, &cmpatm, &cpllnd);

  // end copy
  PUSH_TIMER("Send/receive data from phys comp atm to coupler land, using computed graph")
  if (atmComm != MPI_COMM_NULL ){

     // as always, use nonblocking sends
    // this is for projection to land:
     ierr = iMOAB_SendElementTag(cmpPhAtmPID, "T_ph;u_ph;v_ph;", &atmCouComm, &cpllnd, strlen("T_ph;u_ph;v_ph;"));
     CHECKIERR(ierr, "cannot send tag values towards cpl on land")
  }
  if (couComm != MPI_COMM_NULL){
    // receive on lnd on coupler pes
    ierr = iMOAB_ReceiveElementTag(cplLndPID, "T_proj;u_proj;v_proj;", &atmCouComm, &cmpatm, strlen("T_proj;u_proj;v_proj;"));
    CHECKIERR(ierr, "cannot receive tag values on land on coupler")
  }
  POP_TIMER(MPI_COMM_WORLD, rankInGlobalComm)

  // we can now free the sender buffers
  if (atmComm != MPI_COMM_NULL) {
    ierr = iMOAB_FreeSenderBuffers(cmpPhAtmPID, &cpllnd);
    CHECKIERR(ierr, "cannot free buffers used to resend atm tag towards the land on coupler")
  }


#ifdef VERBOSE
  if (couComm != MPI_COMM_NULL){
    char outputFileTgtLnd[] = "fLndOnCpl.h5m";
    ierr = iMOAB_WriteMesh(cplLndPID, outputFileTgtLnd, fileWriteOptions,
      strlen(outputFileTgtLnd), strlen(fileWriteOptions) );
  }
#endif

  // end land proj
  // send the tags back to land pes, from land mesh on coupler pes
  // send from cplLndPID to cmpLndPID, using common joint comm
      // as always, use nonblocking sends
  // original graph
  // int context_id = -1;
  // the land might not have these tags yet; it should be a different name for land
  // in e3sm we do have different names
  if (lndComm != MPI_COMM_NULL) {
    int tagIndexIn2;
    ierr = iMOAB_DefineTagStorage(cmpLndPID, bottomTempProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndexIn2,  strlen(bottomTempProjectedField) );
    CHECKIERR(ierr, "failed to define the field tag for receiving back the tag a2oTbot_proj on lnd pes");
    ierr = iMOAB_DefineTagStorage(cmpLndPID, bottomUVelProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndexIn2,  strlen(bottomUVelProjectedField) );
    CHECKIERR(ierr, "failed to define the field tag for receiving back the tag a2oUbot_proj on lnd pes");
    ierr = iMOAB_DefineTagStorage(cmpLndPID, bottomVVelProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndexIn2,  strlen(bottomVVelProjectedField) );
    CHECKIERR(ierr, "failed to define the field tag for receiving back the tag a2oVbot_proj on lnd pes");
  }
  if (couComm != MPI_COMM_NULL){
    ierr = iMOAB_SendElementTag(cplLndPID, "T_proj;u_proj;v_proj;", &lndCouComm, &context_id,
        strlen("T_proj;u_proj;v_proj;"));
    CHECKIERR(ierr, "cannot send tag values back to land pes")
  }
  // receive on component 3, land
  if (lndComm != MPI_COMM_NULL)
  {
    ierr = iMOAB_ReceiveElementTag(cmpLndPID, "T_proj;u_proj;v_proj;",
        &lndCouComm, &context_id, strlen("T_proj;u_proj;v_proj;"));
    CHECKIERR(ierr, "cannot receive tag values from land mesh on coupler pes")
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (couComm != MPI_COMM_NULL){
    ierr = iMOAB_FreeSenderBuffers(cplLndPID, &context_id);
  }
  if (lndComm != MPI_COMM_NULL)
  {
    char outputFileLnd[] = "LndWithProj2.h5m";
    ierr = iMOAB_WriteMesh(cmpLndPID, outputFileLnd, fileWriteOptions,
        strlen(outputFileLnd), strlen(fileWriteOptions) );
  }

#endif // ENABLE_ATMLND_COUPLING


#ifdef ENABLE_ATMOCN_COUPLING
  if (couComm != MPI_COMM_NULL){
    ierr = iMOAB_DeregisterApplication(cplAtmOcnPID);
    CHECKIERR(ierr, "cannot deregister app intx AO" )
  }
#endif

#ifdef ENABLE_ATMLND_COUPLING
  if (lndComm != MPI_COMM_NULL) {
      ierr = iMOAB_DeregisterApplication(cmpLndPID);
      CHECKIERR(ierr, "cannot deregister app LND1" )
    }
#endif // ENABLE_ATMLND_COUPLING
  if (atmComm != MPI_COMM_NULL) {
    ierr = iMOAB_DeregisterApplication(cmpPhAtmPID);
    CHECKIERR(ierr, "cannot deregister app PhysAtm " )
  }
#ifdef ENABLE_ATMOCN_COUPLING
  if (ocnComm != MPI_COMM_NULL) {
    ierr = iMOAB_DeregisterApplication(cmpOcnPID);
    CHECKIERR(ierr, "cannot deregister app OCN1" )
  }
#endif // ENABLE_ATMOCN_COUPLING

  if (atmComm != MPI_COMM_NULL) {
    ierr = iMOAB_DeregisterApplication(cmpAtmPID);
    CHECKIERR(ierr, "cannot deregister app ATM1" )
  }

#ifdef ENABLE_ATMLND_COUPLING
  if (couComm != MPI_COMM_NULL){
    ierr = iMOAB_DeregisterApplication(cplLndPID);
    CHECKIERR(ierr, "cannot deregister app LNDX" )
  }
#endif // ENABLE_ATMLND_COUPLING

#ifdef ENABLE_ATMOCN_COUPLING
  if (couComm != MPI_COMM_NULL){
    ierr = iMOAB_DeregisterApplication(cplOcnPID);
    CHECKIERR(ierr, "cannot deregister app OCNX" )
  }
#endif // ENABLE_ATMOCN_COUPLING

  if (couComm != MPI_COMM_NULL){
    ierr = iMOAB_DeregisterApplication(cplAtmPID);
    CHECKIERR(ierr, "cannot deregister app ATMX" )
  }

//#endif
  ierr = iMOAB_Finalize();
  CHECKIERR(ierr, "did not finalize iMOAB" )

  // free atm coupler group and comm
  if (MPI_COMM_NULL != atmCouComm) MPI_Comm_free(&atmCouComm);
  MPI_Group_free(&joinAtmCouGroup);
  if (MPI_COMM_NULL != atmComm) MPI_Comm_free(&atmComm);

#ifdef ENABLE_ATMOCN_COUPLING
  if (MPI_COMM_NULL != ocnComm) MPI_Comm_free(&ocnComm);
  // free ocn - coupler group and comm
  if (MPI_COMM_NULL != ocnCouComm) MPI_Comm_free(&ocnCouComm);
  MPI_Group_free(&joinOcnCouGroup);
#endif

#ifdef ENABLE_ATMLND_COUPLING
  if (MPI_COMM_NULL != lndComm) MPI_Comm_free(&lndComm);
  // free land - coupler group and comm
  if (MPI_COMM_NULL != lndCouComm) MPI_Comm_free(&lndCouComm);
  MPI_Group_free(&joinLndCouGroup);
#endif

  if (MPI_COMM_NULL != couComm)    MPI_Comm_free(&couComm);

  MPI_Group_free(&atmPEGroup);
#ifdef ENABLE_ATMOCN_COUPLING
  MPI_Group_free(&ocnPEGroup);
#endif
#ifdef ENABLE_ATMLND_COUPLING
  MPI_Group_free(&lndPEGroup);
#endif
  MPI_Group_free(&couPEGroup);
  MPI_Group_free(&jgroup);

  MPI_Finalize();
// endif #if 0

  return 0;
}







