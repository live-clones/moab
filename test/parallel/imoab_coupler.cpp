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

// #define VERBOSE

int main( int argc, char* argv[] )
{
  int ierr;
  int rankInAtmComm, rankInOcnComm, rankInLndComm, rankInGlobalComm, numProcesses;
  MPI_Comm globalComm; // will be a copy of the global
  MPI_Group jgroup;
  std::string readopts("PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS");
  std::string readoptsLnd("PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION");

  // Timer data
  moab::CpuTimer timer;
  double timer_ops;
  std::string opName;

  int repartitioner_scheme = 0;
#ifdef MOAB_HAVE_ZOLTAN
  repartitioner_scheme = 2; // use the graph partitioner in that caseS
#endif

  MPI_Init(&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &rankInGlobalComm );
  MPI_Comm_size( MPI_COMM_WORLD, &numProcesses );

  MPI_Comm_dup(MPI_COMM_WORLD, &globalComm);
  MPI_Comm_group(globalComm, &jgroup);

  // on a regular case,  5 ATM, 6 CPLATM (ATMX), 17 OCN     , 18 CPLOCN (OCNX)  ; intx atm/ocn is not in e3sm yet, give a number
  //   6 * 100+ 18 = 618
  // 9 LND, 10 CPLLND
  //   6 * 100 + 10 = 610  atmlndid:
  // cmpatm is for atm on atm pes
  // cmpocn is for ocean, on ocean pe
  // cplatm is for atm on coupler pes
  // cplocn is for ocean on coupelr pes
  // atmocnid is for intx atm / ocn on coupler pes
  //
  int cmpatm=5, cmpocn=17, cplatm=6, cplocn=18, atmocnid=618, cpllnd=10, cmplnd=9, atmlndid=610;  // component ids are unique over all pes, and established in advance;
  
  int nghlay=0; // number of ghost layers for loading the file
  std::vector<int> groupTasks; // at most 2 tasks
  int startG1=0, startG2=0, endG1=numProcesses-1, endG2=numProcesses-1, startG3=startG1, endG3=endG1; // Support launch of imoab_coupler test on any combo of 2*x processes

  // int startG1=0, startG2=0, endG1=0, endG2=0; // Support launch of imoab_coupler test on any combo of 2*x processes
  // int startG1=0, startG2=0, endG1=1, endG2=0; // Support launch of imoab_coupler test on any combo of 2*x processes
  // int startG1=0, startG2=0, endG1=numProcesses-1, endG2=0; // Support launch of imoab_coupler test on any combo of 2*x processes
  // int startG1=0, startG2=0, endG1=0, endG2=1; // Support launch of imoab_coupler test on any combo of 2*x processes
  // int startG1=0, startG2=0, endG1=1, endG2=1; // Support launch of imoab_coupler test on any combo of 2*x processes
  // int startG1=0, startG2=0, endG1=numProcesses-1, endG2=numProcesses-1; // Support launch of imoab_coupler test on any combo of 2*x processes

  // load atm on 2 proc, ocean on 2, migrate both to 2 procs, then compute intx
  // later, we need to compute weight matrix with tempestremap

  std::string atmFilename, ocnFilename, lndFilename;

  atmFilename = TestDir + "/wholeATM_T.h5m";
  ocnFilename = TestDir + "/recMeshOcn.h5m";
  lndFilename = TestDir + "/wholeLnd.h5m";

  ProgOptions opts;
  opts.addOpt<std::string>("atmosphere,t", "atm mesh filename (source)", &atmFilename);
  opts.addOpt<std::string>("ocean,m", "ocean mesh filename (target)", &ocnFilename);
  opts.addOpt<std::string>("land,l", "land mesh filename (target)", &lndFilename);

  opts.addOpt<int>("startAtm,a", "start task for atmosphere layout", &startG1);
  opts.addOpt<int>("endAtm,b", "end task for atmosphere layout", &endG1);
  opts.addOpt<int>("startOcn,c", "start task for ocean layout", &startG2);
  opts.addOpt<int>("endOcn,d", "end task for ocean layout", &endG2);
  opts.addOpt<int>("startLnd,e", "start task for land layout", &startG3);
  opts.addOpt<int>("endLnd,f", "end task for land layout", &endG3);
  opts.addOpt<int>("partitioning,p", "partitioning option for migration", &repartitioner_scheme);

  opts.parseCommandLine(argc, argv);

  if (!rankInGlobalComm)
  {
    std::cout << " atm file: " << atmFilename << "\n   on tasks : " << startG1 << ":"<<endG1 <<
        "\n ocn file: " << ocnFilename << "\n     on tasks : " << startG2 << ":"<<endG2 <<
        "\n lnd file: " << lndFilename << "\n     on tasks : " << startG3 << ":"<<endG3 <<
        "\n  partitioning (0 trivial, 1 graph, 2 geometry) " << repartitioner_scheme << "\n  ";
  }

  // load files on 3 different communicators, groups
  // first groups has task 0, second group tasks 0 and 1
  // coupler will be on joint tasks, will be on a third group (0 and 1, again)
  MPI_Group atmPEGroup, ocnPEGroup, lndPEGroup;
  groupTasks.resize(numProcesses, 0);
  for (int i=startG1; i<=endG1; i++)
    groupTasks [i-startG1] = i;

  ierr = MPI_Group_incl(jgroup, endG1-startG1+1, &groupTasks[0], &atmPEGroup);
  CHECKIERR(ierr, "Cannot create atmPEGroup")

  groupTasks.clear();
  groupTasks.resize(numProcesses, 0);
  for (int i=startG2; i<=endG2; i++)
    groupTasks [i-startG2] = i;

  ierr = MPI_Group_incl(jgroup, endG2-startG2+1, &groupTasks[0], &ocnPEGroup);
  CHECKIERR(ierr, "Cannot create ocnPEGroup")

  groupTasks.clear();
  groupTasks.resize(numProcesses, 0);
  for (int i=startG3; i<=endG3; i++)
    groupTasks [i-startG3] = i;

  ierr = MPI_Group_incl(jgroup, endG3-startG3+1, &groupTasks[0], &lndPEGroup);
  CHECKIERR(ierr, "Cannot create lndPEGroup")

  // create 3 communicators, one for each group
  int ATM_COMM_TAG = 1, OCN_COMM_TAG = 2, LND_COMM_TAG = 3;
  MPI_Comm atmComm, ocnComm, lndComm;
  // atmComm is for atmosphere app;
  ierr = MPI_Comm_create_group(globalComm, atmPEGroup, ATM_COMM_TAG, &atmComm);
  CHECKIERR(ierr, "Cannot create atmComm")

  // ocnComm is for ocean app
  ierr = MPI_Comm_create_group(globalComm, ocnPEGroup, OCN_COMM_TAG, &ocnComm);
  CHECKIERR(ierr, "Cannot create ocnComm")

  // lndComm is for land app
  ierr = MPI_Comm_create_group(globalComm, lndPEGroup, LND_COMM_TAG, &lndComm);
  CHECKIERR(ierr, "Cannot create lndComm")

  ierr = iMOAB_Initialize(argc, argv); // not really needed anything from argc, argv, yet; maybe we should
  CHECKIERR(ierr, "Cannot initialize iMOAB")

  int cmpAtmAppID =-1;
  iMOAB_AppID cmpAtmPID=&cmpAtmAppID; // atm
  int cmpOcnAppID =-1;
  iMOAB_AppID cmpOcnPID=&cmpOcnAppID; // ocn
  int cplAtmAppID=-1, cplOcnAppID=-1, cplAtmOcnAppID=-1;// -1 means it is not initialized
  iMOAB_AppID cplAtmPID=&cplAtmAppID; // atm on coupler PEs
  iMOAB_AppID cplOcnPID=&cplOcnAppID; // ocn on coupler PEs
  iMOAB_AppID cplAtmOcnPID=&cplAtmOcnAppID; // intx atm -ocn on coupler PEs

  int cmpLndAppID =-1;
  iMOAB_AppID cmpLndPID=&cmpLndAppID; // lnd
  int cplLndAppID=-1, cplAtmLndAppID=-1;// -1 means it is not initialized
  iMOAB_AppID cplLndPID=&cplLndAppID; // land on coupler PEs
  iMOAB_AppID cplAtmLndPID=&cplAtmLndAppID; // intx atm - lnd on coupler PEs

  // Register all the applications on the coupler PEs
  ierr = iMOAB_RegisterApplication("ATMX", &globalComm, &cplatm, cplAtmPID); // atm on coupler pes
  CHECKIERR(ierr, "Cannot register ATM over coupler PEs")

  ierr = iMOAB_RegisterApplication("OCNX", &globalComm, &cplocn, cplOcnPID); // ocn on coupler pes
  CHECKIERR(ierr, "Cannot register OCN over coupler PEs")

  ierr = iMOAB_RegisterApplication("LNDX", &globalComm, &cpllnd, cplLndPID); // lnd on coupler pes
  CHECKIERR(ierr, "Cannot register LND over coupler PEs")
  
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
    ierr = iMOAB_SendMesh(cmpAtmPID, &globalComm, &jgroup, &cplatm, &repartitioner_scheme); // send to component 3, on coupler pes
    CHECKIERR(ierr, "cannot send elements" )
    POP_TIMER(atmComm, rankInAtmComm)
  }
  // now, receive meshes, on joint communicator; first mesh 1
  PUSH_TIMER("Receive ATM mesh and resolve shared entities")
  ierr = iMOAB_ReceiveMesh(cplAtmPID, &globalComm, &atmPEGroup, &cmpatm); // receive from component 1
  CHECKIERR(ierr, "cannot receive elements on ATMX app")
  POP_TIMER(globalComm, rankInGlobalComm)

  // we can now free the sender buffers
  if (atmComm != MPI_COMM_NULL) {
    ierr = iMOAB_FreeSenderBuffers(cmpAtmPID, &globalComm, &cplatm);
    CHECKIERR(ierr, "cannot free buffers used to send atm mesh")
  }

  MPI_Barrier(MPI_COMM_WORLD);

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
    ierr = iMOAB_SendMesh(cmpOcnPID, &globalComm, &jgroup, &cplocn, &repartitioner_scheme); // send to component 3, on coupler pes
    CHECKIERR(ierr, "cannot send elements" )
    POP_TIMER(ocnComm, rankInOcnComm)
  }

  PUSH_TIMER("Receive OCN mesh and resolve shared entities")
  ierr = iMOAB_ReceiveMesh(cplOcnPID, &globalComm, &ocnPEGroup, &cmpocn); // receive from component 2, ocn, on coupler pes
  CHECKIERR(ierr, "cannot receive elements on OCNX app")
  POP_TIMER(globalComm, rankInGlobalComm)

  if (ocnComm != MPI_COMM_NULL) {
    ierr = iMOAB_FreeSenderBuffers(cmpOcnPID, &globalComm, &cplocn);
    CHECKIERR(ierr, "cannot free buffers used to send ocn mesh")
  }

  MPI_Barrier(MPI_COMM_WORLD);

  char outputFileTgt3[] = "recvOcn.h5m";
  #ifdef MOAB_HAVE_MPI
    char writeOptions3[] ="PARALLEL=WRITE_PART";
  #else
    char writeOptions3[] ="";
  #endif
  PUSH_TIMER("Write migrated OCN mesh on coupler PEs")
  ierr = iMOAB_WriteMesh(cplOcnPID, outputFileTgt3, writeOptions3,
    strlen(outputFileTgt3), strlen(writeOptions3) );
  CHECKIERR(ierr, "cannot write ocn mesh after receiving")
  POP_TIMER(globalComm, rankInGlobalComm)

  if (lndComm != MPI_COMM_NULL) {
    MPI_Comm_rank( lndComm, &rankInLndComm );
    ierr = iMOAB_RegisterApplication("LND1", &lndComm, &cmplnd, cmpLndPID);
    CHECKIERR(ierr, "Cannot register LND App ")

    // load the next component mesh
    PUSH_TIMER("Load LND mesh")
    ierr = iMOAB_LoadMesh(cmpLndPID, lndFilename.c_str(), readoptsLnd.c_str(), &nghlay, lndFilename.length(), readoptsLnd.length() );
    CHECKIERR(ierr, "Cannot load LND mesh on cplLndPID")
    POP_TIMER(lndComm, rankInLndComm)

    int nverts[3], nelem[3];
    ierr = iMOAB_GetMeshInfo( cmpLndPID, nverts, nelem, 0, 0, 0);
    CHECKIERR(ierr, "failed to get mesh info");
    printf("Land Component Mesh: %d vertices and %d elements\n", nverts[0], nelem[0]);

    // then send mesh to coupler pes, on cplLndPID
    PUSH_TIMER("LND mesh: compute partition and send mesh")
    ierr = iMOAB_SendMesh(cmpLndPID, &globalComm, &jgroup, &cpllnd, &repartitioner_scheme); // send to component 3, on coupler pes
    CHECKIERR(ierr, "cannot send elements" )
    POP_TIMER(lndComm, rankInLndComm)
  }

  PUSH_TIMER("Receive LND mesh and resolve shared entities")
  ierr = iMOAB_ReceiveMesh(cplLndPID, &globalComm, &lndPEGroup, &cmplnd); // receive from component 2, lnd, on coupler pes
  CHECKIERR(ierr, "cannot receive elements on LNDX app")
  POP_TIMER(globalComm, rankInGlobalComm)

  if (ocnComm != MPI_COMM_NULL) {
    ierr = iMOAB_FreeSenderBuffers(cmpLndPID, &globalComm, &cpllnd);
    CHECKIERR(ierr, "cannot free buffers used to send lnd mesh")
  }

  MPI_Barrier(MPI_COMM_WORLD);

  char outputFileLnd[] = "recvLnd.h5m";
  PUSH_TIMER("Write migrated LND mesh on coupler PEs")
  ierr = iMOAB_WriteMesh(cplLndPID, outputFileLnd, writeOptions3,
    strlen(outputFileLnd), strlen(writeOptions3) );
  CHECKIERR(ierr, "cannot write lnd mesh after receiving")
  POP_TIMER(globalComm, rankInGlobalComm)

#ifdef MOAB_HAVE_TEMPESTREMAP
  // now compute intersection between OCNx and ATMx on coupler PEs
  ierr = iMOAB_RegisterApplication("ATMOCN", &globalComm, &atmocnid, cplAtmOcnPID);
  CHECKIERR(ierr, "Cannot register ocn_atm intx over coupler pes ")

  // now compute intersection between OCNx and ATMx on coupler PEs
  ierr = iMOAB_RegisterApplication("ATMLND", &globalComm, &atmlndid, cplAtmLndPID);
  CHECKIERR(ierr, "Cannot register ocn_atm intx over coupler pes ")

  const char* weights_identifiers[2] = {"scalar", "scalar-pc"};
  int disc_orders[3] = {4, 1, 1};
  const char* disc_methods[3] = {"cgll", "fv", "pcloud"};
  const char* dof_tag_names[3] = {"GLOBAL_DOFS", "GLOBAL_ID", "GLOBAL_ID"};

  PUSH_TIMER("Compute ATM-OCN mesh intersection")
  ierr = iMOAB_ComputeMeshIntersectionOnSphere(cplAtmPID, cplOcnPID, cplAtmOcnPID); // coverage mesh was computed here, for cplAtmPID, atm on coupler pes
  // basically, atm was redistributed according to target (ocean) partition, to "cover" the ocean partitions
  // check if intx valid, write some h5m intx file
  CHECKIERR(ierr, "cannot compute intersection" )
  POP_TIMER(globalComm, rankInGlobalComm)

  // the new graph will be for sending data from atm comp to coverage mesh;
  // it involves initial atm app; cmpAtmPID; also migrate atm mesh on coupler pes, cplAtmPID
  // results are in cplAtmOcnPID, intx mesh; remapper also has some info about coverage mesh
  // after this, the sending of tags from atm pes to coupler pes will use the new par comm graph, that has more precise info about
  // what to send for ocean cover ; every time, we will
  //  use the element global id, which should uniquely identify the element
  PUSH_TIMER("Compute OCN coverage graph for ATM mesh")
  ierr = iMOAB_CoverageGraph(&globalComm, cmpAtmPID,  &cmpatm, cplAtmPID,  &cplatm, cplAtmOcnPID, &cplocn); // it happens over joint communicator
  CHECKIERR(ierr, "cannot recompute direct coverage graph for ocean" )
  POP_TIMER(globalComm, rankInGlobalComm)

  PUSH_TIMER("Compute ATM-LND mesh intersection")
  ierr = iMOAB_ComputePointDoFIntersection(cplAtmPID, cplLndPID, cplAtmLndPID);
  CHECKIERR(ierr, "failed to compute point-cloud mapping");
  POP_TIMER(globalComm, rankInGlobalComm)

  // the new graph will be for sending data from atm comp to coverage mesh for land mesh;
  // it involves initial atm app; cmpAtmPID; also migrate atm mesh on coupler pes, cplAtmPID
  // results are in cplAtmLndPID, intx mesh; remapper also has some info about coverage mesh
  // after this, the sending of tags from atm pes to coupler pes will use the new par comm graph, that has more precise info about
  // what to send (specifically for land cover); every time,
  /// we will use the element global id, which should uniquely identify the element
  PUSH_TIMER("Compute LND coverage graph for ATM mesh")
  ierr = iMOAB_CoverageGraph(&globalComm, cmpAtmPID,  &cmpatm, cplAtmPID,  &cplatm, cplAtmLndPID, &cpllnd); // it happens over joint communicator
  CHECKIERR(ierr, "cannot recompute direct coverage graph for land" )
  POP_TIMER(globalComm, rankInGlobalComm)

  MPI_Barrier(MPI_COMM_WORLD);

#ifdef VERBOSE
  std::stringstream outf;
  outf<<"intx_0" << rankInGlobalComm<<".h5m";
  std::string intxfile=outf.str(); // write in serial the intx file, for debugging
  char writeOptions[] ="";
  ierr = iMOAB_WriteMesh(cplAtmOcnPID, (char*)intxfile.c_str(), writeOptions, (int)intxfile.length(), strlen(writeOptions));
  CHECKIERR(ierr, "cannot write intx file result" )
#endif


  int fMonotoneTypeID=0, fVolumetric=0, fValidate=1, fNoConserve=0;
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
  POP_TIMER(globalComm, rankInGlobalComm)

  MPI_Barrier(MPI_COMM_WORLD);

  /* Compute the weights to preoject the solution from ATM component to LND compoenent */
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
  POP_TIMER(globalComm, rankInGlobalComm)

  const char* bottomTempField = "a2oTbot";
  const char* bottomTempProjectedField = "a2oTbot_proj";
  int tagIndex[2];
  int tagTypes[2] = { DENSE_DOUBLE, DENSE_DOUBLE } ;
  int atmCompNDoFs = disc_orders[0]*disc_orders[0], ocnCompNDoFs = 1/*FV*/;

  ierr = iMOAB_DefineTagStorage(cplAtmPID, bottomTempField, &tagTypes[0], &atmCompNDoFs, &tagIndex[0],  strlen(bottomTempField) );
  CHECKIERR(ierr, "failed to define the field tag a2oTbot");
  ierr = iMOAB_DefineTagStorage(cplOcnPID, bottomTempProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndex[1],  strlen(bottomTempProjectedField) );
  CHECKIERR(ierr, "failed to define the field tag a2oTbot_proj");

  // two more fields for testing : U, V
  const char* bottomUVelField = "a2oUbot";
  const char* bottomUVelProjectedField = "a2oUbot_proj";
  const char* bottomVVelField = "a2oVbot";
  const char* bottomVVelProjectedField = "a2oVbot_proj";

  // Define more fields
  ierr = iMOAB_DefineTagStorage(cplAtmPID, bottomUVelField, &tagTypes[0], &atmCompNDoFs, &tagIndex[0],  strlen(bottomUVelField) );
  CHECKIERR(ierr, "failed to define the field tag a2oUbot");
  ierr = iMOAB_DefineTagStorage(cplOcnPID, bottomUVelProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndex[1],  strlen(bottomUVelProjectedField) );
  CHECKIERR(ierr, "failed to define the field tag a2oUbot_proj");
  ierr = iMOAB_DefineTagStorage(cplAtmPID, bottomVVelField, &tagTypes[0], &atmCompNDoFs, &tagIndex[0],  strlen(bottomVVelField) );
  CHECKIERR(ierr, "failed to define the field tag a2oUbot");
  ierr = iMOAB_DefineTagStorage(cplOcnPID, bottomVVelProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndex[1],  strlen(bottomVVelProjectedField) );
  CHECKIERR(ierr, "failed to define the field tag a2oUbot_proj");

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

  PUSH_TIMER("Send/receive data from component to coupler")
  if (atmComm != MPI_COMM_NULL ){
     // as always, use nonblocking sends
    // this is for projection to ocean:
     ierr = iMOAB_SendElementTag(cmpAtmPID, &cmpatm, &cplatm, "a2oTbot;a2oUbot;a2oVbot;", &globalComm, &cplocn, strlen("a2oTbot;a2oUbot;a2oVbot;"));
     CHECKIERR(ierr, "cannot send tag values")
  }
  // receive on atm on coupler pes, that was redistributed according to coverage
  ierr = iMOAB_ReceiveElementTag(cplAtmPID, &cmpatm, &cplatm, "a2oTbot;a2oUbot;a2oVbot;", &globalComm, &cplocn, strlen("a2oTbot;a2oUbot;a2oVbot;"));
  CHECKIERR(ierr, "cannot receive tag values")
  POP_TIMER(globalComm, rankInGlobalComm)

  // we can now free the sender buffers
  if (atmComm != MPI_COMM_NULL) {
    ierr = iMOAB_FreeSenderBuffers(cmpAtmPID, &globalComm, &cplatm);
    CHECKIERR(ierr, "cannot free buffers used to resend atm mesh tag towards the coverage mesh")
  }
#ifdef VERBOSE
    char outputFileRecvd[] = "recvAtmCoupOcn.h5m";
    ierr = iMOAB_WriteMesh(cplAtmPID, outputFileRecvd, writeOptions3,
        strlen(outputFileRecvd), strlen(writeOptions3) );
#endif

  /* We have the remapping weights now. Let us apply the weights onto the tag we defined
     on the source mesh and get the projection on the target mesh */
  PUSH_TIMER("Apply Scalar projection weights")
  const char * concat_fieldname = "a2oTbot;a2oUbot;a2oVbot;";
  const char * concat_fieldnameT = "a2oTbot_proj;a2oUbot_proj;a2oVbot_proj;";
  ierr = iMOAB_ApplyScalarProjectionWeights ( cplAtmOcnPID, weights_identifiers[0],
                                            concat_fieldname,
                                            concat_fieldnameT,
                                            strlen(weights_identifiers[0]),
                                            strlen(concat_fieldname),
                                            strlen(concat_fieldnameT)
                                            );
  CHECKIERR(ierr, "failed to compute projection weight application");
  POP_TIMER(globalComm, rankInGlobalComm)

  char outputFileTgt[] = "fOcnOnCpl.h5m";
  char writeOptions2[] ="PARALLEL=WRITE_PART";

  ierr = iMOAB_WriteMesh(cplOcnPID, outputFileTgt, writeOptions2,
    strlen(outputFileTgt), strlen(writeOptions2) );

  // send the projected tag back to ocean pes, with send/receive tag
  if (ocnComm != MPI_COMM_NULL) {
    int tagIndexIn2;
    ierr = iMOAB_DefineTagStorage(cmpOcnPID, bottomTempProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndexIn2,  strlen(bottomTempProjectedField) );
    CHECKIERR(ierr, "failed to define the field tag for receiving back the tag a2oTbot_proj on ocn pes");
    ierr = iMOAB_DefineTagStorage(cmpOcnPID, bottomUVelProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndexIn2,  strlen(bottomUVelProjectedField) );
    CHECKIERR(ierr, "failed to define the field tag for receiving back the tag a2oUbot_proj on ocn pes");
    ierr = iMOAB_DefineTagStorage(cmpOcnPID, bottomVVelProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndexIn2,  strlen(bottomVVelProjectedField) );
    CHECKIERR(ierr, "failed to define the field tag for receiving back the tag a2oVbot_proj on ocn pes");

  }
  // send the tag to ocean pes, from ocean mesh on coupler pes
  // first, send from cplocn to cplocn, from ocnComm, using common joint comm
      // as always, use nonblocking sends
  // original graph
  int other = -1;
  ierr = iMOAB_SendElementTag(cplOcnPID, &cplocn, &cmpocn, "a2oTbot_proj;a2oUbot_proj;a2oVbot_proj;", &globalComm, &other,
      strlen("a2oTbot_proj;a2oUbot_proj;a2oVbot_proj;"));
  CHECKIERR(ierr, "cannot send tag values back to ocean pes")

  // receive on component 2, ocean
  if (ocnComm != MPI_COMM_NULL)
  {
    ierr = iMOAB_ReceiveElementTag(cmpOcnPID, &cplocn, &cmpocn, "a2oTbot_proj;a2oUbot_proj;a2oVbot_proj;",
        &globalComm, &other, strlen("a2oTbot_proj;a2oUbot_proj;a2oVbot_proj;"));
    CHECKIERR(ierr, "cannot receive tag values from ocean mesh on coupler pes")
  }

  MPI_Barrier(MPI_COMM_WORLD);

  ierr = iMOAB_FreeSenderBuffers(cplOcnPID, &globalComm, &cmpocn);
  if (ocnComm != MPI_COMM_NULL)
  {
    char outputFileOcn[] = "OcnWithProj.h5m";
    ierr = iMOAB_WriteMesh(cmpOcnPID, outputFileOcn, writeOptions2,
        strlen(outputFileOcn), strlen(writeOptions2) );
  }

// start land proj:

  PUSH_TIMER("Send/receive data from component atm to coupler, in land context")
  if (atmComm != MPI_COMM_NULL ){

     // as always, use nonblocking sends
    // this is for projection to land:
     ierr = iMOAB_SendElementTag(cmpAtmPID, &cmpatm, &cplatm, "a2oTbot;a2oUbot;a2oVbot;", &globalComm, &cpllnd, strlen("a2oTbot;a2oUbot;a2oVbot;"));
     CHECKIERR(ierr, "cannot send tag values")
  }
  // receive on atm on coupler pes, that was redistributed according to coverage, for land context
  ierr = iMOAB_ReceiveElementTag(cplAtmPID, &cmpatm, &cplatm, "a2oTbot;a2oUbot;a2oVbot;", &globalComm, &cpllnd, strlen("a2oTbot;a2oUbot;a2oVbot;"));
  CHECKIERR(ierr, "cannot receive tag values")
  POP_TIMER(globalComm, rankInGlobalComm)

  // we can now free the sender buffers
  if (atmComm != MPI_COMM_NULL) {
    ierr = iMOAB_FreeSenderBuffers(cmpAtmPID, &globalComm, &cplatm);
    CHECKIERR(ierr, "cannot free buffers used to resend atm mesh tag towards the coverage mesh")
  }
#ifdef VERBOSE
    char outputFileRecvd[] = "recvAtmCoupLnd.h5m";
    ierr = iMOAB_WriteMesh(cplAtmPID, outputFileRecvd, writeOptions3,
        strlen(outputFileRecvd), strlen(writeOptions3) );
#endif

  /* We have the remapping weights now. Let us apply the weights onto the tag we defined
     on the source mesh and get the projection on the target mesh */
  PUSH_TIMER("Apply Scalar projection weights for land")
  ierr = iMOAB_ApplyScalarProjectionWeights ( cplAtmLndPID, weights_identifiers[1],
                                            concat_fieldname,
                                            concat_fieldnameT,
                                            strlen(weights_identifiers[1]),
                                            strlen(concat_fieldname),
                                            strlen(concat_fieldnameT)
                                            );
  CHECKIERR(ierr, "failed to compute projection weight application");
  POP_TIMER(globalComm, rankInGlobalComm)

  char outputFileTgt2[] = "fLndOnCpl.h5m";

  ierr = iMOAB_WriteMesh(cplLndPID, outputFileTgt2, writeOptions2,
    strlen(outputFileTgt2), strlen(writeOptions2) );

  // end land proj

  ierr = iMOAB_DeregisterApplication(cplAtmOcnPID);
  CHECKIERR(ierr, "cannot deregister app intx AO" )

#endif

  if (ocnComm != MPI_COMM_NULL) {
    ierr = iMOAB_DeregisterApplication(cmpOcnPID);
    CHECKIERR(ierr, "cannot deregister app OCN1" )
  }
  if (atmComm != MPI_COMM_NULL) {
    ierr = iMOAB_DeregisterApplication(cmpAtmPID);
    CHECKIERR(ierr, "cannot deregister app ATM1" )
  }

  ierr = iMOAB_DeregisterApplication(cplOcnPID);
  CHECKIERR(ierr, "cannot deregister app OCNX" )

  ierr = iMOAB_DeregisterApplication(cplAtmPID);
  CHECKIERR(ierr, "cannot deregister app OCNX" )

  ierr = iMOAB_Finalize();
  CHECKIERR(ierr, "did not finalize iMOAB" )

  if (MPI_COMM_NULL != atmComm) MPI_Comm_free(&atmComm);
  if (MPI_COMM_NULL != ocnComm) MPI_Comm_free(&ocnComm);

  MPI_Group_free(&atmPEGroup);
  MPI_Group_free(&ocnPEGroup);
  MPI_Group_free(&jgroup);
  MPI_Comm_free(&globalComm);

  return 0;
}


