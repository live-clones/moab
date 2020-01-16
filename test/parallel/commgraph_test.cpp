/*
 * This commgraph_test primary goal is to setup a communication framework between
 * two components, based only on a marker associated to the data being migrated
 * between the components
 * One use case for us is the physics grid in the atmosphere and spectral dynamics grid in the
 * atmosphere
 * These could be in theory on different sets of PES, although these 2 models are on the same PES
 *
 * In our case, spectral element data is hooked using the GLOBAL_DOFS tag of 4x4 ids, associated to
 * element
 * PhysGrid is coming from an AtmPhys.h5m point cloud grid partitioned in 16
 * Spectral mesh is our regular wholeATM_T_01.h5m, which is after one time step
 *
 * phys grid is very sprinkled in the partition, spectral mesh is more compact; For correct
 * identification/correspondence in parallel, it would make sense to use boxes for the spectral mesh
 *
 * We employ our friends the crystal router, in which we use rendezvous algorithm, to set
 * up the communication pattern
 *
 * input: wholeATM_T.h5m file, on 128 procs, the same file that is used by imoab_coupler test
 * input: AtmPhys.h5m file, which contains the physics grid, distributed on 16 processes
 * input 2: wholeLND.h5m, which is land distributed on 16 processes too
 *
 * The communication pattern will be established using a rendezvous method, based on the marker
 * (in this case, GLOBAL_ID on vertices on phys grid and GLOBAL_DOFS tag on spectral elements)
 *
 * in the end, we need to modify tag migrate to move data between these types of components, by
 * ID
 *
 * wholeLnd.h5m has holes in the ID space, that we need to account for;
 * In the end, this could be used to send data directly from Dyn atm to land; or to lnd on coupler ?
 *
 *
 */

#include "moab/Core.hpp"

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

#define CHECKIERR(rc, message)  if (0!=rc) { printf ("%s. ErrorCode = %d\n", message, rc); CHECK(0);}

using namespace moab;

// #define VERBOSE



//declare some variables outside main method
// easier to pass them around to the test
int ierr;
int rankInGlobalComm, numProcesses;
MPI_Group jgroup;
std::string atmFilename = TestDir + "/wholeATM_T.h5m";
// on a regular case,  5 ATM
// cmpatm is for atm on atm pes ! it has the spectral mesh
// cmpphys is for atm on atm phys pes ! it has the point cloud , phys grid

//
int rankInAtmComm = -1;
// it is the spectral mesh unique comp id
int cmpatm=605;  // component ids are unique over all pes, and established in advance;

std::string atmPhysFilename = TestDir + "/AtmPhys.h5m";
std::string atmPhysOutFilename = "outPhys.h5m";
std::string atmFilename2 = "wholeATM_new.h5m";
int rankInPhysComm = -1;
// this will be the physics atm com id; it should be actually 5
int physatm = 5;  // component ids are unique over all pes, and established in advance;

//int rankInJointComm = -1;

int nghlay=0; // number of ghost layers for loading the file

std::vector<int> groupTasks;
int startG1=0, startG2=0, endG1=numProcesses-1, endG2=numProcesses-1;
int typeA = 1; // spectral mesh, with GLOBAL_DOFS tags on cells
int typeB = 2; // point cloud mesh, with GLOBAL_ID tag on vertices

std::string readopts("PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS");
std::string readoptsPC("PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION");
std::string fileWriteOptions("PARALLEL=WRITE_PART");
std::string tagT("a2oTbot");
std::string tagU("a2oUbot");
std::string tagV("a2oVbot");
std::string separ(";");
std::string tagT1("a2oTbot_1");
std::string tagU1("a2oUbot_1");
std::string tagV1("a2oVbot_1");
std::string tagT2("a2oTbot_2"); // just one send



int commgraphtest();

void testspectral_phys()
{
  // no changes
  commgraphtest();
}

void testspectral_lnd()
{
  // first model is spectral, second is land
  atmPhysFilename= TestDir + "/wholeLnd.h5m";
  atmPhysOutFilename = std::string("outLnd.h5m");
  atmFilename2 = std::string("wholeATM_lnd.h5m");
  commgraphtest();
}

void testphysatm_lnd()
{
  //use for first file the output "outPhys.h5m" from first test
  atmFilename = std::string("outPhys.h5m");
  atmPhysFilename= std::string("outLnd.h5m");
  atmPhysOutFilename = std::string("physAtm_lnd.h5m");
  atmFilename2 = std::string("physBack_lnd.h5m");
  tagT = tagT1;
  tagU = tagU1;
  tagV = tagV1;
  tagT1 = std::string("newT");
  tagT2 = std::string("newT2");
  typeA = 2;
  commgraphtest();

}
int commgraphtest()
{

  if (!rankInGlobalComm)
  {
    std::cout << " first  file: " << atmFilename << "\n   on tasks : " << startG1 << ":"<<endG1 <<
        "\n second file: " << atmPhysFilename << "\n     on tasks : " << startG2 << ":" << endG2 << "\n  ";
  }

  // load files on 2 different communicators, groups
  // will create a joint comm for rendezvous
  MPI_Group atmPEGroup;
  groupTasks.resize(numProcesses, 0);
  for (int i=startG1; i<=endG1; i++)
    groupTasks [i-startG1] = i;

  ierr = MPI_Group_incl(jgroup, endG1-startG1+1, &groupTasks[0], &atmPEGroup);
  CHECKIERR(ierr, "Cannot create atmPEGroup")

  groupTasks.clear();
  groupTasks.resize(numProcesses, 0);
  MPI_Group atmPhysGroup;
  for (int i=startG2; i<=endG2; i++)
    groupTasks [i-startG2] = i;

  ierr = MPI_Group_incl(jgroup, endG2-startG2+1, &groupTasks[0], &atmPhysGroup);
  CHECKIERR(ierr, "Cannot create atmPhysGroup")

  // create 2 communicators, one for each group
  int ATM_COMM_TAG = 1;
  MPI_Comm atmComm;
  // atmComm is for atmosphere app;
  ierr = MPI_Comm_create_group(MPI_COMM_WORLD, atmPEGroup, ATM_COMM_TAG, &atmComm);
  CHECKIERR(ierr, "Cannot create atmComm")

  int PHYS_COMM_TAG = 2;
  MPI_Comm physComm;
  // physComm is for phys atm app
  ierr = MPI_Comm_create_group(MPI_COMM_WORLD, atmPhysGroup, PHYS_COMM_TAG, &physComm);
  CHECKIERR(ierr, "Cannot create physComm")

  // now, create the joint communicator atm physatm

  //
  MPI_Group joinAtmPhysAtmGroup;
  ierr = MPI_Group_union(atmPEGroup, atmPhysGroup, &joinAtmPhysAtmGroup);
  CHECKIERR(ierr, "Cannot create joint atm - phys atm group")
  int JOIN_COMM_TAG = 5;
  MPI_Comm joinComm;
  ierr = MPI_Comm_create_group(MPI_COMM_WORLD, joinAtmPhysAtmGroup, JOIN_COMM_TAG, &joinComm);
  CHECKIERR(ierr, "Cannot create joint atm cou communicator")

  ierr = iMOAB_Initialize(0, NULL); // not really needed anything from argc, argv, yet; maybe we should
  CHECKIERR(ierr, "Cannot initialize iMOAB")

  int cmpAtmAppID =-1;
  iMOAB_AppID cmpAtmPID=&cmpAtmAppID; // atm
  int physAtmAppID=-1;// -1 means it is not initialized
  iMOAB_AppID physAtmPID=&physAtmAppID; // phys atm on phys pes

  // load atm mesh
  if (atmComm != MPI_COMM_NULL)
  {
    MPI_Comm_rank( atmComm, &rankInAtmComm );
    ierr = iMOAB_RegisterApplication("ATM1", &atmComm, &cmpatm, cmpAtmPID);
    CHECKIERR(ierr, "Cannot register ATM App")

    // load first model
    std::string rdopts = readopts;
    if (typeA==2)
      rdopts = readoptsPC; // point cloud
    ierr = iMOAB_LoadMesh(cmpAtmPID, atmFilename.c_str(), rdopts.c_str(), &nghlay, atmFilename.length(), rdopts.length() );
    CHECKIERR(ierr, "Cannot load ATM mesh")
  }

  // load atm phys mesh
  if (physComm != MPI_COMM_NULL)
  {
    MPI_Comm_rank( physComm, &rankInPhysComm );
    ierr = iMOAB_RegisterApplication("PhysATM", &physComm, &physatm, physAtmPID);
    CHECKIERR(ierr, "Cannot register PHYS ATM App")

    // load phys atm mesh all tests  this is PC
    ierr = iMOAB_LoadMesh(physAtmPID, atmPhysFilename.c_str(), readoptsPC.c_str(), &nghlay, atmPhysFilename.length(), readoptsPC.length() );
    CHECKIERR(ierr, "Cannot load Phys ATM mesh")
  }


  if (MPI_COMM_NULL != joinComm)
  {
    ierr = iMOAB_ComputeCommGraph(cmpAtmPID, physAtmPID, &joinComm, &atmPEGroup, &atmPhysGroup,
          &typeA, &typeB, &cmpatm, &physatm);
    // it will generate parcomm graph between atm and atmPhys models
    // 2 meshes, that are distributed in parallel
    CHECKIERR(ierr, "Cannot compute comm graph between the two apps ")
  }

  if (atmComm != MPI_COMM_NULL)
  {
    // call send tag;
    std::string tags=tagT+separ+tagU+separ+tagV+separ;
    ierr = iMOAB_SendElementTag(cmpAtmPID, tags.c_str(), &joinComm, &physatm, tags.length());
    CHECKIERR(ierr, "cannot send tag values")
  }

  if (physComm != MPI_COMM_NULL)
  {
    // need to define tag storage
    std::string tags1=tagT1+separ+tagU1+separ+tagV1+separ;
    int tagType = DENSE_DOUBLE;
    int ndof=1;
    if (typeB==1)
      ndof=16;
    int tagIndex = 0;
    ierr = iMOAB_DefineTagStorage(physAtmPID, tagT1.c_str(), &tagType, &ndof, &tagIndex,  tagT1.length() );
    CHECKIERR(ierr, "failed to define the field tag a2oTbot");

    ierr = iMOAB_DefineTagStorage(physAtmPID, tagU1.c_str(), &tagType, &ndof, &tagIndex,  tagU1.length() );
    CHECKIERR(ierr, "failed to define the field tag a2oUbot");

    ierr = iMOAB_DefineTagStorage(physAtmPID, tagV1.c_str(), &tagType, &ndof, &tagIndex,  tagV1.length() );
    CHECKIERR(ierr, "failed to define the field tag a2oVbot");


    ierr = iMOAB_ReceiveElementTag(physAtmPID, tags1.c_str(), &joinComm, &cmpatm, tags1.length());
        CHECKIERR(ierr, "cannot receive tag values")
  }

  // we can now free the sender buffers
  if (atmComm != MPI_COMM_NULL) {
    ierr = iMOAB_FreeSenderBuffers(cmpAtmPID, &physatm);
    CHECKIERR(ierr, "cannot free buffers")
  }

  if (physComm != MPI_COMM_NULL)
  {
    ierr = iMOAB_WriteMesh(physAtmPID, (char*)atmPhysOutFilename.c_str(), (char*)fileWriteOptions.c_str(),
          atmPhysOutFilename.length(), fileWriteOptions.length() );
  }
  if (physComm != MPI_COMM_NULL)
  {
    // send back first tag only
    ierr = iMOAB_SendElementTag(physAtmPID, tagT1.c_str(), &joinComm, &cmpatm, tagT1.length());
    CHECKIERR(ierr, "cannot send tag values")
  }
  // receive it in a different tag
  if (atmComm != MPI_COMM_NULL) {
    // need to define tag storage
    int tagType = DENSE_DOUBLE;
    int ndof=16;
    if (typeA==2)
      ndof=1;
    int tagIndex = 0;
    ierr = iMOAB_DefineTagStorage(cmpAtmPID, tagT2.c_str(), &tagType, &ndof, &tagIndex,  tagT2.length() );
    CHECKIERR(ierr, "failed to define the field tag a2oTbot_2");

    ierr = iMOAB_ReceiveElementTag(cmpAtmPID, tagT2.c_str(), &joinComm, &physatm, tagT2.length());
    CHECKIERR(ierr, "cannot receive tag values a2oTbot_2")

  }
  // now send back one tag , into a different tag, and see if we get the same values back
  // we can now free the sender buffers
  if (physComm != MPI_COMM_NULL) {
    ierr = iMOAB_FreeSenderBuffers(physAtmPID, &cmpatm);
    CHECKIERR(ierr, "cannot free buffers ")
  }
  if (atmComm != MPI_COMM_NULL)
  {
    ierr = iMOAB_WriteMesh(cmpAtmPID, (char*)atmFilename2.c_str(), (char*)fileWriteOptions.c_str(),
        atmFilename2.length(), fileWriteOptions.length() );
  }

  // unregister in reverse order
  if (physComm != MPI_COMM_NULL){
    ierr = iMOAB_DeregisterApplication(physAtmPID);
    CHECKIERR(ierr, "cannot deregister second app model" )
  }

  if (atmComm != MPI_COMM_NULL) {
    ierr = iMOAB_DeregisterApplication(cmpAtmPID);
    CHECKIERR(ierr, "cannot deregister first app model" )
  }

  ierr = iMOAB_Finalize();
  CHECKIERR(ierr, "did not finalize iMOAB" )

  // free atm group and comm
  if (MPI_COMM_NULL != atmComm) MPI_Comm_free(&atmComm);
  MPI_Group_free(&atmPEGroup);

  // free atm phys group and comm
  if (MPI_COMM_NULL != physComm) MPI_Comm_free(&physComm);
  MPI_Group_free(&atmPhysGroup);

  // free atm phys group and comm
  if (MPI_COMM_NULL != joinComm) MPI_Comm_free(&joinComm);
  MPI_Group_free(&joinAtmPhysAtmGroup);

  return 0;
}
int main( int argc, char* argv[] )
{


  MPI_Init(&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &rankInGlobalComm );
  MPI_Comm_size( MPI_COMM_WORLD, &numProcesses );

  MPI_Comm_group(MPI_COMM_WORLD, &jgroup);// all processes in global group

  // default: load atm on 2 proc, phys grid on 2 procs, establish comm graph, then migrate
  // data from atm pes to phys pes, and back
  startG1=0, startG2=0, endG1=numProcesses-1, endG2=numProcesses-1;

  ProgOptions opts;
  opts.addOpt<std::string>("modelA,m", "first model file ", &atmFilename);
  opts.addOpt<int>("typeA,t", " type of first model ", &typeA);


  opts.addOpt<std::string>("modelB,n", "second model file", &atmPhysFilename);
  opts.addOpt<int>("typeB,v", " type of the second model ", &typeB);

  opts.addOpt<int>("startAtm,a", "start task for first model layout", &startG1);
  opts.addOpt<int>("endAtm,b", "end task for first model layout", &endG1);

  opts.addOpt<int>("startPhys,c", "start task for second model layout", &startG2);
  opts.addOpt<int>("endPhys,d", "end task for second model layout", &endG2);

  opts.addOpt<std::string>("output,o", "output filename", &atmPhysOutFilename);

  opts.addOpt<std::string>("output,o", "output filename", &atmFilename2);

  opts.parseCommandLine(argc, argv);

  int num_err = 0;
  num_err += RUN_TEST( testspectral_phys );

  //
  if (argc == 1)
  {
    num_err += RUN_TEST( testspectral_lnd );
    num_err += RUN_TEST( testphysatm_lnd );
  }
  MPI_Group_free(&jgroup);

  MPI_Finalize();

  return num_err;
}

