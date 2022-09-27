/*
 * This imoab_setwithGid test will test the method
 * iMOAB_SetDoubleTagStorageWithGid, that is used in parallel
 * basically, 2 moab apps have different distribution of the same mesh, and use the method to
 * set the tag from one app to the other
 * 2 meshes are distributed differently, and each will be loaded in an app. GlobalId tag are the same
 * for both
 * example came from mct instance, and need to set the frac tag according to our moab variant of the
 * coupler distribution of the same mesh.

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
// CHECKIERR
#include "imoab_coupler_utils.hpp"

using namespace moab;

int main( int argc, char* argv[] )
{
    int ierr;
    int rankInGlobalComm, numProcesses;
    MPI_Group jgroup;
    std::string readoptsLnd( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION" );
    std::string readopts( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS" );

    // Timer data
    moab::CpuTimer timer;
    double timer_ops;
    std::string opName;

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rankInGlobalComm );
    MPI_Comm_size( MPI_COMM_WORLD, &numProcesses );

    MPI_Comm dup_comm_world;
    MPI_Comm_dup( MPI_COMM_WORLD, &dup_comm_world );

    std::string lndMoab = TestDir + "unittest/recLand.h5m";
    std::string lndMct = TestDir + "unittest/WHOLE_cx_lnd.h5m";
    std::string tagname("frac");
    int entType = 0;
    int ntags = 1;
    int nghlay = 0; // no ghosts

    ierr = iMOAB_Initialize( argc, argv );  // not really needed anything from argc, argv, yet; maybe we should
    CHECKIERR( ierr, "Cannot initialize iMOAB" )

    ProgOptions opts;
    opts.addOpt< std::string >( "source,s", "mct mesh filename (source)", &lndMct );

    opts.addOpt< std::string >( "target,t", "moab mesh filename (target)", &lndMoab );

    opts.addOpt< std::string >( "tags,g", "list of tags, colon separated", &tagname );

    opts.addOpt< int >( "enttype,e", "entity type for moab mesh (0 point cloud, 1 cells)", &entType );

    opts.addOpt< int >( "ntags,n", "number of tags to set", &ntags );

    opts.parseCommandLine( argc, argv );

    int cplLndAppID = -1, cplLnd2AppID = -1;   // -1 means it is not initialized
    iMOAB_AppID cplLndPID    = &cplLndAppID;     // land on coupler PEs, moab dist
    iMOAB_AppID cplLnd2PID    = &cplLnd2AppID;     // land on coupler PEs, MCT dist
    int cpllnd = 9, cpllnd2 = 109;
    ierr = iMOAB_RegisterApplication( "LNDX", &dup_comm_world, &cpllnd,
                                              cplLndPID );  // mesh on coupler pes
    CHECKIERR( ierr, "Cannot register LNDX over coupler PEs" )
    ierr = iMOAB_RegisterApplication( "LNDX2", &dup_comm_world, &cpllnd2,
                                                  cplLnd2PID );  // mesh on coupler pes, MCT distr, point cloud
    CHECKIERR( ierr, "Cannot register LNDX2 over coupler PEs" )

    if (entType == 1)
    {
        ierr = iMOAB_LoadMesh( cplLndPID, lndMoab.c_str(), readopts.c_str(), &nghlay ); // moab mesh can be cells
    }
    else // entType == 0
    {
        ierr = iMOAB_LoadMesh( cplLndPID, lndMoab.c_str(), readoptsLnd.c_str(), &nghlay ); // moab mesh is point cloud
    }

    CHECKIERR( ierr, "Cannot load moab mesh on coupler pes" )

    ierr = iMOAB_LoadMesh( cplLnd2PID, lndMct.c_str(), readoptsLnd.c_str(), &nghlay );
    CHECKIERR( ierr, "Cannot load mct point cloud on coupler pes" )

    int nverts[3], nelem[3];
    int tagType[2] = {  DENSE_DOUBLE, DENSE_INTEGER};
    int sizeTag = 1;
    int tagIndex = -1;
	/*
	 * Each process in the communicator will have access to a local mesh instance, which will
	 * contain the original cells in the local partition and ghost entities. Number of vertices,
	 * primary cells, visible blocks, number of sidesets and nodesets boundary conditions will be
	 * returned in size 3 arrays, for local, ghost and total numbers.
	 */
	ierr = iMOAB_GetMeshInfo( cplLnd2PID, nverts, nelem, 0, 0, 0 );
	CHECKIERR( ierr, "Cannot get info on mct mesh on coupler pes" )
	std::vector<int> gids(nverts[0]);
	int nvals = nverts[0]*ntags;
	std::vector<double> fracts(nvals);
	ierr = iMOAB_DefineTagStorage( cplLnd2PID, tagname.c_str(), &tagType[0], &sizeTag, &tagIndex );
	CHECKIERR( ierr, "Cannot define frac tag" )

	ierr = iMOAB_DefineTagStorage( cplLnd2PID, "GLOBAL_ID", &tagType[1], &sizeTag, &tagIndex );
    CHECKIERR( ierr, "Cannot define gid tag" )

	int mctEntType = 0; // vertex
	ierr  = iMOAB_GetDoubleTagStorage( cplLnd2PID, tagname.c_str(), &nvals, &mctEntType, &fracts[0] );
	CHECKIERR( ierr, "Cannot get frac tag on lnd2 on coupler pes" )

	ierr  = iMOAB_GetIntTagStorage( cplLnd2PID, "GLOBAL_ID", &nverts[0], &mctEntType, &gids[0] );
    CHECKIERR( ierr, "Cannot get global id tag on lnd2 on coupler pes" )

	// now on moab mesh
	ierr = iMOAB_DefineTagStorage( cplLndPID, tagname.c_str(), &tagType[0], &sizeTag, &tagIndex );
    CHECKIERR( ierr, "Cannot define list of tags" )

	ierr = iMOAB_SetDoubleTagStorageWithGid(cplLndPID, tagname.c_str(), &nvals, &entType, &fracts[0], &gids[0] );
    CHECKIERR( ierr, "Cannot set double tag with global ids " )

    char fileWriteOptions[] = "PARALLEL=WRITE_PART";
	char outputFileLnd[] = "recvLnd4.h5m";
	ierr                 = iMOAB_WriteMesh( cplLndPID, outputFileLnd, fileWriteOptions );
	CHECKIERR( ierr, "cannot write lnd mesh after setting tag" )

    ierr = iMOAB_DeregisterApplication( cplLnd2PID );
    CHECKIERR( ierr, "cannot deregister app LNDX2" )
    ierr = iMOAB_DeregisterApplication( cplLndPID );
    CHECKIERR( ierr, "cannot deregister app LNDX" )

    ierr = iMOAB_Finalize();
    CHECKIERR( ierr, "did not finalize iMOAB" )

    MPI_Comm_free(&dup_comm_world );
    MPI_Finalize();
    return 0;
}
