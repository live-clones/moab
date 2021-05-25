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

    std::string atmFilename =
            "../../sandbox/MeshFiles/e3sm/ne4pg2_o240/ne4pg2_p8.h5m";
    return 0;
}
