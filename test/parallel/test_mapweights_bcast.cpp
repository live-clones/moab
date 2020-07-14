/**
 * @file Read and distribute mapping weights to multiple processes
 * @author Vijay Mahadevan (mahadevan@anl.gov)
 * @brief
 * @version 0.1
 * @date 2020-07-13
 *
 */

#include <iostream>
#include <vector>
#include "moab/Core.hpp"

#ifndef MOAB_HAVE_TEMPESTREMAP
#error Need TempestRemap dependency for this test
#endif

#include "OfflineMap.h"
#include "moab/Remapping/TempestRemapper.hpp"
#include "moab/ParallelComm.hpp"

#ifndef IS_BUILDING_MB
#define IS_BUILDING_MB
#endif

#include "Internals.hpp"
#include "TestRunner.hpp"

using namespace moab;

const static double radius               = 1.0;
const double MOAB_PI                     = 3.1415926535897932384626433832795028841971693993751058209749445923;
const static double surface_area         = 4.0 * MOAB_PI * radius * radius;
const static std::string outFilenames[3] = { "outCS5ICOD5_map.nc", "outTempestCS.g", "outTempestRLL.g" };

void test_tempest_create_meshes();
void test_tempest_map_bcast();

int main( int argc, char** argv )
{
    // REGISTER_TEST( test_tempest_create_meshes );
    REGISTER_TEST( test_tempest_map_bcast );

#ifdef MOAB_HAVE_MPI
    MPI_Init( &argc, &argv );
#endif
    int result = RUN_TESTS( argc, argv );
#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif
    return result;
}

void test_tempest_create_meshes()
{
    NcError error( NcError::verbose_nonfatal );
    const int blockSize    = 5;
    const bool computeDual = true;
    int ierr;

    std::cout << "Creating TempestRemap Cubed-Sphere Mesh ...\n";
    Mesh src_tempest_mesh, dest_tempest_mesh;
    ierr = GenerateCSMesh( src_tempest_mesh, blockSize, outFilenames[1], "NetCDF4" );
    CHECK_EQUAL( ierr, 0 );
    ierr = GenerateICOMesh( dest_tempest_mesh, blockSize, computeDual, outFilenames[2], "NetCDF4" );
    CHECK_EQUAL( ierr, 0 );

    // Compute the surface area of CS mesh
    const double ssphere_area = src_tempest_mesh.CalculateFaceAreas( false );
    CHECK_REAL_EQUAL( ssphere_area, surface_area, 1e-10 );
    const double tsphere_area = dest_tempest_mesh.CalculateFaceAreas( false );
    CHECK_REAL_EQUAL( tsphere_area, surface_area, 1e-10 );
}

void test_tempest_map_bcast()
{
    /*
     *
     * Load mapping file generated with CS5 and ICOD5 meshes with FV-FV
     * projection for field data. It has following attributes:
     *
     *      Analyzing map
     *      ..Per-dof consistency  (tol 1.00000e-08).. PASS
     *      ..Per-dof conservation (tol 1.00000e-08).. PASS
     *      ..Weights within range [-10,+10].. Done
     *      ..
     *      ..  Total nonzero entries: 871
     *      ..   Column index min/max: 1 / 150 (150 source dofs)
     *      ..      Row index min/max: 1 / 252 (252 target dofs)
     *      ..      Source area / 4pi: 9.999999999999986e-01
     *      ..    Source area min/max: 7.758193626682276e-02 / 9.789674024213996e-02
     *      ..      Target area / 4pi: 9.999999999999988e-01
     *      ..    Target area min/max: 3.418848585325412e-02 / 5.362041648788487e-02
     *      ..    Source frac min/max: 9.999999999999984e-01 / 9.999999999999997e-01
     *      ..    Target frac min/max: 9.999999999999981e-01 / 1.000000000000000e+00
     *      ..    Map weights min/max: 2.377132934153108e-07 / 9.999999999999988e-01
     *      ..       Row sums min/max: 9.999999999999981e-01 / 1.000000000000000e+00
     *      ..   Consist. err min/max: -1.887379141862766e-15 / 0.000000000000000e+00
     *      ..    Col wt.sums min/max: 9.999999999999983e-01 / 9.999999999999997e-01
     *      ..   Conserv. err min/max: -1.665334536937735e-15 / -3.330669073875470e-16
     *      ..
     *      ..Histogram of nonzero entries in sparse matrix
     *      ....Column 1: Number of nonzero entries (bin minimum)
     *      ....Column 2: Number of columns with that many nonzero values
     *      ....Column 3: Number of rows with that many nonzero values
     *      ..[[1, 0, 2],[2, 0, 51],[3, 0, 48],[4, 15, 138],[5, 39, 7],[6, 64, 6],[7, 24, 0],[8, 8, 0]]
     *      ..
     *      ..Histogram of weights
     *      ....Column 1: Lower bound on weights
     *      ....Column 2: Upper bound on weights
     *      ....Column 3: # of weights in that bin
     *      ..[[0, 1, 871]]
     *      ------------------------------------------------------------
     *      0  Writing remap weights with size [252 X 150] and NNZ = 871
     *
     */
    NcError error( NcError::silent_nonfatal );
    MPI_Comm commW                = MPI_COMM_WORLD;
    const int rootProc            = 0;
    const std::string outFilename = TestDir + "/" + outFilenames[0];
    double dTolerance             = 1e-08;
    int mpierr                    = 0;

    int nprocs, rank;
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    const double dFillValueOverride = 0.0;
    OfflineMap mapRemap;
    if( rank == 0 )
    {
        std::cout << "Reading file: " << outFilename << std::endl;
        mapRemap.Read( outFilename );
        mapRemap.SetFillValueOverrideDbl( dFillValueOverride );
        mapRemap.SetFillValueOverride( static_cast< float >( dFillValueOverride ) );

        int isConsistent = mapRemap.IsConsistent( dTolerance );
        CHECK_EQUAL( isConsistent, 0 );
        int isConservative = mapRemap.IsConservative( dTolerance );
        CHECK_EQUAL( isConservative, 0 );
        int isMonotone = mapRemap.IsMonotone( dTolerance );
        CHECK_EQUAL( isMonotone, 4600 );

        // verify that the source and target mesh areas are equal
        DataArray1D< double >& srcAreas = mapRemap.GetSourceAreas();
        DataArray1D< double >& tgtAreas = mapRemap.GetTargetAreas();

        double totSrcArea = 0.0;
        for( size_t i = 0; i < srcAreas.GetRows(); ++i )
            totSrcArea += srcAreas[i];

        double totTgtArea = 0.0;
        for( size_t i = 0; i < tgtAreas.GetRows(); ++i )
            totTgtArea += tgtAreas[i];

        CHECK_REAL_EQUAL( totSrcArea, totTgtArea, 1e-10 );
    }

    MPI_Barrier( commW );

    SparseMatrix< double >& sparseMatrix = mapRemap.GetSparseMatrix();

    // Get the row, col and entry data tupl
    // Retrieve all data from the map operator on root process
    DataArray1D< int > dataRowsGlobal, dataColsGlobal;
    DataArray1D< double > dataEntriesGlobal;
    if( rank == 0 ) sparseMatrix.GetEntries( dataRowsGlobal, dataColsGlobal, dataEntriesGlobal );

    int *rg, *cg;
    double *de;

    int nMapGlobalRowCols[3] = { sparseMatrix.GetRows(), sparseMatrix.GetColumns(),
                                 static_cast< int >( dataEntriesGlobal.GetRows() ) };

    /* split the rows and send to processes in chunks */
    std::vector< int > rowAllocation;
    std::vector< int > sendCounts;
    std::vector< int > displs;
    const int NDATA = 3;
    if( rank == 0 )
    {
        rg = dataRowsGlobal;
        cg = dataColsGlobal;
        de = dataEntriesGlobal;

        // Let us do a trivial partitioning of the row data
        rowAllocation.resize( nprocs * NDATA );
        displs.resize( nprocs );
        sendCounts.resize(nprocs);

        int nRowPerPart = nMapGlobalRowCols[0] / nprocs;
        int remainder = nMapGlobalRowCols[0] % nprocs;  // Keep the remainder in root
        // printf( "nRowPerPart = %d, remainder = %d\n", nRowPerPart, remainder );
        int sum = 0;
        for( int i = 0; i < nprocs; ++i )
        {
            rowAllocation[i * NDATA] = nRowPerPart + ( i == 0 ? remainder : 0 );
            displs[i] = sum;
            sum += rowAllocation[i * NDATA];
            rowAllocation[i * NDATA + 1] = std::find( rg, rg + dataRowsGlobal.GetRows(), displs[i] ) - rg;
        }
        for( int i = 0; i < nprocs-1; ++i )
        {
            rowAllocation[i * NDATA + 2] = rowAllocation[( i + 1 ) * NDATA + 1] - rowAllocation[i * NDATA + 1];
            sendCounts[i]                = rowAllocation[i * NDATA + 2];
        }
        rowAllocation[(nprocs-1) * NDATA + 2] = nMapGlobalRowCols[2] - displs[nprocs - 1];
        sendCounts[nprocs-1]                  = rowAllocation[( nprocs - 1 ) * NDATA + 2];

        // for( int i = 0; i < nprocs; i++ )
        //     printf( "sendcounts[%d] = %d, %d, %d\tdispls[%d] = %d\n", i, rowAllocation[i * NDATA],
        //             rowAllocation[i * NDATA + 1], rowAllocation[i * NDATA + 2], i, displs[i] );
    }

    // std::cout << nMapGlobalRowCols[0] << ", " << nMapGlobalRowCols[1] << ", " << nMapGlobalRowCols[2] << "\n";
    mpierr = MPI_Bcast( nMapGlobalRowCols, 3, MPI_INTEGER, rootProc, commW );
    CHECK_EQUAL( mpierr, 0 );

    int allocationSize[NDATA] = { 0, 0, 0 };
    mpierr =
        MPI_Scatter( rowAllocation.data(), NDATA, MPI_INTEGER, &allocationSize, NDATA, MPI_INTEGER, rootProc, commW );
    CHECK_EQUAL( mpierr, 0 );

    // create buffers to receive the data from root process
    DataArray1D< int > dataRows, dataCols;
    DataArray1D< double > dataEntries;
    dataRows.Allocate( allocationSize[2] );
    dataCols.Allocate( allocationSize[2] );
    dataEntries.Allocate( allocationSize[2] );

    // Scatter the rows, cols and entries to the processes according to the partition
    mpierr = MPI_Scatterv( rg, sendCounts.data(), displs.data(), MPI_INTEGER, (int*)( dataRows ),
                           dataRows.GetByteSize(), MPI_INTEGER, rootProc, commW );
    CHECK_EQUAL( mpierr, 0 );

    mpierr = MPI_Scatterv( cg, sendCounts.data(), displs.data(), MPI_INTEGER, (int*)( dataCols ),
                           dataCols.GetByteSize(), MPI_INTEGER, rootProc, commW );
    CHECK_EQUAL( mpierr, 0 );

    mpierr = MPI_Scatterv( de, sendCounts.data(), displs.data(), MPI_INTEGER,
                           (double*)( dataEntries ), dataEntries.GetByteSize(), MPI_INTEGER, rootProc, commW );
    CHECK_EQUAL( mpierr, 0 );

    // Now set the data received from root onto the sparse matrix
    sparseMatrix.SetEntries( dataRows, dataCols, dataEntries );

    // Let us run some more tests to verify that the map data has been distributed correctly
}
