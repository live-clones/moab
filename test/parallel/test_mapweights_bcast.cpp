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
     *      ..  Total nonzero entries: 8976
     *      ..   Column index min/max: 1 / 150 (150 source dofs)
     *      ..      Row index min/max: 1 / 252 (252 target dofs)
     *      ..      Source area / 4pi: 9.999999999991666e-01
     *      ..    Source area min/max: 7.758193626656797e-02 / 9.789674024202324e-02
     *      ..      Target area / 4pi: 9.999999999999928e-01
     *      ..    Target area min/max: 3.418848585325412e-02 / 5.362041648788460e-02
     *      ..    Source frac min/max: 9.999999999997851e-01 / 1.000000000003865e+00
     *      ..    Target frac min/max: 9.999999999993099e-01 / 1.000000000000784e+00
     *      ..    Map weights min/max: -2.062659472710135e-01 / 1.143815352351317e+00
     *      ..       Row sums min/max: 9.999999999993099e-01 / 1.000000000000784e+00
     *      ..   Consist. err min/max: -6.901146321069973e-13 / 7.838174553853605e-13
     *      ..    Col wt.sums min/max: 9.999999999997854e-01 / 1.000000000003865e+00
     *      ..   Conserv. err min/max: -2.146061106600428e-13 / 3.865352482534945e-12
     *      ..
     *      ..Histogram of nonzero entries in sparse matrix
     *      ....Column 1: Number of nonzero entries (bin minimum)
     *      ....Column 2: Number of columns with that many nonzero values
     *      ....Column 3: Number of rows with that many nonzero values
     *      ..[[25, 0, 2],[29, 0, 8],[30, 0, 10],[31, 150, 232]]
     *      ..
     *      ..Histogram of weights
     *      ....Column 1: Lower bound on weights
     *      ....Column 2: Upper bound on weights
     *      ....Column 3: # of weights in that bin
     *      ..[[-1, 0, 4577],[0, 1, 4365],[1, 2, 34]]
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
        CHECK_EQUAL( 0, isConsistent );
        int isConservative = mapRemap.IsConservative( dTolerance );
        CHECK_EQUAL( 0, isConservative );
        int isMonotone = mapRemap.IsMonotone( dTolerance );
        CHECK_EQUAL( 4600, isMonotone );

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
    if( rank == 0 )
    {
        sparseMatrix.GetEntries( dataRowsGlobal, dataColsGlobal, dataEntriesGlobal );
        std::cout << "Total NNZ in remap matrix = " << dataRowsGlobal.GetRows() << std::endl;
    }

    int *rg, *cg;
    double* de;
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
        sendCounts.resize( nprocs );

        int nRowPerPart = nMapGlobalRowCols[0] / nprocs;
        int remainder   = nMapGlobalRowCols[0] % nprocs;  // Keep the remainder in root
        // printf( "nRowPerPart = %d, remainder = %d\n", nRowPerPart, remainder );
        int sum = 0;
        for( int i = 0; i < nprocs; ++i )
        {
            rowAllocation[i * NDATA]     = nRowPerPart + ( i == 0 ? remainder : 0 );
            rowAllocation[i * NDATA + 1] = std::find( rg, rg + dataRowsGlobal.GetRows(), sum ) - rg;
            sum += rowAllocation[i * NDATA];
            displs[i] = rowAllocation[i * NDATA + 1];
        }
        for( int i = 0; i < nprocs - 1; ++i )
        {
            rowAllocation[i * NDATA + 2] = rowAllocation[( i + 1 ) * NDATA + 1] - rowAllocation[i * NDATA + 1];
            sendCounts[i]                = rowAllocation[i * NDATA + 2];
        }
        rowAllocation[( nprocs - 1 ) * NDATA + 2] = nMapGlobalRowCols[2] - displs[nprocs - 1];
        sendCounts[nprocs - 1]                    = rowAllocation[( nprocs - 1 ) * NDATA + 2];

        // for( int i = 0; i < nprocs; i++ )
        //     printf( "sendcounts[%d] = %d, %d, %d\tdispls[%d] = %d, %d\n", i, rowAllocation[i * NDATA],
        //             rowAllocation[i * NDATA + 1], rowAllocation[i * NDATA + 2], i, displs[i], sendCounts[i] );
    }

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

    /* TODO: combine row and column scatters together. Save one communication by better packing */
    // Scatter the rows, cols and entries to the processes according to the partition
    mpierr = MPI_Scatterv( rg, sendCounts.data(), displs.data(), MPI_INTEGER, (int*)( dataRows ),
                           dataRows.GetByteSize(), MPI_INTEGER, rootProc, commW );
    CHECK_EQUAL( mpierr, 0 );

    mpierr = MPI_Scatterv( cg, sendCounts.data(), displs.data(), MPI_INTEGER, (int*)( dataCols ),
                           dataCols.GetByteSize(), MPI_INTEGER, rootProc, commW );
    CHECK_EQUAL( mpierr, 0 );

    mpierr = MPI_Scatterv( de, sendCounts.data(), displs.data(), MPI_DOUBLE, (double*)( dataEntries ),
                           dataEntries.GetByteSize(), MPI_DOUBLE, rootProc, commW );
    CHECK_EQUAL( mpierr, 0 );

    // Compute an offset for the rows and columns by creating a local to global mapping
    std::vector< int > rowMap, colMap;
    int rindex = 0, cindex = 0;
    for( int i = 0; i < allocationSize[2]; ++i )
    {
        std::vector< int >::iterator riter = std::find( rowMap.begin(), rowMap.end(), dataRows[i] );
        if( riter == rowMap.end() )
        {
            rowMap.push_back( dataRows[i] );
            rindex = rowMap.size() - 1;
        }
        else
            rindex = riter - rowMap.begin();
        dataRows[i] = rindex;

        std::vector< int >::iterator citer = std::find( colMap.begin(), colMap.end(), dataCols[i] );
        if( citer == colMap.end() )
        {
            colMap.push_back( dataCols[i] );
            cindex = colMap.size() - 1;
        }
        else
            cindex = citer - colMap.begin();
        dataCols[i] = cindex;

        // if (rank == 1)
        //     std::cout << i << " -- row,col,val = " << rindex << ", " << rowMap[rindex] << ", " << cindex << ", " <<
        //     colMap[cindex] << ", " << dataEntries[i]
        //               << std::endl;
    }

    // Now set the data received from root onto the sparse matrix
    sparseMatrix.SetEntries( dataRows, dataCols, dataEntries );

    /**
     * Perform a series of checks to ensure that the local maps in each process still preserve map properties
     * Let us run our usual tests to verify that the map data has been distributed correctly
     *   1) Consistency check -- row sum = 1.0
     *   2) Disable conservation checks that require global column sums (and hence more communication/reduction)
     *   3) Perform monotonicity check and ensure failures globally are same as serial case
     */

    // consistency: row sums
    int isConsistentP = mapRemap.IsConsistent( dTolerance );
    CHECK_EQUAL( 0, isConsistentP );

    // conservation: we will disable this conservation check for now.
    // int isConservativeP = mapRemap.IsConservative( dTolerance );
    // CHECK_EQUAL( 0, isConservativeP );

    // monotonicity: local failures
    int isMonotoneP = mapRemap.IsMonotone( dTolerance );
    // Accumulate sum of failures to ensure that it is exactly same as serial case
    int isMonotoneG;
    mpierr = MPI_Allreduce( &isMonotoneP, &isMonotoneG, 1, MPI_INTEGER, MPI_SUM, commW );
    CHECK_EQUAL( mpierr, 0 );
    // 4600 fails in serial. What about parallel ? Check and confirm.
    CHECK_EQUAL( 4600, isMonotoneG );

    return;
}
