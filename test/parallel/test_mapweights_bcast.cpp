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

#include <netcdf.h>
#include <netcdf_par.h>

using namespace moab;

const static double radius               = 1.0;
const double MOAB_PI                     = 3.1415926535897932384626433832795028841971693993751058209749445923;
const static double surface_area         = 4.0 * MOAB_PI * radius * radius;
std::string inMapFilename = "outCS5ICOD5_map.nc";

void test_tempest_create_meshes();
void test_tempest_map_bcast();
void read_buffered_map();

int main( int argc, char** argv )
{
    if( argc > 1 )
        inMapFilename = TestDir + "/" + std::string( argv[1] );
    else
        inMapFilename = TestDir + "/" + inMapFilename;

    // REGISTER_TEST( test_tempest_create_meshes );
    REGISTER_TEST( test_tempest_map_bcast );
    REGISTER_TEST( read_buffered_map );

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
    const static std::string outFilenames[2] = { "outTempestCS.g", "outTempestRLL.g" };

    std::cout << "Creating TempestRemap Cubed-Sphere Mesh ...\n";
    Mesh src_tempest_mesh, dest_tempest_mesh;
    ierr = GenerateCSMesh( src_tempest_mesh, blockSize, outFilenames[0], "NetCDF4" );
    CHECK_EQUAL( ierr, 0 );
    ierr = GenerateICOMesh( dest_tempest_mesh, blockSize, computeDual, outFilenames[1], "NetCDF4" );
    CHECK_EQUAL( ierr, 0 );

    // Compute the surface area of CS mesh
    const double ssphere_area = src_tempest_mesh.CalculateFaceAreas( false );
    CHECK_REAL_EQUAL( ssphere_area, surface_area, 1e-10 );
    const double tsphere_area = dest_tempest_mesh.CalculateFaceAreas( false );
    CHECK_REAL_EQUAL( tsphere_area, surface_area, 1e-10 );
}

#define NCERROREXIT( err, MSG )                                                     \
    if( err )                                                                       \
    {                                                                               \
        std::cout << "Error (" << err << "): " << MSG;                              \
        std::cout << "\nAborting with message: " << nc_strerror( err ) << "... \n"; \
        return;                                                                     \
    }

// void read_nc_par(const std::string& filename, MPI_Comm& comm)
// {
//     int err;
//     int ncid;
//     int na_id, nb_id, ns_id;
//     int n_a, n_b, n_s;
//     MPI_Info info = MPI_INFO_NULL;
//     int omode = NC_NOWRITE | NC_PNETCDF | NC_MPIIO;
//     err           = nc_open_par( filename.c_str(), omode, comm, info, &ncid );NCERROREXIT( err, "reading the NetCDF Map file: " << filename );

//     /* free info object */
//     if( info != MPI_INFO_NULL ) MPI_Info_free( &info );

//     /* inquire dimension IDs and lengths */
//     err = nc_inq_dimid( ncid, "n_a", &na_id );NCERROREXIT( err, "reading the n_s dimension" );
//     err = nc_inq_dimid( ncid, "n_b", &nb_id );NCERROREXIT( err, "reading the n_s dimension" );
//     err = nc_inq_dimid( ncid, "n_s", &ns_id );NCERROREXIT( err, "reading the n_s dimension" );

//     err = nc_inq_dimlen( ncid, na_id, &n_a );NCERROREXIT( err, "reading the n_s value" );
//     err = nc_inq_dimlen( ncid, nb_id, &n_b );NCERROREXIT( err, "reading the n_s value" );
//     err = nc_inq_dimlen( ncid, ns_id, &n_s );NCERROREXIT( err, "reading the n_s value" );

//     if (verbose)
//         printf( "global dimension of map: %zd X %zd, with %zd non-zero elements.\n", n_a, n_b, n_s );

//     err = nc_close( ncid );NCERROREXIT( err, "closing file failed" );

//     return;
// }

void read_buffered_map( )
{
    NcError error( NcError::silent_nonfatal );

    // Define our total buffer size to use
    NcFile* ncMap = NULL;
    NcVar *varRow = NULL, *varCol = NULL, *varS = NULL;
    const int nBufferSize = 1 * 1024 * 1024;  // 64 KB
    int nA = 0, nB = 0, nVA = 0, nVB = 0, nS = 0;

    MPI_Comm commW                = MPI_COMM_WORLD;
    const int rootProc            = 0;
    // const std::string outFilename = TestDir + "/" + outFilenames[0];
    double dTolerance             = 1e-08;
    int mpierr                    = 0;

    int nprocs, rank;
    MPI_Comm_size( commW, &nprocs );
    MPI_Comm_rank( commW, &rank );

    if (rank == 0)
    {
        ncMap = new NcFile( inMapFilename.c_str(), NcFile::ReadOnly );
        if( !ncMap->is_valid() ) { _EXCEPTION1( "Unable to open input map file \"%s\"", inMapFilename.c_str() ); }

        // Source and Target mesh resolutions
        NcDim* dimNA = ncMap->get_dim( "n_a" );
        if( dimNA == NULL ) { _EXCEPTIONT( "Input map missing dimension \"n_a\"" ); }

        NcDim* dimNB = ncMap->get_dim( "n_b" );
        if( dimNB == NULL ) { _EXCEPTIONT( "Input map missing dimension \"n_b\"" ); }

        nA = dimNA->size();
        nB = dimNB->size();

        NcDim* dimNVA = ncMap->get_dim( "nv_a" );
        if( dimNA == NULL ) { _EXCEPTIONT( "Input map missing dimension \"nv_a\"" ); }

        NcDim* dimNVB = ncMap->get_dim( "nv_b" );
        if( dimNB == NULL ) { _EXCEPTIONT( "Input map missing dimension \"nv_b\"" ); }

        nVA = dimNVA->size();
        nVB = dimNVB->size();

        // Read SparseMatrix entries
        NcDim* dimNS = ncMap->get_dim( "n_s" );
        if( dimNS == NULL ) { _EXCEPTION1( "Map file \"%s\" does not contain dimension \"n_s\"", inMapFilename.c_str() ); }

        // store total number of nonzeros
        nS = dimNS->size();

        varRow = ncMap->get_var( "row" );
        if( varRow == NULL ) { _EXCEPTION1( "Map file \"%s\" does not contain variable \"row\"", inMapFilename.c_str() ); }

        varCol = ncMap->get_var( "col" );
        if( varRow == NULL ) { _EXCEPTION1( "Map file \"%s\" does not contain variable \"col\"", inMapFilename.c_str() ); }

        varS = ncMap->get_var( "S" );
        if( varRow == NULL ) { _EXCEPTION1( "Map file \"%s\" does not contain variable \"S\"", inMapFilename.c_str() ); }

    }

    const int nNNZBytes      = 2 * sizeof( int ) + sizeof( double );
    const int nTotalBytes    = nS * nNNZBytes;
    int nBytesRemaining      = nS;
    int nBufferedReads = static_cast< int >( ceil( (1.0 * nTotalBytes) / (1.0 * nBufferSize) ) );
    if (rank == 0) printf( "nBufferedReads: nTotalBytes = %d, nBufferSize = %d, nBufferedReads = %d\n", nTotalBytes, nBufferSize, nBufferedReads);

    int runData[6] = { nA, nB, nVA, nVB, nS, nBufferedReads };

    mpierr = MPI_Bcast( runData, 6, MPI_INTEGER, rootProc, commW );
    CHECK_EQUAL( mpierr, 0 );

    if (rank)
    {
        nA             = runData[0];
        nB             = runData[1];
        nVA            = runData[2];
        nVB            = runData[3];
        nS             = runData[4];
        nBufferedReads = runData[5];
    }

    // Let us declare the map object for every process
    OfflineMap mapRemap;
    // SparseMatrix< double > sparseMatrix;
    SparseMatrix< double >& sparseMatrix = mapRemap.GetSparseMatrix();

    if( rank == 0 )
        printf( "Parameters: nNNZBytes = %d, nS = %d, nTotalBytes = %d, nBufferedReads = %d\n", nNNZBytes, nS,
                nTotalBytes, nBufferedReads );

    long offset = 0;
    for( int iRead = 0; iRead < nBufferedReads ; ++iRead)
    {
        /* split the rows and send to processes in chunks */
        std::vector< int > rowAllocation;
        std::vector< int > sendCounts;
        std::vector< int > displs;
        const int NDATA          = 3;
        int nMapGlobalRowCols[3] = { 0, 0, 0 };

        // pointers to the data
        DataArray1D< int > vecRow;
        DataArray1D< int > vecCol;
        DataArray1D< double > vecS;
        if( rank == 0 )
        {
            int nLocSize = std::min( nBytesRemaining, static_cast< int >( nBufferSize / nNNZBytes ) );
            // int nLocSize    = nLocBufSize / nNNZBytes;

            printf( "Reading file: bytes %ld to %ld\n", offset, offset + nLocSize );

            vecRow.Allocate( nLocSize );
            vecCol.Allocate( nLocSize );
            vecS.Allocate( nLocSize );

            varRow->set_cur( (long)( offset * sizeof( int ) ) );
            varRow->get( &( vecRow[0] ), nLocSize );

            varCol->set_cur( (long)( offset * sizeof( int ) ) );
            varCol->get( &( vecCol[0] ), nLocSize );

            varS->set_cur( (long)( offset * sizeof( double ) ) );
            varS->get( &( vecS[0] ), nLocSize );

            // Decrement vecRow and vecCol
            std::vector<int> rowUn, colUn;
            for( size_t i = 0; i < vecRow.GetRows(); i++ )
            {
                vecRow[i]--;
                if( std::find( rowUn.begin(), rowUn.end(), vecRow[i] ) == rowUn.end() ) rowUn.push_back( vecRow[i] );
                vecCol[i]--;
                if( std::find( colUn.begin(), colUn.end(), vecRow[i] ) == colUn.end() ) colUn.push_back( vecRow[i] );
            }

            offset += nLocSize;
            nBytesRemaining -= nLocSize;

            // Set the entries of the map
            // sparseMatrix.SetEntries( vecRow, vecCol, vecS );

            nMapGlobalRowCols[0] = rowUn.size();
            nMapGlobalRowCols[1] = colUn.size();
            nMapGlobalRowCols[2] = static_cast< int >( vecS.GetRows() );
            printf( "Rank 0: nMapGlobalRowCols = %d, %d, %d, vecRow = %zu\n", nMapGlobalRowCols[0], nMapGlobalRowCols[1],
                    nMapGlobalRowCols[2], vecRow.GetRows() );
            rowUn.clear();
            colUn.clear();

            // Now let us perform the communication of the data
            // Form the dataset to send to other processes
            int* rg    = vecRow;

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
                rowAllocation[i * NDATA + 1] = std::find( rg, rg + vecS.GetRows(), sum ) - rg;
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
            //             rowAllocation[i * NDATA + 1], rowAllocation[i * NDATA + 2], i, displs[i], sendCounts[i]
            //             );
        }

        MPI_Barrier(commW);

        printf( "Rank %d: nMapGlobalRowCols = %d, %d, %d\n", rank, nMapGlobalRowCols[0], nMapGlobalRowCols[1],
                nMapGlobalRowCols[2] );
        mpierr = MPI_Bcast( nMapGlobalRowCols, 3, MPI_INTEGER, rootProc, commW );
        CHECK_EQUAL( mpierr, 0 );

        int allocationSize[NDATA] = { 0, 0, 0 };
        mpierr = MPI_Scatter( rowAllocation.data(), NDATA, MPI_INTEGER, &allocationSize, NDATA, MPI_INTEGER,
                                rootProc, commW );
        CHECK_EQUAL( mpierr, 0 );

        // create buffers to receive the data from root process
        DataArray1D< int > dataRows, dataCols;
        DataArray1D< double > dataEntries;
        dataRows.Allocate( allocationSize[2] );
        dataCols.Allocate( allocationSize[2] );
        dataEntries.Allocate( allocationSize[2] );

        /* TODO: combine row and column scatters together. Save one communication by better packing */
        // Scatter the rows, cols and entries to the processes according to the partition
        mpierr = MPI_Scatterv( vecRow, sendCounts.data(), displs.data(), MPI_INTEGER, (int*)( dataRows ),
                                dataRows.GetByteSize(), MPI_INTEGER, rootProc, commW );
        CHECK_EQUAL( mpierr, 0 );

        mpierr = MPI_Scatterv( vecCol, sendCounts.data(), displs.data(), MPI_INTEGER, (int*)( dataCols ),
                                dataCols.GetByteSize(), MPI_INTEGER, rootProc, commW );
        CHECK_EQUAL( mpierr, 0 );

        mpierr = MPI_Scatterv( vecS, sendCounts.data(), displs.data(), MPI_DOUBLE, (double*)( dataEntries ),
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
            //     std::cout << i << " -- row,col,val = " << rindex << ", " << rowMap[rindex] << ", " << cindex <<
            //     ", "
            //     << colMap[cindex] << ", " << dataEntries[i]
            //               << std::endl;
        }

        // Now set the data received from root onto the sparse matrix
        sparseMatrix.SetEntries( dataRows, dataCols, dataEntries );

        printf( "Rank %d: sparsematrix size = %d x %d\n", rank, allocationSize[0], allocationSize[2] );
    }

    assert( nBytesRemaining == 0 );

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

    delete ncMap;

    return;
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
    const std::string outFilename = inMapFilename;
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
