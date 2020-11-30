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

const static double radius       = 1.0;
const double MOAB_PI             = 3.1415926535897932384626433832795028841971693993751058209749445923;
const static double surface_area = 4.0 * MOAB_PI * radius * radius;
std::string inMapFilename        = "outCS5ICOD5_map.nc";

void test_tempest_create_meshes();
void test_tempest_map_bcast();
void read_buffered_map();

int main( int argc, char** argv )
{
    if( argc > 1 )
        inMapFilename = std::string( argv[1] );
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


void read_buffered_map()
{
    NcError error( NcError::silent_nonfatal );

    // Define our total buffer size to use
    NcFile* ncMap = NULL;
    NcVar *varRow = NULL, *varCol = NULL, *varS = NULL;

#ifdef MOAB_HAVE_MPI
    int mpierr = 0;
#endif
    // Main inputs;
    //   1. Map Filename
    //   2. Ownership info from the processes to be gathered on root
    const int nBufferSize = 64 * 1024;  // 64 KB
    const int nNNZBytes   = 2 * sizeof( int ) + sizeof( double );
    const int nMaxEntries = nBufferSize / nNNZBytes;
    int nA = 0, nB = 0, nVA = 0, nVB = 0, nS = 0;
    double dTolerance  = 1e-08;

    const int rootProc = 0;
    int nprocs = 1, rank = 0;
#ifdef MOAB_HAVE_MPI
    MPI_Comm commW     = MPI_COMM_WORLD;
    MPI_Comm_size( commW, &nprocs );
    MPI_Comm_rank( commW, &rank );
#endif

    if( rank == rootProc )
    {
        ncMap = new NcFile( inMapFilename.c_str(), NcFile::ReadOnly );
        if( !ncMap->is_valid() ) { _EXCEPTION1( "Unable to open input map file \"%s\"", inMapFilename.c_str() ); }

        // Source and Target mesh resolutions
        NcDim* dimNA = ncMap->get_dim( "n_a" );
        if( dimNA == NULL ) { _EXCEPTIONT( "Input map missing dimension \"n_a\"" ); }
        else nA = dimNA->size();

        NcDim* dimNB = ncMap->get_dim( "n_b" );
        if( dimNB == NULL ) { _EXCEPTIONT( "Input map missing dimension \"n_b\"" ); }
        else nB = dimNB->size();

        NcDim* dimNVA = ncMap->get_dim( "nv_a" );
        if( dimNA == NULL ) { _EXCEPTIONT( "Input map missing dimension \"nv_a\"" ); }
        else nVA = dimNVA->size();

        NcDim* dimNVB = ncMap->get_dim( "nv_b" );
        if( dimNB == NULL ) { _EXCEPTIONT( "Input map missing dimension \"nv_b\"" ); }
        else nVB = dimNVB->size();

        // Read SparseMatrix entries
        NcDim* dimNS = ncMap->get_dim( "n_s" );
        if( dimNS == NULL )
        { _EXCEPTION1( "Map file \"%s\" does not contain dimension \"n_s\"", inMapFilename.c_str() ); }

        // store total number of nonzeros
        nS = dimNS->size();

        varRow = ncMap->get_var( "row" );
        if( varRow == NULL )
        { _EXCEPTION1( "Map file \"%s\" does not contain variable \"row\"", inMapFilename.c_str() ); }

        varCol = ncMap->get_var( "col" );
        if( varRow == NULL )
        { _EXCEPTION1( "Map file \"%s\" does not contain variable \"col\"", inMapFilename.c_str() ); }

        varS = ncMap->get_var( "S" );
        if( varRow == NULL )
        { _EXCEPTION1( "Map file \"%s\" does not contain variable \"S\"", inMapFilename.c_str() ); }
    }

    // const int nTotalBytes = nS * nNNZBytes;
    int nEntriesRemaining = nS;
    int nBufferedReads    = static_cast< int >( ceil( ( 1.0 * nS ) / ( 1.0 * nMaxEntries ) ) );
    int runData[6]        = { nA, nB, nVA, nVB, nS, nBufferedReads };

    mpierr = MPI_Bcast( runData, 6, MPI_INTEGER, rootProc, commW );
    CHECK_EQUAL( mpierr, 0 );

    if( rank != rootProc )
    {
        nA             = runData[0];
        nB             = runData[1];
        nVA            = runData[2];
        nVB            = runData[3];
        nS             = runData[4];
        nBufferedReads = runData[5];
    }
    else
        printf( "Global parameters: nA = %d, nB = %d, nS = %d\n", nA, nB, nS );

    std::vector< std::pair< int, int > > rowOwnership( nprocs );
    int nGRowPerPart   = nB / nprocs;
    int nGRowRemainder = nB % nprocs;  // Keep the remainder in root
    rowOwnership[0]    = std::make_pair( 0, nGRowPerPart + nGRowRemainder );
    int roffset        = rowOwnership[0].second;
    // printf("Rank %d ownership: %d -- %d\n", 0, rowOwnership[0].first, rowOwnership[0].second );
    for( int ip = 1; ip < nprocs; ++ip )
    {
        rowOwnership[ip] = std::make_pair( roffset, roffset + nGRowPerPart );
        roffset          = rowOwnership[ip].second;
        // printf( "Rank %d ownership: %d -- %d\n", ip, rowOwnership[ip].first, rowOwnership[ip].second );
    }

    // Let us declare the map object for every process
    OfflineMap mapRemap;
    SparseMatrix< double >& sparseMatrix = mapRemap.GetSparseMatrix();

    if( rank == rootProc )
        printf( "Parameters: nA=%d, nB=%d, nVA=%d, nVB=%d, nS=%d, nNNZBytes = %d, nBufferedReads = %d\n", nA, nB, nVA,
                nVB, nS, nNNZBytes, nBufferedReads );

    std::map< int, int > rowMap, colMap;
    int rindexMax = 0, cindexMax = 0;
    long offset = 0;

    /* Split the rows and send to processes in chunks
       Let us start the buffered read */
    for( int iRead = 0; iRead < nBufferedReads; ++iRead )
    {
        // pointers to the data
        DataArray1D< int > vecRow;
        DataArray1D< int > vecCol;
        DataArray1D< double > vecS;
        std::vector< std::vector< int > > dataPerProcess( nprocs );
        std::vector< int > nDataPerProcess( nprocs, 0 );
        for( int ip = 0; ip < nprocs; ++ip )
            dataPerProcess[ip].reserve( std::max( 100, nMaxEntries / nprocs ) );

        if( rank == rootProc )
        {
            int nLocSize = std::min( nEntriesRemaining, static_cast< int >( ceil( nBufferSize * 1.0 / nNNZBytes ) ) );

            printf( "Reading file: elements %ld to %ld\n", offset, offset + nLocSize );

            // Allocate and resize based on local buffer size
            vecRow.Allocate( nLocSize );
            vecCol.Allocate( nLocSize );
            vecS.Allocate( nLocSize );

            varRow->set_cur( (long)( offset ) );
            varRow->get( &( vecRow[0] ), nLocSize );

            varCol->set_cur( (long)( offset ) );
            varCol->get( &( vecCol[0] ), nLocSize );

            varS->set_cur( (long)( offset ) );
            varS->get( &( vecS[0] ), nLocSize );

            // Decrement vecRow and vecCol
            for( size_t i = 0; i < vecRow.GetRows(); i++ )
            {
                vecRow[i]--;
                vecCol[i]--;

                // TODO: Fix this linear search to compute process ownership
                int pOwner = -1;
                for( int ip = 0; ip < nprocs; ++ip )
                {
                    if( rowOwnership[ip].first <= vecRow[i] && rowOwnership[ip].second > vecRow[i] )
                    {
                        pOwner = ip;
                        break;
                    }
                }

                assert( pOwner >= 0 && pOwner < nprocs );
                dataPerProcess[pOwner].push_back( i );
            }

            offset += nLocSize;
            nEntriesRemaining -= nLocSize;

            for( int ip = 0; ip < nprocs; ++ip )
                nDataPerProcess[ip] = dataPerProcess[ip].size();
        }

        // Communicate the number of DoF data to expect from root
        int nEntriesComm = 0;
        mpierr = MPI_Scatter( nDataPerProcess.data(), 1, MPI_INT, &nEntriesComm, 1, MPI_INT, rootProc, commW );
        CHECK_EQUAL( mpierr, 0 );

        if( rank == rootProc )
        {
            std::vector< MPI_Request > cRequests;
            std::vector< MPI_Status > cStats( 2 * nprocs - 2 );
            cRequests.reserve( 2 * nprocs - 2 );

            // First send out all the relevant data buffers to other process
            std::vector< std::vector< int > > dataRowColsAll( nprocs );
            std::vector< std::vector< double > > dataEntriesAll( nprocs );
            for( int ip = 1; ip < nprocs; ++ip )
            {
                const int nDPP = nDataPerProcess[ip];
                if( nDPP > 0 )
                {
                    std::vector< int >& dataRowCols = dataRowColsAll[ip];
                    dataRowCols.resize( 2 * nDPP );
                    std::vector< double >& dataEntries = dataEntriesAll[ip];
                    dataEntries.resize( nDPP );

                    for( int ij = 0; ij < nDPP; ++ij )
                    {
                        dataRowCols[ij * 2]     = vecRow[dataPerProcess[ip][ij]];
                        dataRowCols[ij * 2 + 1] = vecCol[dataPerProcess[ip][ij]];
                        dataEntries[ij]         = vecS[dataPerProcess[ip][ij]];
                    }

                    MPI_Request rcsend, dsend;
                    mpierr = MPI_Isend( dataRowCols.data(), 2 * nDPP, MPI_INT, ip, iRead * 1000, commW, &rcsend );
                    CHECK_EQUAL( mpierr, 0 );
                    mpierr = MPI_Isend( dataEntries.data(), nDPP, MPI_DOUBLE, ip, iRead * 1000 + 1, commW, &dsend );
                    CHECK_EQUAL( mpierr, 0 );

                    cRequests.push_back( rcsend );
                    cRequests.push_back( dsend );
                }
            }

            // Next perform all necessary local work while we wait for the buffers to be sent out
            // Compute an offset for the rows and columns by creating a local to global mapping for rootProc
            int rindex = 0, cindex = 0;
            assert( dataPerProcess[0].size() - nEntriesComm == 0 );  // sanity check
            for( int i = 0; i < nEntriesComm; ++i )
            {
                const int& vecRowValue = vecRow[dataPerProcess[0][i]];
                const int& vecColValue = vecCol[dataPerProcess[0][i]];

                std::map< int, int >::iterator riter = rowMap.find( vecRowValue );
                if( riter == rowMap.end() )
                {
                    rowMap[vecRowValue] = rindexMax;
                    rindex              = rindexMax;
                    rindexMax++;
                }
                else
                    rindex = riter->second;

                std::map< int, int >::iterator citer = colMap.find( vecColValue );
                if( citer == colMap.end() )
                {
                    colMap[vecColValue] = cindexMax;
                    cindex              = cindexMax;
                    cindexMax++;
                }
                else
                    cindex = citer->second;

                sparseMatrix( rindex, cindex ) = vecS[i];
            }

            // Wait until all communication is pushed out
            mpierr = MPI_Waitall( cRequests.size(), cRequests.data(), cStats.data() );
            CHECK_EQUAL( mpierr, 0 );
        }  // if( rank == rootProc )
        else
        {
            if( nEntriesComm > 0 )
            {
                MPI_Request cRequests[2];
                MPI_Status cStats[2];
                /* code */
                // create buffers to receive the data from root process
                std::vector< int > dataRowCols( 2 * nEntriesComm );
                vecRow.Allocate( nEntriesComm );
                vecCol.Allocate( nEntriesComm );
                vecS.Allocate( nEntriesComm );

                /* TODO: combine row and column scatters together. Save one communication by better packing */
                // Scatter the rows, cols and entries to the processes according to the partition
                mpierr = MPI_Irecv( dataRowCols.data(), 2 * nEntriesComm, MPI_INT, rootProc, iRead * 1000, commW,
                                    &cRequests[0] );
                CHECK_EQUAL( mpierr, 0 );

                mpierr = MPI_Irecv( vecS, nEntriesComm, MPI_DOUBLE, rootProc, iRead * 1000 + 1, commW, &cRequests[1] );
                CHECK_EQUAL( mpierr, 0 );

                // Wait until all communication is pushed out
                mpierr = MPI_Waitall( 2, cRequests, cStats );
                CHECK_EQUAL( mpierr, 0 );

                // Compute an offset for the rows and columns by creating a local to global mapping
                int rindex = 0, cindex = 0;
                for( int i = 0; i < nEntriesComm; ++i )
                {
                    const int& vecRowValue = dataRowCols[i * 2];
                    const int& vecColValue = dataRowCols[i * 2 + 1];

                    std::map< int, int >::iterator riter = rowMap.find( vecRowValue );
                    if( riter == rowMap.end() )
                    {
                        rowMap[vecRowValue] = rindexMax;
                        rindex              = rindexMax;
                        rindexMax++;
                    }
                    else
                        rindex = riter->second;
                    vecRow[i] = rindex;

                    std::map< int, int >::iterator citer = colMap.find( vecColValue );
                    if( citer == colMap.end() )
                    {
                        colMap[vecColValue] = cindexMax;
                        cindex              = cindexMax;
                        cindexMax++;
                    }
                    else
                        cindex = citer->second;
                    vecCol[i] = cindex;

                    sparseMatrix( rindex, cindex ) = vecS[i];
                }

                // clear the local buffer
                dataRowCols.clear();
            }  // if( nEntriesComm > 0 )
        }      // if( rank != rootProc )

        MPI_Barrier( commW );
    }

    if( rank == rootProc )
    {
        assert( nEntriesRemaining == 0 );
        ncMap->close();
    }

    /**
     * Perform a series of checks to ensure that the local maps in each process still preserve map properties
     * Let us run our usual tests to verify that the map data has been distributed correctly
     *   1) Consistency check -- row sum = 1.0
     *   2) Disable conservation checks that require global column sums (and hence more communication/reduction)
     *   3) Perform monotonicity check and ensure failures globally are same as serial case
     */
    printf( "Rank %d: sparsematrix size = %d x %d\n", rank, sparseMatrix.GetRows(), sparseMatrix.GetColumns() );

    // consistency: row sums
    {
        int isConsistentP = mapRemap.IsConsistent( dTolerance );
        if( 0 != isConsistentP ) { printf( "Rank %d failed consistency checks...\n", rank ); }
        CHECK_EQUAL( 0, isConsistentP );
    }

    // conservation: we will disable this conservation check for now.
    // int isConservativeP = mapRemap.IsConservative( dTolerance );
    // CHECK_EQUAL( 0, isConservativeP );

    // monotonicity: local failures
    {
        int isMonotoneP = mapRemap.IsMonotone( dTolerance );
        // Accumulate sum of failures to ensure that it is exactly same as serial case
        int isMonotoneG;
        mpierr = MPI_Allreduce( &isMonotoneP, &isMonotoneG, 1, MPI_INTEGER, MPI_SUM, commW );
        CHECK_EQUAL( mpierr, 0 );

        // 4600 monotonicity fails in serial and parallel for the test map file
        CHECK_EQUAL( 4600, isMonotoneG );
    }

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
