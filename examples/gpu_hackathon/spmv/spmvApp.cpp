/** @example spmvApp.cpp
 * \brief This example reads two map files to project a field in forward and reverse
 * direction by computing projections via sparse Matrix-Vector products between source
 * component and target component. This is strictly an example to work on single
 * node and we do not care about MPI parallelism in this experiment.
 *
 * Usage:
 *      ./spmvApp srcRemapWeightFile -t -n iterations
 *
 * Note: Some datasets for forward and reverse maps have been uploaded to:
 *       https://ftp.mcs.anl.gov/pub/fathom/MeshFiles/maps/
 */

#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/ProgOptions.hpp"

#ifndef MOAB_HAVE_NETCDF
#error Require NetCDF to read the map files
#else
#include "netcdfcpp.h"
#endif

#ifndef MOAB_HAVE_EIGEN3
#error Require Eigen3 headers in order to apply SpMV routines
#else
#include <Eigen/Dense>
#include <Eigen/Sparse>
#endif

// C++ includes
#include <iostream>
#include <chrono>

using namespace moab;
using namespace std;

// Some global typedefs
typedef int MOABSInt;
typedef size_t MOABUInt;
typedef double MOABReal;

typedef std::chrono::high_resolution_clock Clock;
typedef std::chrono::high_resolution_clock::time_point Timer;
using std::chrono::duration_cast;

Timer start;
std::map< std::string, std::chrono::nanoseconds > timeLog;

moab::ErrorCode ReadParallelMap( const std::string& strMapFile, std::vector< MOABSInt >& vecRow,
                                 std::vector< MOABSInt >& vecCol, std::vector< MOABReal >& vecS, MOABSInt& nRows,
                                 MOABSInt& nCols, MOABSInt& nNZs );

#ifdef MOAB_HAVE_EIGEN3
void SetEigen3Matrix( Eigen::SparseMatrix< MOABReal >& mapOperator, const MOABSInt nRows, const MOABSInt nCols,
                      const std::vector< MOABSInt >& vecRow, const std::vector< MOABSInt >& vecCol,
                      const std::vector< MOABReal >& vecS );
#endif

#define PUSH_TIMER()          \
    {                         \
        start = Clock::now(); \
    }

#define POP_TIMER( EventName )                                                                                \
    {                                                                                                         \
        std::chrono::nanoseconds elapsed = duration_cast< std::chrono::nanoseconds >( Clock::now() - start ); \
        timeLog[EventName]               = elapsed;                                                           \
    }

#define PRINT_TIMER( EventName )                                                                                       \
    {                                                                                                                  \
        std::cout << "[ " << EventName                                                                                 \
                  << " ]: elapsed = " << static_cast< double >( timeLog[EventName].count() / 1e6 ) << " milli-seconds" \
                  << std::endl;                                                                                        \
    }

int main( int argc, char** argv )
{
    std::string remap_operator_filename = "";
    bool is_target_transposed     = false;
    int n_remap_iterations        = 100;
    ProgOptions opts( "Remap SpMV Mini-App" );
    opts.addRequiredArg< std::string >( "MapOperator", "Remap weights filename to optimize for SpMV", &remap_operator_filename );
    // opts.addOpt< std::string >( "srcmap,s", "Source map file name for projection", &remap_operator_filename );
    opts.addOpt< void >( "transpose,t", "Compute the tranpose operator application as well", &is_target_transposed );

    // Need option handling here for input filename
    opts.addOpt< int >( "iterations,n",
                        "Number of iterations to perform to get the average performance profile (default=100)",
                        &n_remap_iterations );

    opts.parseCommandLine( argc, argv );

    // Print problem parameter details
    std::cout << "    SpMV-Remap Application" << std::endl;
    std::cout << "-------------------------------" << std::endl;
    std::cout << "Source map             = " << remap_operator_filename << std::endl;
    std::cout << "Compute transpose map  = " << ( is_target_transposed  ? "Yes" : "No" ) << std::endl;
    std::cout << "Number of iterations   = " << n_remap_iterations << std::endl;
    std::cout << std::endl;

#ifdef MOAB_HAVE_EIGEN3
    Eigen::SparseMatrix< MOABReal > srcMapOperator;
#endif

    ErrorCode rval;
    MOABSInt nOpRows, nOpCols, nOpNNZs;

    // compute source data
    {
        std::vector< MOABSInt > opNNZRows, opNNZCols;
        std::vector< MOABReal > opNNZVals;
        PUSH_TIMER()
        rval =
            ReadParallelMap( remap_operator_filename, opNNZRows, opNNZCols, opNNZVals,
                                nOpRows, nOpCols, nOpNNZs );MB_CHK_ERR( rval );
        POP_TIMER( "ReadRemapOperator" )
        PRINT_TIMER( "ReadRemapOperator" )

#ifdef MOAB_HAVE_EIGEN3
        PUSH_TIMER()
        SetEigen3Matrix( srcMapOperator, nOpRows, nOpCols, opNNZRows, opNNZCols, opNNZVals );
        POP_TIMER( "SetRemapOperator" )
        PRINT_TIMER( "SetRemapOperator" )
#endif
    }

    std::cout << "Map operator size = [" << nOpRows << " x " << nOpCols << "] with NNZs = " << nOpNNZs << std::endl;

    // First let us perform SpMV from Source to Target
    {
        Eigen::VectorXd srcTgt = Eigen::VectorXd::Random( srcMapOperator.cols() );
        Eigen::VectorXd tgtSrc = Eigen::VectorXd::Zero( srcMapOperator.rows() );

        PUSH_TIMER()
        for( int iR = 0; iR < n_remap_iterations; ++iR )
        {
            // Project data from source to target
            tgtSrc = srcMapOperator * srcTgt;
        }
        POP_TIMER( "RemapTotalSpMV" )

        std::cout << "Average time (milli-secs) taken for " << n_remap_iterations
                  << " RemapOperator SpMV application = "
                  << static_cast< double >( timeLog["RemapTotalSpMV"].count() ) / ( 1E6 * n_remap_iterations )
                  << std::endl;
    }

    // Now let us repeat SpMV from Target to Source if requested
    if( is_target_transposed )
    {
        Eigen::VectorXd srcTgt = Eigen::VectorXd::Zero( srcMapOperator.cols() );
        Eigen::VectorXd tgtSrc = Eigen::VectorXd::Random( srcMapOperator.rows() );

        PUSH_TIMER()
        for( int iR = 0; iR < n_remap_iterations; ++iR )
        {
            // Project data from target to source through transpose application
            srcTgt = srcMapOperator.transpose() * tgtSrc;
        }
        POP_TIMER( "RemapTransposeTotalSpMV" )

        std::cout << "Average time (milli-secs) taken for " << n_remap_iterations
                  << " RemapOperator (tranpose) SpMV application = "
                  << static_cast< double >( timeLog["RemapTransposeTotalSpMV"].count() ) / ( 1E6 * n_remap_iterations )
                  << std::endl;
    }

    return 0;
}

void SetEigen3Matrix( Eigen::SparseMatrix< MOABReal >& mapOperator, const MOABSInt nRows, const MOABSInt nCols,
                      const std::vector< MOABSInt >& vecRow, const std::vector< MOABSInt >& vecCol,
                      const std::vector< MOABReal >& vecS )
{
    const size_t nS = vecS.size();
    // Let us populate the map object for every process
    mapOperator.resize( nRows, nCols );
    mapOperator.reserve( nS );
    // create a triplet vector
    typedef Eigen::Triplet< MOABReal > SparseEntry;
    std::vector< SparseEntry > tripletList( nS );

    // loop over nnz and populate the sparse matrix operator
    for( size_t innz = 0; innz < nS; ++innz )
    {
        // mapOperator.insert( vecRow[innz]-1, vecCol[innz]-1 ) = vecS[innz];
        tripletList[innz] = SparseEntry( vecRow[innz] - 1, vecCol[innz] - 1, vecS[innz] );
    }

    mapOperator.setFromTriplets( tripletList.begin(), tripletList.end() );

    mapOperator.makeCompressed();

    return;
}

moab::ErrorCode ReadParallelMap( const std::string& strMapFile, std::vector< MOABSInt >& vecRow,
                                 std::vector< MOABSInt >& vecCol, std::vector< MOABReal >& vecS, MOABSInt& nRows,
                                 MOABSInt& nCols, MOABSInt& nNZs )
{
    NcError error( NcError::silent_nonfatal );

    NcVar *varRow = NULL, *varCol = NULL, *varS = NULL;
    int nS = 0, nA = 0, nB = 0;

    // Create the NetCDF C++ interface
    NcFile ncMap( strMapFile.c_str(), NcFile::ReadOnly );

    // Read SparseMatrix entries
    NcDim* dimNA = ncMap.get_dim( "n_a" );
    if( dimNA == NULL ) {
        MB_CHK_SET_ERR( MB_FAILURE, "Map file " << strMapFile << " does not contain dimension 'nA'" );
    }

    NcDim* dimNB = ncMap.get_dim( "n_b" );
    if( dimNB == NULL ) {
        MB_CHK_SET_ERR( MB_FAILURE, "Map file " << strMapFile << " does not contain dimension 'nB'" );
    }

    NcDim* dimNS = ncMap.get_dim( "n_s" );
    if( dimNS == NULL ) {
        MB_CHK_SET_ERR( MB_FAILURE, "Map file " << strMapFile << " does not contain dimension 'nS'" );
    }

    // store total number of nonzeros
    nS = dimNS->size();
    nA = dimNA->size();
    nB = dimNB->size();

    varRow = ncMap.get_var( "row" );
    if( varRow == NULL ) {
        MB_CHK_SET_ERR( MB_FAILURE, "Map file " << strMapFile << " does not contain variable 'row'" );
    }

    varCol = ncMap.get_var( "col" );
    if( varCol == NULL ) {
        MB_CHK_SET_ERR( MB_FAILURE, "Map file " << strMapFile << " does not contain variable 'col'" );
    }

    varS = ncMap.get_var( "S" );
    if( varS == NULL ) {
        MB_CHK_SET_ERR( MB_FAILURE, "Map file " << strMapFile << " does not contain variable 'S'" );
    }

    // Resize the vectors to hold row/col indices and nnz values
    const long offsetRead = 0;
    vecRow.resize( nS );
    vecCol.resize( nS );
    vecS.resize( nS );

    varRow->set_cur( offsetRead );
    varRow->get( &( vecRow[0] ), nS );

    varCol->set_cur( offsetRead );
    varCol->get( &( vecCol[0] ), nS );

    varS->set_cur( offsetRead );
    varS->get( &( vecS[0] ), nS );

    ncMap.close();

    nRows = nB;
    nCols = nA;
    nNZs = nS;

    return moab::MB_SUCCESS;
}
