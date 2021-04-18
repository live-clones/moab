/** @example spmvApp.cpp
 * Description: read two map files to project a solution forward and in reverse direction
 * Then compute the sparse Matrix-Vector products to compute the projections from
 * source component to target component. This is strictly an example to work on single
 * node and we do not care about MPI parallelism in this experiment.
 *
 * To run: ./spmvApp srcMapFile tgtMapFile
 */

#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/ProgOptions.hpp"
#include <iostream>

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

using namespace moab;
using namespace std;

// Some global typedefs
typedef int MOABSInt;
typedef size_t MOABUInt;
typedef double MOABReal;

moab::ErrorCode ReadParallelMap( const std::string& strMapFile, std::vector< MOABSInt >& vecRow,
                                 std::vector< MOABSInt >& vecCol, std::vector< MOABReal >& vecS,
                                 Eigen::SparseMatrix< MOABReal >& mapOperator );

int main( int argc, char** argv )
{
    std::string src_map_file_name = string( MESH_DIR ) + string( "/src_map.nc" );
    std::string tgt_map_file_name = string( MESH_DIR ) + string( "/tgt_map.nc" );
    int n_remap_iterations        = 100;
    ProgOptions opts;
    opts.addOpt< std::string >( "srcmap,s", "Source map file name for projection)", &src_map_file_name );
    opts.addOpt< std::string >( "tgtmap,t", "Target map file name for projection)", &tgt_map_file_name );
    // Need option handling here for input filename
    opts.addOpt< int >( "iterations,i",
                        "Number of iterations to perform to get the average performance profile (default=100)",
                        &n_remap_iterations );

    opts.parseCommandLine( argc, argv );

    // Get MOAB instance
    Core* mb = new( std::nothrow ) Core;
    if( mb == nullptr ) return 1;

    ErrorCode rval;
    std::vector< MOABSInt > srcNNZRows, srcNNZCols;
    std::vector< MOABReal > srcNNZVals;
    Eigen::SparseMatrix< MOABReal > srcMapOperator;
    rval = ReadParallelMap( src_map_file_name, srcNNZRows, srcNNZCols, srcNNZVals, srcMapOperator );MB_CHK_ERR(rval);

    std::cout << "Source Map [" << srcMapOperator.rows() << " x " << srcMapOperator.cols() << "]: NNZs ("
              << srcMapOperator.nonZeros() << ")" << std::endl;

    std::vector< MOABSInt > tgtNNZRows, tgtNNZCols;
    std::vector< MOABReal > tgtNNZVals;
    Eigen::SparseMatrix< MOABReal > tgtMapOperator;
    rval = ReadParallelMap( tgt_map_file_name, tgtNNZRows, tgtNNZCols, tgtNNZVals, tgtMapOperator );MB_CHK_ERR( rval );

    std::cout << "Target Map [" << tgtMapOperator.rows() << " x " << tgtMapOperator.cols() << "]: NNZs ("
              << tgtMapOperator.nonZeros() << ")" << std::endl;

    delete mb;

    return 0;
}

moab::ErrorCode ReadParallelMap( const std::string& strMapFile, std::vector< MOABSInt >& vecRow,
                                 std::vector< MOABSInt >& vecCol, std::vector< MOABReal >& vecS,
                                 Eigen::SparseMatrix< MOABReal >& mapOperator )
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

    // Let us populate the map object for every process
    mapOperator.resize( nB, nA );
    mapOperator.reserve( nS );

    // loop over nnz and populate the sparse matrix operator
    for( int innz = 0; innz < nS; ++innz )
    {
        mapOperator.insert( vecRow[innz]-1, vecCol[innz]-1 ) = vecS[innz];
    }

    return moab::MB_SUCCESS;
}
