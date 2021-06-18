/*
 * compareMaps.cpp
 * this tool will take 2 existing map files in nc format, and compare their sparse matrices
 * we can compare xc, yc, areas, fractions with ncdiff from nco
 * maybe there is another utility in nco, need to ask Charlie Zender
 *
 * example of usage:
 * ./mbcmpmaps map1.nc map2.nc
 * will look for row, col, S entries, and use eigen3 sparse matrix constructor
 *
 * can be built only if netcdf and eigen3 are available
 *
 *
 */
#include "moab/MOABConfig.h"

#ifndef MOAB_HAVE_EIGEN3
#error compareMaps tool requires eigen3 configuration
#endif

#include "moab/ProgOptions.hpp"

#include "netcdf.h"
#include <cmath>
#include <sstream>
#include <Eigen/Sparse>
#define ERR_NC( e )                                \
    {                                              \
        printf( "Error: %s\n", nc_strerror( e ) ); \
        exit( 2 );                                 \
    }

// copy from ReadNCDF.cpp some useful macros for reading from a netcdf file (exodus?)
// ncFile is an integer initialized when opening the nc file in read mode

int ncFile1;

#define INS_ID( stringvar, prefix, id ) sprintf( stringvar, prefix, id )

#define GET_DIM1( ncdim, name, val )                            \
    {                                                           \
        int gdfail = nc_inq_dimid( ncFile1, name, &( ncdim ) ); \
        if( NC_NOERR == gdfail )                                \
        {                                                       \
            size_t tmp_val;                                     \
            gdfail = nc_inq_dimlen( ncFile1, ncdim, &tmp_val ); \
            if( NC_NOERR != gdfail ) { ERR_NC( gdfail ) }       \
            else                                                \
                ( val ) = tmp_val;                              \
        }                                                       \
        else                                                    \
            ( val ) = 0;                                        \
    }

#define GET_DIMB1( ncdim, name, varname, id, val ) \
    INS_ID( name, varname, id );                   \
    GET_DIM1( ncdim, name, val );

#define GET_VAR1( name, id, dims )                                     \
    {                                                                  \
        ( id )     = -1;                                               \
        int gvfail = nc_inq_varid( ncFile1, name, &( id ) );           \
        if( NC_NOERR == gvfail )                                       \
        {                                                              \
            int ndims;                                                 \
            gvfail = nc_inq_varndims( ncFile1, id, &ndims );           \
            if( NC_NOERR == gvfail )                                   \
            {                                                          \
                ( dims ).resize( ndims );                              \
                gvfail = nc_inq_vardimid( ncFile1, id, &( dims )[0] ); \
                if( NC_NOERR != gvfail ) { ERR_NC( gvfail ) }          \
            }                                                          \
        }                                                              \
    }

#define GET_1D_INT_VAR1( name, id, vals )                                               \
    {                                                                                   \
        GET_VAR1( name, id, vals );                                                     \
        if( -1 != ( id ) )                                                              \
        {                                                                               \
            size_t ntmp;                                                                \
            int ivfail = nc_inq_dimlen( ncFile1, ( vals )[0], &ntmp );                  \
            if( NC_NOERR != ivfail ) { ERR_NC( ivfail ) }                               \
            ( vals ).resize( ntmp );                                                    \
            size_t ntmp1 = 0;                                                           \
            ivfail       = nc_get_vara_int( ncFile1, id, &ntmp1, &ntmp, &( vals )[0] ); \
            if( NC_NOERR != ivfail ) { ERR_NC( ivfail ) }                               \
        }                                                                               \
    }

#define GET_1D_DBL_VAR1( name, id, vals )                                                  \
    {                                                                                      \
        std::vector< int > dum_dims;                                                       \
        GET_VAR1( name, id, dum_dims );                                                    \
        if( -1 != ( id ) )                                                                 \
        {                                                                                  \
            size_t ntmp;                                                                   \
            int dvfail = nc_inq_dimlen( ncFile1, dum_dims[0], &ntmp );                     \
            if( NC_NOERR != dvfail ) { ERR_NC( dvfail ) }                                  \
            ( vals ).resize( ntmp );                                                       \
            size_t ntmp1 = 0;                                                              \
            dvfail       = nc_get_vara_double( ncFile1, id, &ntmp1, &ntmp, &( vals )[0] ); \
            if( NC_NOERR != dvfail ) { ERR_NC( dvfail ) }                                  \
        }                                                                                  \
    }

int ncFile2;

#define GET_DIM2( ncdim, name, val )                            \
    {                                                           \
        int gdfail = nc_inq_dimid( ncFile2, name, &( ncdim ) ); \
        if( NC_NOERR == gdfail )                                \
        {                                                       \
            size_t tmp_val;                                     \
            gdfail = nc_inq_dimlen( ncFile2, ncdim, &tmp_val ); \
            if( NC_NOERR != gdfail ) { ERR_NC( gdfail ) }       \
            else                                                \
                ( val ) = tmp_val;                              \
        }                                                       \
        else                                                    \
            ( val ) = 0;                                        \
    }

#define GET_DIMB2( ncdim, name, varname, id, val ) \
    INS_ID( name, varname, id );                   \
    GET_DIM2( ncdim, name, val );

#define GET_VAR2( name, id, dims )                                     \
    {                                                                  \
        ( id )     = -1;                                               \
        int gvfail = nc_inq_varid( ncFile2, name, &( id ) );           \
        if( NC_NOERR == gvfail )                                       \
        {                                                              \
            int ndims;                                                 \
            gvfail = nc_inq_varndims( ncFile2, id, &ndims );           \
            if( NC_NOERR == gvfail )                                   \
            {                                                          \
                ( dims ).resize( ndims );                              \
                gvfail = nc_inq_vardimid( ncFile2, id, &( dims )[0] ); \
                if( NC_NOERR != gvfail ) { ERR_NC( gvfail ) }          \
            }                                                          \
        }                                                              \
    }

#define GET_1D_INT_VAR2( name, id, vals )                                               \
    {                                                                                   \
        GET_VAR2( name, id, vals );                                                     \
        if( -1 != ( id ) )                                                              \
        {                                                                               \
            size_t ntmp;                                                                \
            int ivfail = nc_inq_dimlen( ncFile2, ( vals )[0], &ntmp );                  \
            if( NC_NOERR != ivfail ) { ERR_NC( ivfail ) }                               \
            ( vals ).resize( ntmp );                                                    \
            size_t ntmp1 = 0;                                                           \
            ivfail       = nc_get_vara_int( ncFile2, id, &ntmp1, &ntmp, &( vals )[0] ); \
            if( NC_NOERR != ivfail ) { ERR_NC( ivfail ) }                               \
        }                                                                               \
    }

#define GET_1D_DBL_VAR2( name, id, vals )                                                  \
    {                                                                                      \
        std::vector< int > dum_dims;                                                       \
        GET_VAR2( name, id, dum_dims );                                                    \
        if( -1 != ( id ) )                                                                 \
        {                                                                                  \
            size_t ntmp;                                                                   \
            int dvfail = nc_inq_dimlen( ncFile2, dum_dims[0], &ntmp );                     \
            if( NC_NOERR != dvfail ) { ERR_NC( dvfail ) }                                  \
            ( vals ).resize( ntmp );                                                       \
            size_t ntmp1 = 0;                                                              \
            dvfail       = nc_get_vara_double( ncFile2, id, &ntmp1, &ntmp, &( vals )[0] ); \
            if( NC_NOERR != dvfail ) { ERR_NC( dvfail ) }                                  \
        }                                                                                  \
    }

#define GET_2D_DBL_VAR1( name, id, vals )                                                   \
    {                                                                                       \
        std::vector< int > dum_dims;                                                        \
        GET_VAR1( name, id, dum_dims );                                                     \
        if( -1 != ( id ) )                                                                  \
        {                                                                                   \
            size_t ntmp[2];                                                                 \
            int dvfail = nc_inq_dimlen( ncFile1, dum_dims[0], &ntmp[0] );                   \
            if( NC_NOERR != dvfail ) { ERR_NC( dvfail ) }                                   \
            dvfail = nc_inq_dimlen( ncFile1, dum_dims[1], &ntmp[1] );                       \
            if( NC_NOERR != dvfail ) { ERR_NC( dvfail ) }                                   \
            ( vals ).resize( ntmp[0] * ntmp[1] );                                           \
            size_t ntmp1[2] = { 0, 0 };                                                     \
            dvfail          = nc_get_vara_double( ncFile1, id, ntmp1, ntmp, &( vals )[0] ); \
            if( NC_NOERR != dvfail ) { ERR_NC( dvfail ) }                                   \
        }                                                                                   \
    }

#define GET_2D_DBL_VAR2( name, id, vals )                                                   \
    {                                                                                       \
        std::vector< int > dum_dims;                                                        \
        GET_VAR2( name, id, dum_dims );                                                     \
        if( -1 != ( id ) )                                                                  \
        {                                                                                   \
            size_t ntmp[2];                                                                 \
            int dvfail = nc_inq_dimlen( ncFile2, dum_dims[0], &ntmp[0] );                   \
            if( NC_NOERR != dvfail ) { ERR_NC( dvfail ) }                                   \
            dvfail = nc_inq_dimlen( ncFile2, dum_dims[1], &ntmp[1] );                       \
            if( NC_NOERR != dvfail ) { ERR_NC( dvfail ) }                                   \
            ( vals ).resize( ntmp[0] * ntmp[1] );                                           \
            size_t ntmp1[2] = { 0, 0 };                                                     \
            dvfail          = nc_get_vara_double( ncFile2, id, ntmp1, ntmp, &( vals )[0] ); \
            if( NC_NOERR != dvfail ) { ERR_NC( dvfail ) }                                   \
        }                                                                                   \
    }

typedef Eigen::Map< Eigen::VectorXd > EigenV;

void diff_vect( const char* var_name, int n )
{
    // compare frac_a between maps
    std::vector< double > fraca1( n ), fraca2( n );
    int idfa1, idfa2;
    GET_1D_DBL_VAR1( var_name, idfa1, fraca1 );
    EigenV fa1( fraca1.data(), n );
    GET_1D_DBL_VAR2( var_name, idfa2, fraca2 );
    EigenV fa2( fraca2.data(), n );

    EigenV diff( fraca2.data(), n );
    diff = fa1 - fa2;

    int imin, imax;
    double minV = diff.minCoeff( &imin );
    double maxV = diff.maxCoeff( &imax );
    std::cout << var_name << " diff norm: " << diff.norm() << " min at " << imin << " : " << minV << " max at " << imax
              << " : " << maxV << "\n";
    return;
}

void diff_2d_vect( const char* var_name, int n )
{
    // compare frac_a between maps
    std::vector< double > fraca1( n ), fraca2( n );
    int idfa1, idfa2;
    GET_2D_DBL_VAR1( var_name, idfa1, fraca1 );
    EigenV fa1( fraca1.data(), n );
    GET_2D_DBL_VAR2( var_name, idfa2, fraca2 );
    EigenV fa2( fraca2.data(), n );
    std::cout << var_name << " diff norm: " << ( fa1 - fa2 ).norm() << "\n";
    return;
}
int main( int argc, char* argv[] )
{

    ProgOptions opts;

    std::string inputfile1, inputfile2;
    opts.addOpt< std::string >( "firstMap,i", "input filename 1", &inputfile1 );
    opts.addOpt< std::string >( "secondMap,j", "input second map", &inputfile2 );

    opts.parseCommandLine( argc, argv );

    // Open netcdf/exodus file
    int fail = nc_open( inputfile1.c_str(), 0, &ncFile1 );
    if( NC_NOWRITE != fail ) { ERR_NC( fail ) }

    std::cout << " opened " << inputfile1 << " for map 1 \n";

    fail = nc_open( inputfile2.c_str(), 0, &ncFile2 );
    if( NC_NOWRITE != fail ) { ERR_NC( fail ) }

    std::cout << " opened " << inputfile2 << " for map 2 \n";
    int temp_dim;
    int na1, nb1, ns1, na2, nb2, ns2, nv_a, nv_b, nv_a2, nv_b2;
    GET_DIM1( temp_dim, "n_a", na1 );
    GET_DIM2( temp_dim, "n_a", na2 );
    GET_DIM1( temp_dim, "n_b", nb1 );
    GET_DIM2( temp_dim, "n_b", nb2 );
    GET_DIM1( temp_dim, "n_s", ns1 );
    GET_DIM2( temp_dim, "n_s", ns2 );
    GET_DIM1( temp_dim, "nv_a", nv_a );
    GET_DIM2( temp_dim, "nv_a", nv_a2 );
    GET_DIM1( temp_dim, "nv_b", nv_b );
    GET_DIM2( temp_dim, "nv_b", nv_b2 );
    if( nv_a != nv_a2 || nv_b != nv_b2 )
    {
        std::cout << " different nv dimensions:" << nv_a << " == " << nv_a2 << "  or " << nv_b << " == " << nv_b2
                  << "  bail out \n";
        return 1;
    }
    std::cout << " n_a, n_b, n_s : " << na1 << ", " << nb1 << ", " << ns1 << " for map 1 \n";

    std::cout << " n_a, n_b, n_s : " << na2 << ", " << nb2 << ", " << ns2 << " for map 2 \n";
    if( na1 != na2 || nb1 != nb2 )
    {
        std::cout << " different dimensions bail out \n";
        return 1;
    }
    std::vector< int > col1( ns1 ), row1( ns1 );
    std::vector< int > col2( ns2 ), row2( ns2 );

    std::vector< double > val1( ns1 ), val2( ns2 );

    int idrow1, idcol1, idrow2, idcol2, ids1, ids2;
    GET_1D_INT_VAR1( "row", idrow1, row1 );
    GET_1D_INT_VAR1( "col", idcol1, col1 );
    GET_1D_DBL_VAR1( "S", ids1, val1 );
    GET_1D_INT_VAR2( "row", idrow2, row2 );
    GET_1D_INT_VAR2( "col", idcol2, col2 );
    GET_1D_DBL_VAR2( "S", ids2, val2 );

    // first matrix
    typedef Eigen::Triplet< double > Triplet;
    std::vector< Triplet > tripletList;
    tripletList.reserve( ns1 );
    for( int iv = 0; iv < ns1; iv++ )
    {
        tripletList.push_back( Triplet( row1[iv] - 1, col1[iv] - 1, val1[iv] ) );
    }
    Eigen::SparseMatrix< double > weight1( nb1, na1 );

    weight1.setFromTriplets( tripletList.begin(), tripletList.end() );
    weight1.makeCompressed();

    if( ns1 != ns2 ) tripletList.resize( ns2 );
    for( int iv = 0; iv < ns2; iv++ )
    {
        tripletList[iv] = Triplet( row2[iv] - 1, col2[iv] - 1, val2[iv] );
    }
    Eigen::SparseMatrix< double > weight2( nb2, na2 );

    weight2.setFromTriplets( tripletList.begin(), tripletList.end() );
    weight2.makeCompressed();

    Eigen::SparseMatrix< double > diff = weight1 - weight2;
    double maxv                        = diff.coeffs().maxCoeff();
    double minv                        = diff.coeffs().minCoeff();
    std::cout << " euclidian norm for difference: " << diff.norm()
              << " \n squared norm for difference: " << diff.squaredNorm() << "\n"
              << " minv: " << minv << " maxv: " << maxv << "\n";

    // compare frac_a between maps
    diff_vect( "frac_a", na1 );
    diff_vect( "frac_b", nb1 );
    diff_vect( "area_a", na1 );
    diff_vect( "area_b", nb1 );
    diff_vect( "yc_a", na1 );
    diff_vect( "yc_b", nb1 );
    diff_vect( "xc_a", na1 );
    diff_vect( "xc_b", nb1 );
    diff_2d_vect( "xv_a", na1 * nv_a );
    diff_2d_vect( "yv_a", na1 * nv_a );
    diff_2d_vect( "xv_b", nb1 * nv_b );
    diff_2d_vect( "yv_b", nb1 * nv_b );

    return 0;
}
