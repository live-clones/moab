/*
 * visuMapVtk.cpp
 * this tool will take a source file, target file (h5m) and a map file in nc format,
 *
 * example of usage:
 * ./mbvisumap -s source.h5m -t target.h5m -m map.nc -b startSourceID -e endSourceID  -c startTargetID -f endTargetID
 * will associate row i in map with a partial mesh
 *  will associate row j in map with a partial mesh
 *
 * can be built only if netcdf and hdf5 and eigen3 are available
 *
 *
 */
#include "moab/MOABConfig.h"

#ifndef MOAB_HAVE_EIGEN3
#error compareMaps tool requires eigen3 configuration
#endif

#include "moab/ProgOptions.hpp"
#include "moab/Core.hpp"
#include "moab/Range.hpp"

#include "netcdf.h"
#include <cmath>
#include <sstream>
#include <map>
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
            if( NC_NOERR != gdfail )                            \
            {                                                   \
                ERR_NC( gdfail )                                \
            }                                                   \
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
                if( NC_NOERR != gvfail )                               \
                {                                                      \
                    ERR_NC( gvfail )                                   \
                }                                                      \
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
            if( NC_NOERR != ivfail )                                                    \
            {                                                                           \
                ERR_NC( ivfail )                                                        \
            }                                                                           \
            ( vals ).resize( ntmp );                                                    \
            size_t ntmp1 = 0;                                                           \
            ivfail       = nc_get_vara_int( ncFile1, id, &ntmp1, &ntmp, &( vals )[0] ); \
            if( NC_NOERR != ivfail )                                                    \
            {                                                                           \
                ERR_NC( ivfail )                                                        \
            }                                                                           \
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
            if( NC_NOERR != dvfail )                                                       \
            {                                                                              \
                ERR_NC( dvfail )                                                           \
            }                                                                              \
            ( vals ).resize( ntmp );                                                       \
            size_t ntmp1 = 0;                                                              \
            dvfail       = nc_get_vara_double( ncFile1, id, &ntmp1, &ntmp, &( vals )[0] ); \
            if( NC_NOERR != dvfail )                                                       \
            {                                                                              \
                ERR_NC( dvfail )                                                           \
            }                                                                              \
        }                                                                                  \
    }

using namespace moab;
using namespace std;

int main( int argc, char* argv[] )
{

    ProgOptions opts;
    int dimSource = 2;  // for FV meshes is 2; for SE meshes, use fine mesh, dim will be 0
    int dimTarget = 2;  //
    int otype     = 0;

    std::string inputfile1, inputSource, inputTarget;
    opts.addOpt< std::string >( "map,m", "input map ", &inputfile1 );
    opts.addOpt< std::string >( "source,s", "source mesh", &inputSource );
    opts.addOpt< std::string >( "target,t", "target mesh", &inputTarget );
    opts.addOpt< int >( "dimSource,d", "dimension of source  ", &dimSource );
    opts.addOpt< int >( "dimTarget,g", "dimension of target  ", &dimTarget );
    opts.addOpt< int >( "typeOutput,o", " output type vtk(0), h5m(1)", &otype );

    int startSourceID = -1, endSourceID = -1, startTargetID = -1, endTargetID = -1;
    opts.addOpt< int >( "startSourceID,b", "start source id ", &startSourceID );
    opts.addOpt< int >( "endSourceID,e", "end source id ", &endSourceID );
    opts.addOpt< int >( "startTargetID,c", "start target id ", &startTargetID );
    opts.addOpt< int >( "endTargetID,f", "end target id ", &endTargetID );
    //  -b startSourceID -e endSourceID  -c startTargetID -f endTargetID

    opts.parseCommandLine( argc, argv );

    std::string extension = ".vtk";
    if( 1 == otype ) extension = ".h5m";
    // Open netcdf/exodus file
    int fail = nc_open( inputfile1.c_str(), 0, &ncFile1 );
    if( NC_NOWRITE != fail )
    {
        ERR_NC( fail )
    }

    std::cout << " opened " << inputfile1 << " for map 1 \n";

    int temp_dim;
    int na1, nb1, ns1;
    GET_DIM1( temp_dim, "n_a", na1 );
    GET_DIM1( temp_dim, "n_b", nb1 );
    GET_DIM1( temp_dim, "n_s", ns1 );
    std::cout << " n_a, n_b, n_s : " << na1 << ", " << nb1 << ", " << ns1 << " for map 1 \n";
    std::vector< int > col1( ns1 ), row1( ns1 );
    std::vector< double > val1( ns1 );
    int idrow1, idcol1, ids1;
    GET_1D_INT_VAR1( "row", idrow1, row1 );
    GET_1D_INT_VAR1( "col", idcol1, col1 );
    GET_1D_DBL_VAR1( "S", ids1, val1 );
    // first matrix
    typedef Eigen::Triplet< double > Triplet;
    std::vector< Triplet > tripletList;
    tripletList.reserve( ns1 );
    for( int iv = 0; iv < ns1; iv++ )
    {
        // all row and column indices are 1-based in the map file.
        // inside Eigen3, we will use 0-based; then we will have to add back 1 when
        //     we dump out the files
        tripletList.push_back( Triplet( row1[iv] - 1, col1[iv] - 1, val1[iv] ) );
    }
    Eigen::SparseMatrix< double > weight1( nb1, na1 );

    weight1.setFromTriplets( tripletList.begin(), tripletList.end() );
    weight1.makeCompressed();
    // we read the matrix; now read moab source and target
    Core core;
    Interface* mb = &core;
    ErrorCode rval;
    Tag gtag = mb->globalId_tag();
    EntityHandle sourceSet, targetSet;
    // those are the maps from global ids to the moab entity handles corresponding to those global ids
    //   which are corresponding to the global DOFs
    map< int, EntityHandle > sourceHandles;
    map< int, EntityHandle > targetHandles;
    rval = mb->create_meshset( MESHSET_SET, sourceSet );MB_CHK_SET_ERR( rval, "can't create source mesh set" );
    rval = mb->create_meshset( MESHSET_SET, targetSet );MB_CHK_SET_ERR( rval, "can't create target mesh set" );
    const char* readopts = "";
    rval                 = mb->load_file( inputSource.c_str(), &sourceSet, readopts );MB_CHK_SET_ERR( rval, "Failed to read" );
    rval = mb->load_file( inputTarget.c_str(), &targetSet, readopts );MB_CHK_SET_ERR( rval, "Failed to read" );
    Range sRange;
    rval = mb->get_entities_by_dimension( sourceSet, dimSource, sRange );MB_CHK_SET_ERR( rval, "Failed to get sRange" );
    vector< int > sids;
    sids.resize( sRange.size() );
    rval = mb->tag_get_data( gtag, sRange, &sids[0] );MB_CHK_SET_ERR( rval, "Failed to get ids for srange" );
    // all global ids are 1 based in the source file, and they correspond to the dofs in the map file
    for( size_t i = 0; i < sids.size(); i++ )
    {
        EntityHandle eh    = sRange[i];
        int gid            = sids[i];
        sourceHandles[gid] = eh;
    }
    Range tRange;
    rval = mb->get_entities_by_dimension( targetSet, dimTarget, tRange );MB_CHK_SET_ERR( rval, "Failed to get tRange" );
    vector< int > tids;
    tids.resize( tRange.size() );
    rval = mb->tag_get_data( gtag, tRange, &tids[0] );MB_CHK_SET_ERR( rval, "Failed to get ids for trange" );
    // all global ids are 1 based in the target file, and they correspond to the dofs in the map file
    for( size_t i = 0; i < tids.size(); i++ )
    {
        EntityHandle eh    = tRange[i];
        int gid            = tids[i];
        targetHandles[gid] = eh;
    }
    // a dense tag for weights
    Tag wtag;
    double defVal = 0;

    std::string name_map = inputfile1;
    // strip last 3 chars (.nc extension)
    name_map.erase( name_map.begin() + name_map.length() - 3, name_map.end() );
    // if path , remove from name
    size_t pos = name_map.rfind( '/', name_map.length() );
    if( pos != std::string::npos ) name_map = name_map.erase( 0, pos + 1 );

    rval = mb->tag_get_handle( "weight", 1, MB_TYPE_DOUBLE, wtag, MB_TAG_CREAT | MB_TAG_DENSE, &defVal );MB_CHK_SET_ERR( rval, "Failed to create weight" );
    EntityHandle partialSet;
    rval = mb->create_meshset( MESHSET_SET, partialSet );MB_CHK_SET_ERR( rval, "can't create partial set" );
    // how to get a complete row in sparse matrix? Or a complete column ?
    for( int col = startSourceID - 1; col <= endSourceID - 1; col++ )
    {
        Range targetEnts;  // will find entries for column col-1 in sparse matrix weight1, and its entries
        // will assign a dense tag with values, and write out the file
        if( col < 0 ) continue;
        Eigen::SparseVector< double > colVect = weight1.col( col );
        // the row indices correspond to target cells
        for( Eigen::SparseVector< double >::InnerIterator it( colVect ); it; ++it )
        {
            double weight = it.value();  // == vec[ it.index() ]
                                         // we add back the 1 that we subtract
            int globalIdRow = it.index() + 1;
            EntityHandle th = targetHandles[globalIdRow];
            targetEnts.insert( th );
            rval = mb->tag_set_data( wtag, &th, 1, &weight );MB_CHK_SET_ERR( rval, "Failed to set weight tag on target" );
        }

        if( dimTarget == 0 )
        {
            Range adjCells;
            rval = mb->get_adjacencies( targetEnts, 2, false, adjCells, Interface::UNION );MB_CHK_SET_ERR( rval, " can't get adj cells " );
            targetEnts.merge( adjCells );
        }

        rval = mb->add_entities( partialSet, targetEnts );MB_CHK_SET_ERR( rval, "Failed to add target entities to partial set" );
        // write now the set in a numbered file
        std::stringstream fff;
        fff << name_map << "_column" << col + 1 << extension;
        rval = mb->write_mesh( fff.str().c_str(), &partialSet, 1 );MB_CHK_ERR( rval );
        // remove from partial set the entities it has
        rval = mb->clear_meshset( &partialSet, 1 );MB_CHK_SET_ERR( rval, "Failed to empty partial set" );
    }

    // how to get a complete row in sparse matrix?
    for( int row = startTargetID - 1; row <= endTargetID - 1; row++ )
    {
        Range sourceEnts;  // will find entries for row in sparse matrix weight1, and its entries
        // will assign a dense tag with values, and write out the file
        if( row < 0 ) continue;
        Eigen::SparseVector< double > rowVect = weight1.row( row );
        // the row indices correspond to target cells
        for( Eigen::SparseVector< double >::InnerIterator it( rowVect ); it; ++it )
        {
            double weight   = it.value();  // == vec[ it.index() ]
            int globalIdCol = it.index() + 1;
            EntityHandle sh = sourceHandles[globalIdCol];
            sourceEnts.insert( sh );
            rval = mb->tag_set_data( wtag, &sh, 1, &weight );MB_CHK_SET_ERR( rval, "Failed to set weight tag on source" );
        }
        if( dimSource == 0 )
        {
            Range adjCells;
            rval = mb->get_adjacencies( sourceEnts, 2, false, adjCells, Interface::UNION );MB_CHK_SET_ERR( rval, " can't get adj cells " );
            sourceEnts.merge( adjCells );
        }
        rval = mb->add_entities( partialSet, sourceEnts ); MB_CHK_SET_ERR( rval, "Failed to add source entities" );
        // write now the set in a numbered file
        std::stringstream fff;
        fff << name_map << "_row" << row + 1 << extension;
        rval = mb->write_mesh( fff.str().c_str(), &partialSet, 1 );MB_CHK_ERR( rval );
        rval = mb->clear_meshset( &partialSet, 1 );MB_CHK_SET_ERR( rval, "Failed to empty partial set" );
    }
    return 0;
}
