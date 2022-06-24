/*
 * NCHelperScrip.cpp
 */

#include "NCHelperScrip.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/IntxMesh/IntxUtils.hpp"
#ifdef MOAB_HAVE_MPI
#include "moab/ParallelMergeMesh.hpp"
#endif
#ifdef MOAB_HAVE_ZOLTAN
#include "moab/ZoltanPartitioner.hpp"
#endif

namespace moab
{

bool NCHelperScrip::can_read_file( ReadNC* readNC, int /*fileId*/ )
{
    std::vector< std::string >& dimNames = readNC->dimNames;

    // If dimension names "grid_size" AND "grid_corners" AND "grid_rank" exist then it should be the Scrip grid
    if( ( std::find( dimNames.begin(), dimNames.end(), std::string( "grid_size" ) ) != dimNames.end() ) &&
        ( std::find( dimNames.begin(), dimNames.end(), std::string( "grid_corners" ) ) != dimNames.end() ) &&
        ( std::find( dimNames.begin(), dimNames.end(), std::string( "grid_rank" ) ) != dimNames.end() ) )
    {

        return true;
    }

    return false;
}
ErrorCode NCHelperScrip::init_mesh_vals()
{
    Interface*& mbImpl                   = _readNC->mbImpl;
    std::vector< std::string >& dimNames = _readNC->dimNames;
    std::vector< int >& dimLens          = _readNC->dimLens;

    unsigned int idx;
    std::vector< std::string >::iterator vit;

    // get grid_size
    if( ( vit = std::find( dimNames.begin(), dimNames.end(), "grid_size" ) ) != dimNames.end() )
    {
        idx       = vit - dimNames.begin();
        grid_size = dimLens[idx];
    }

    // get grid_corners
    if( ( vit = std::find( dimNames.begin(), dimNames.end(), "grid_corners" ) ) != dimNames.end() )
    {
        idx          = vit - dimNames.begin();
        grid_corners = dimLens[idx];
    }
    // do not need conventional tags
    Tag convTagsCreated = 0;
    int def_val         = 0;
    ErrorCode rval      = mbImpl->tag_get_handle( "__CONV_TAGS_CREATED", 1, MB_TYPE_INTEGER, convTagsCreated,
                                             MB_TAG_SPARSE | MB_TAG_CREAT, &def_val );MB_CHK_SET_ERR( rval, "Trouble getting _CONV_TAGS_CREATED tag" );
    int create_conv_tags_flag = 1;
    rval                      = mbImpl->tag_set_data( convTagsCreated, &_fileSet, 1, &create_conv_tags_flag );MB_CHK_SET_ERR( rval, "Trouble setting _CONV_TAGS_CREATED tag" );

    // decide now the units, by looking at grid_center_lon
    int xCellVarId;
    int success = NCFUNC( inq_varid )( _fileId, "grid_center_lon", &xCellVarId );
    if( success ) MB_CHK_SET_ERR( MB_FAILURE, "Trouble getting grid_center_lon" );
    std::map< std::string, ReadNC::VarData >& varInfo = _readNC->varInfo;
    auto vmit                                         = varInfo.find( "grid_center_lon" );
    if( varInfo.end() == vmit )
        MB_SET_ERR( MB_FAILURE, "Couldn't find variable "
                                    << "grid_center_lon" );
    ReadNC::VarData& glData = vmit->second;
    auto attIt              = glData.varAtts.find( "units" );
    if( attIt != glData.varAtts.end() )
    {
        unsigned int sz = attIt->second.attLen;
        std::string att_data;
        att_data.resize( sz + 1 );
        att_data[sz] = '\000';
        success =
            NCFUNC( get_att_text )( _fileId, attIt->second.attVarId, attIt->second.attName.c_str(), &att_data[0] );
        if( 0 == success && att_data.find( "radians" ) != std::string::npos ) degrees = false;
    }

    return MB_SUCCESS;
}
ErrorCode NCHelperScrip::create_mesh( Range& faces )
{
    Interface*& mbImpl  = _readNC->mbImpl;
    DebugOutput& dbgOut = _readNC->dbgOut;
    Tag& mGlobalIdTag   = _readNC->mGlobalIdTag;
    ErrorCode rval;

#ifdef MOAB_HAVE_MPI
    int rank         = 0;
    int procs        = 1;
    bool& isParallel = _readNC->isParallel;
    if( isParallel )
    {
        ParallelComm*& myPcomm = _readNC->myPcomm;
        rank                   = myPcomm->proc_config().proc_rank();
        procs                  = myPcomm->proc_config().proc_size();
    }

    if( procs >= 2 )
    {
        // Shift rank to obtain a rotated trivial partition
        int shifted_rank           = rank;
        int& trivialPartitionShift = _readNC->trivialPartitionShift;
        if( trivialPartitionShift > 0 ) shifted_rank = ( rank + trivialPartitionShift ) % procs;

        // Compute the number of local cells on this proc
        nLocalCells = int( std::floor( 1.0 * grid_size / procs ) );

        // The starting global cell index in the MPAS file for this proc
        int start_cell_idx = shifted_rank * nLocalCells;

        // Number of extra cells after equal split over procs
        int iextra = grid_size % procs;

        // Allocate extra cells over procs
        if( shifted_rank < iextra ) nLocalCells++;
        start_cell_idx += std::min( shifted_rank, iextra );

        start_cell_idx++;  // 0 based -> 1 based

        // Redistribute local cells after trivial partition (e.g. apply Zoltan partition)
        ErrorCode rval = redistribute_local_cells( start_cell_idx );MB_CHK_SET_ERR( rval, "Failed to redistribute local cells after trivial partition" );
    }
    else
    {
        nLocalCells = grid_size;
        localGidCells.insert( 1, nLocalCells );
    }
#else
    nLocalCells = grid_size;
    localGidCells.insert( 1, nLocalCells );
#endif
    dbgOut.tprintf( 1, " localGidCells.psize() = %d\n", (int)localGidCells.psize() );
    dbgOut.tprintf( 1, " localGidCells.size() = %d\n", (int)localGidCells.size() );

    // double grid_corner_lat(grid_size, grid_corners) ;
    // double grid_corner_lon(grid_size, grid_corners) ;
    int xvId, yvId;
    int success = NCFUNC( inq_varid )( _fileId, "grid_corner_lon", &xvId );
    if( success ) MB_SET_ERR( MB_FAILURE, "Failed to get variable id of grid_corner_lon" );
    success = NCFUNC( inq_varid )( _fileId, "grid_corner_lat", &yvId );
    if( success ) MB_SET_ERR( MB_FAILURE, "Failed to get variable id of grid_corner_lat" );

    // important upgrade: read masks if they exist, and save them as tags
    int gmId           = -1;
    int sizeMasks      = 0;
#ifdef MOAB_HAVE_PNETCDF
    int factorRequests = 2;  // we would read in general only 2 variables, xv and yv
#endif
    success     = NCFUNC( inq_varid )( _fileId, "grid_imask", &gmId );
    Tag maskTag = 0;  // not sure yet if we have the masks or not
    if( success )
    {
        gmId = -1;  // we do not have masks
    }
    else
    {
        sizeMasks      = nLocalCells;
#ifdef MOAB_HAVE_PNETCDF
        factorRequests = 3;  // we also need to read masks distributed
#endif
        // create the maskTag GRID_IMASK, with default value of 1
        int def_val = 1;
        rval =
            mbImpl->tag_get_handle( "GRID_IMASK", 1, MB_TYPE_INTEGER, maskTag, MB_TAG_DENSE | MB_TAG_CREAT, &def_val );MB_CHK_SET_ERR( rval, "Trouble creating GRID_IMASK tag" );
    }

    std::vector< double > xv( nLocalCells * grid_corners );
    std::vector< double > yv( nLocalCells * grid_corners );
    std::vector< int > masks( sizeMasks );
#ifdef MOAB_HAVE_PNETCDF
    size_t nb_reads = localGidCells.psize();
    std::vector< int > requests( nb_reads * factorRequests );
    std::vector< int > statuss( nb_reads * factorRequests );
    size_t idxReq = 0;
#endif
    size_t indexInArray     = 0;
    size_t indexInMaskArray = 0;
    for( Range::pair_iterator pair_iter = localGidCells.pair_begin(); pair_iter != localGidCells.pair_end();
         ++pair_iter )
    {
        EntityHandle starth      = pair_iter->first;
        EntityHandle endh        = pair_iter->second;
        NCDF_SIZE read_starts[2] = { static_cast< NCDF_SIZE >( starth - 1 ), 0 };
        NCDF_SIZE read_counts[2] = { static_cast< NCDF_SIZE >( endh - starth + 1 ),
                                     static_cast< NCDF_SIZE >( grid_corners ) };

        // Do a partial read in each subrange
#ifdef MOAB_HAVE_PNETCDF
        success = NCFUNCREQG( _vara_double )( _fileId, xvId, read_starts, read_counts, &( xv[indexInArray] ),
                                              &requests[idxReq++] );
#else
        success = NCFUNCAG( _vara_double )( _fileId, xvId, read_starts, read_counts, &( xv[indexInArray] ) );
#endif
        if( success ) MB_SET_ERR( MB_FAILURE, "Failed to read grid_corner_lon data in a loop" );

            // Do a partial read in each subrange
#ifdef MOAB_HAVE_PNETCDF
        success = NCFUNCREQG( _vara_double )( _fileId, yvId, read_starts, read_counts, &( yv[indexInArray] ),
                                              &requests[idxReq++] );
#else
        success = NCFUNCAG( _vara_double )( _fileId, yvId, read_starts, read_counts, &( yv[indexInArray] ) );
#endif
        if( success ) MB_SET_ERR( MB_FAILURE, "Failed to read grid_corner_lat data in a loop" );
        // Increment the index for next subrange
        indexInArray += ( endh - starth + 1 ) * grid_corners;

        if( gmId >= 0 )  // it means we need to read masks too, distributed:
        {
            NCDF_SIZE read_st = static_cast< NCDF_SIZE >( starth - 1 );
            NCDF_SIZE read_ct = static_cast< NCDF_SIZE >( endh - starth + 1 );
            // Do a partial read in each subrange, for mask variable:
#ifdef MOAB_HAVE_PNETCDF
            success = NCFUNCREQG( _vara_int )( _fileId, gmId, &read_st, &read_ct, &( masks[indexInMaskArray] ),
                                               &requests[idxReq++] );
#else
            success = NCFUNCAG( _vara_int )( _fileId, gmId, &read_st, &read_ct, &( masks[indexInMaskArray] ) );
#endif
            if( success ) MB_SET_ERR( MB_FAILURE, "Failed on mask read " );
            indexInMaskArray += endh - starth + 1;
        }
    }

#ifdef MOAB_HAVE_PNETCDF
    // Wait outside the loop
    success = NCFUNC( wait_all )( _fileId, requests.size(), &requests[0], &statuss[0] );
    if( success ) MB_SET_ERR( MB_FAILURE, "Failed on wait_all" );
#endif

    // so we read xv, yv for all corners in the local mesh, and masks if they exist

    // Create vertices; first identify different ones, with a tolerance
    std::map< Node3D, EntityHandle > vertex_map;

    // Set vertex coordinates
    // will read all xv, yv, but use only those with correct mask on

    int elem_index = 0;   // local index in netcdf arrays
    double pideg   = 1.;  // radians
    if( degrees ) pideg = acos( -1.0 ) / 180.0;

    for( ; elem_index < nLocalCells; elem_index++ )
    {
        // set area and fraction on those elements too
        for( int k = 0; k < grid_corners; k++ )
        {
            int index_v_arr = grid_corners * elem_index + k;
            double x, y;
            x             = xv[index_v_arr];
            y             = yv[index_v_arr];
            double cosphi = cos( pideg * y );
            double zmult  = sin( pideg * y );
            double xmult  = cosphi * cos( x * pideg );
            double ymult  = cosphi * sin( x * pideg );
            Node3D pt( xmult, ymult, zmult );
            vertex_map[pt] = 0;
        }
    }
    int nLocalVertices = (int)vertex_map.size();
    std::vector< double* > arrays;
    EntityHandle start_vertex, vtx_handle;
    rval = _readNC->readMeshIface->get_node_coords( 3, nLocalVertices, 0, start_vertex, arrays );MB_CHK_SET_ERR( rval, "Failed to create local vertices" );

    vtx_handle = start_vertex;
    // Copy vertex coordinates into entity sequence coordinate arrays
    // and copy handle into vertex_map.
    double *x = arrays[0], *y = arrays[1], *z = arrays[2];
    for( auto i = vertex_map.begin(); i != vertex_map.end(); ++i )
    {
        i->second = vtx_handle;
        ++vtx_handle;
        *x = i->first.coords[0];
        ++x;
        *y = i->first.coords[1];
        ++y;
        *z = i->first.coords[2];
        ++z;
    }

    EntityHandle start_cell;
    int nv              = grid_corners;
    EntityType mdb_type = MBVERTEX;
    if( nv == 3 )
        mdb_type = MBTRI;
    else if( nv == 4 )
        mdb_type = MBQUAD;
    else if( nv > 4 )  // (nv > 4)
        mdb_type = MBPOLYGON;

    Range tmp_range;
    EntityHandle* conn_arr;

    rval = _readNC->readMeshIface->get_element_connect( nLocalCells, nv, mdb_type, 0, start_cell, conn_arr );MB_CHK_SET_ERR( rval, "Failed to create local cells" );
    tmp_range.insert( start_cell, start_cell + nLocalCells - 1 );

    elem_index = 0;

    for( ; elem_index < nLocalCells; elem_index++ )
    {
        for( int k = 0; k < nv; k++ )
        {
            int index_v_arr = nv * elem_index + k;
            if( nv > 1 )
            {
                double x      = xv[index_v_arr];
                double y      = yv[index_v_arr];
                double cosphi = cos( pideg * y );
                double zmult  = sin( pideg * y );
                double xmult  = cosphi * cos( x * pideg );
                double ymult  = cosphi * sin( x * pideg );
                Node3D pt( xmult, ymult, zmult );
                conn_arr[elem_index * nv + k] = vertex_map[pt];
            }
        }
        EntityHandle cell = start_cell + elem_index;
        // set other tags, like xc, yc, frac, area
        /*rval = mbImpl->tag_set_data( xcTag, &cell, 1, &xc[elem_index] );MB_CHK_SET_ERR( rval, "Failed to set xc tag" );
        rval = mbImpl->tag_set_data( ycTag, &cell, 1, &yc[elem_index] );MB_CHK_SET_ERR( rval, "Failed to set yc tag" );
        rval = mbImpl->tag_set_data( areaTag, &cell, 1, &area[elem_index] );MB_CHK_SET_ERR( rval, "Failed to set area tag" );
        rval = mbImpl->tag_set_data( fracTag, &cell, 1, &frac[elem_index] );MB_CHK_SET_ERR( rval, "Failed to set frac tag" );
*/
        // set the global id too:
        int globalId = localGidCells[elem_index];

        rval = mbImpl->tag_set_data( mGlobalIdTag, &cell, 1, &globalId );MB_CHK_SET_ERR( rval, "Failed to set global id tag" );
        if( gmId >= 0 )
        {
            int localMask = masks[elem_index];
            rval          = mbImpl->tag_set_data( maskTag, &cell, 1, &localMask );MB_CHK_SET_ERR( rval, "Failed to set mask tag" );
        }
    }

    rval = mbImpl->add_entities( _fileSet, tmp_range );MB_CHK_SET_ERR( rval, "Failed to add new cells to current file set" );

    // modify local file set, to merge coincident vertices, and to correct repeated vertices in elements
    std::vector< Tag > tagList;
    tagList.push_back( mGlobalIdTag );
    if( gmId >= 0 ) tagList.push_back( maskTag );
    rval = IntxUtils::remove_padded_vertices( mbImpl, _fileSet, tagList );MB_CHK_SET_ERR( rval, "Failed to remove duplicate vertices" );

    rval = mbImpl->get_entities_by_dimension( _fileSet, 2, faces );MB_CHK_ERR( rval );
    Range all_verts;
    rval = mbImpl->get_connectivity( faces, all_verts );MB_CHK_ERR( rval );
    rval = mbImpl->add_entities( _fileSet, all_verts );MB_CHK_ERR( rval );
#ifdef MOAB_HAVE_MPI
    ParallelComm*& myPcomm = _readNC->myPcomm;
    if( myPcomm )
    {
        double tol = 1.e-12;  // this is the same as static tolerance in NCHelper
        ParallelMergeMesh pmm( myPcomm, tol );
        rval = pmm.merge( _fileSet,
                          /* do not do local merge*/ false,
                          /*  2d cells*/ 2 );MB_CHK_SET_ERR( rval, "Failed to merge vertices in parallel" );

        // assign global ids only for vertices, cells have them fine
        rval = myPcomm->assign_global_ids( _fileSet, /*dim*/ 0 );MB_CHK_ERR( rval );
        // remove all sets, edges and vertices from the file set
        Range edges, vertices;
        rval = mbImpl->get_entities_by_dimension(_fileSet, 1, edges, /*recursive*/ true);MB_CHK_ERR( rval );
        rval = mbImpl->get_entities_by_dimension(_fileSet, 0, vertices, /*recursive*/ true);MB_CHK_ERR( rval );
        rval = mbImpl->remove_entities(_fileSet, edges);MB_CHK_ERR( rval );
        rval = mbImpl->remove_entities(_fileSet, vertices);MB_CHK_ERR( rval );

        Range intfSets = myPcomm->interface_sets();
        // empty intf sets
        rval = mbImpl->clear_meshset(intfSets);MB_CHK_ERR( rval );
        // delete the sets without shame :)
        //sets.merge(intfSets);
        //rval = myPcomm->delete_entities(sets);MB_CHK_ERR( rval ); // will also clean shared ents !
        rval = myPcomm->delete_entities(edges);MB_CHK_ERR( rval ); // will also clean shared ents !
    }
#else
    rval = mbImpl->remove_entities( _fileSet, all_verts );MB_CHK_ERR( rval );
#endif

    return MB_SUCCESS;
}

#ifdef MOAB_HAVE_MPI
ErrorCode NCHelperScrip::redistribute_local_cells( int start_cell_idx )
{
    // If possible, apply Zoltan partition
#ifdef MOAB_HAVE_ZOLTAN
    if( ScdParData::RCBZOLTAN == _readNC->partMethod )
    {
        // Read grid_center_lat coordinates of cell centers
        int xCellVarId;
        int success = NCFUNC( inq_varid )( _fileId, "grid_center_lon", &xCellVarId );
        if( success ) MB_SET_ERR( MB_FAILURE, "Failed to get variable id of grid_center_lon" );
        std::vector< double > xc( nLocalCells );
        NCDF_SIZE read_start = static_cast< NCDF_SIZE >( start_cell_idx - 1 );
        NCDF_SIZE read_count = static_cast< NCDF_SIZE >( nLocalCells );
        success              = NCFUNCAG( _vara_double )( _fileId, xCellVarId, &read_start, &read_count, &xc[0] );
        if( success ) MB_SET_ERR( MB_FAILURE, "Failed to read grid_center_lat data" );

        // Read grid_center_lon coordinates of cell centers
        int yCellVarId;
        success = NCFUNC( inq_varid )( _fileId, "grid_center_lat", &yCellVarId );
        if( success ) MB_SET_ERR( MB_FAILURE, "Failed to get variable id of grid_center_lat" );
        std::vector< double > yc( nLocalCells );
        success = NCFUNCAG( _vara_double )( _fileId, yCellVarId, &read_start, &read_count, &yc[0] );
        if( success ) MB_SET_ERR( MB_FAILURE, "Failed to read grid_center_lon data" );

        // Zoltan partition using RCB; maybe more studies would be good, as to which partition
        // is better
        Interface*& mbImpl         = _readNC->mbImpl;
        DebugOutput& dbgOut        = _readNC->dbgOut;
        ZoltanPartitioner* mbZTool = new ZoltanPartitioner( mbImpl, false, 0, NULL );
        std::vector< double > xCell( nLocalCells );
        std::vector< double > yCell( nLocalCells );
        std::vector< double > zCell( nLocalCells );
        double pideg = 1.;  // radians
        if( degrees ) pideg = acos( -1.0 ) / 180.0;
        double x, y, cosphi;
        for( int i = 0; i < nLocalCells; i++ )
        {
            x        = xc[i];
            y        = yc[i];
            cosphi   = cos( pideg * y );
            zCell[i] = sin( pideg * y );
            xCell[i] = cosphi * cos( x * pideg );
            yCell[i] = cosphi * sin( x * pideg );
        }
        ErrorCode rval = mbZTool->repartition( xCell, yCell, zCell, start_cell_idx, "RCB", localGidCells );MB_CHK_SET_ERR( rval, "Error in Zoltan partitioning" );
        delete mbZTool;

        dbgOut.tprintf( 1, "After Zoltan partitioning, localGidCells.psize() = %d\n", (int)localGidCells.psize() );
        dbgOut.tprintf( 1, "                           localGidCells.size() = %d\n", (int)localGidCells.size() );

        // This is important: local cells are now redistributed, so nLocalCells might be different!
        nLocalCells = localGidCells.size();

        return MB_SUCCESS;
    }
#endif

    // By default, apply trivial partition
    localGidCells.insert( start_cell_idx, start_cell_idx + nLocalCells - 1 );

    return MB_SUCCESS;
}
#endif

} /* namespace moab */
