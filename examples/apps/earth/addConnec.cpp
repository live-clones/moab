/*
 *
 *
/** @example addConnec.cpp  Add connectivity to point cloud data, for better view in visit
 * this tool will take an existing h5m fine point cloud data (phys grid or land pc)
 *   and add 2d cells from a fine atm mesh
 *
 * example of usage:
 * ./addConnec -i wholeFineATM.h5m -s wholeLND_proj01.h5m -o wholeLndM.h5m
 *
 * Basically, will output a new h5m file (wholeLndM.h5m), which has also cell connectivity copied from fine atm
 *  matching is based on the global ids
 * between what we think is the order on the original file (wholeFineATM.h5m)
 *
 *  file  wholeFineATM.h5m is obtained from a coupled run in e3sm, with the ne 11, np 4,
 *
 */
#include "moab/ProgOptions.hpp"
#include "moab/Core.hpp"

#include <cmath>
#include <sstream>

using namespace moab;

int main( int argc, char* argv[] )
{

    ProgOptions opts;

    std::string inputfile, outfile( "out.h5m" ), sourcefile;

    opts.addOpt< std::string >( "input,i", "input mesh filename", &inputfile );
    opts.addOpt< std::string >( "source,s", "h5m file aligned with the mesh input file", &sourcefile );
    opts.addOpt< std::string >( "output,o", "output mesh filename", &outfile );

    opts.parseCommandLine( argc, argv );

    std::cout << "input mesh file: " << inputfile << "\n";
    std::cout << "source file: " << sourcefile << "\n";
    std::cout << "output file: " << outfile << "\n";

    if( inputfile.empty() )
    {
        opts.printHelp();
        return 0;
    }
    ErrorCode rval;
    Core* mb = new Core();

    rval = mb->load_file( inputfile.c_str() );MB_CHK_SET_ERR( rval, "can't load input file" );

    Core* mb2 = new Core();
    rval      = mb2->load_file( sourcefile.c_str() );MB_CHK_SET_ERR( rval, "can't load source file" );

    // get vertices on ini mesh; get global id on ini mesh
    // get global id on source mesh
    Tag gid;
    Tag gid2;
    rval = mb->tag_get_handle( "GLOBAL_ID", gid );MB_CHK_SET_ERR( rval, "can't get GLOBAL_ID tag on ini mesh " );
    rval = mb2->tag_get_handle( "GLOBAL_ID", gid2 );MB_CHK_SET_ERR( rval, "can't get GLOBAL_ID tag on source mesh " );

    // get vertices on ini mesh; build
    Range iniVerts;
    rval = mb->get_entities_by_dimension( 0, 0, iniVerts );MB_CHK_SET_ERR( rval, "can't get verts on initial mesh " );

    std::vector< int > gids;
    gids.resize( iniVerts.size() );
    rval = mb->tag_get_data( gid, iniVerts, &( gids[0] ) );MB_CHK_SET_ERR( rval, "can't get gid on initial verts " );
    // build now the map
    std::map< int, EntityHandle > fromGidToEh;
    int i = 0;
    for( Range::iterator vit = iniVerts.begin(); vit != iniVerts.end(); ++vit, i++ )
    {
        fromGidToEh[gids[i]] = *vit;
    }
    // now get the source verts, and tags, and set it on new mesh

    Range sourceVerts;
    rval = mb2->get_entities_by_dimension( 0, 0, sourceVerts );MB_CHK_SET_ERR( rval, "can't get verts on source mesh " );
    std::vector< int > gids2;
    gids2.resize( sourceVerts.size() );
    rval = mb2->tag_get_data( gid2, sourceVerts, &( gids2[0] ) );MB_CHK_SET_ERR( rval, "can't get gid2 on cloud mesh " );
    std::map< int, EntityHandle > fromGid2ToEh;
    i = 0;
    for( Range::iterator vit = sourceVerts.begin(); vit != sourceVerts.end(); ++vit, i++ )
    {
        fromGid2ToEh[gids2[i]] = *vit;
    }

    Range usedVerts;
    for( size_t i = 0; i < gids2.size(); i++ )
    {
        int globalId          = gids2[i];
        EntityHandle usedVert = fromGidToEh[globalId];
        usedVerts.insert( usedVert );
    }
    Range cells;
    rval = mb->get_adjacencies( usedVerts, 2, false, cells, Interface::UNION );MB_CHK_SET_ERR( rval, "can't get adj cells " );

    // now create a new cell in mb2 for each one in mb
    for( Range::iterator cit = cells.begin(); cit != cells.end(); ++cit )
    {
        EntityHandle cell        = *cit;
        int nnodes               = 0;
        const EntityHandle* conn = NULL;
        rval                     = mb->get_connectivity( cell, conn, nnodes );MB_CHK_SET_ERR( rval, "can't get connectivity " );

        std::vector< EntityHandle > nconn( nnodes );
        bool goodCell = true;
        for( int i = 0; i < nnodes; i++ )
        {
            EntityHandle v = conn[i];
            int index      = iniVerts.index( v );
            if( index < 0 ) goodCell = false;
            int id = gids[index];
            if( fromGid2ToEh.find( id ) == fromGid2ToEh.end() )
                goodCell = false;
            else
                nconn[i] = fromGid2ToEh[id];
        }
        if( !goodCell ) continue;
        EntityType type = MBTRI;
        if( nnodes == 3 )
            type = MBTRI;
        else if( nnodes == 4 )
            type = MBQUAD;
        else if( nnodes > 4 )
            type = MBPOLYGON;
        EntityHandle nel;
        rval = mb2->create_element( type, &nconn[0], nnodes, nel );MB_CHK_SET_ERR( rval, "can't create new cell " );
    }
    // save file
    rval = mb2->write_file( outfile.c_str() );MB_CHK_SET_ERR( rval, "can't write file" );

    delete mb;
    delete mb2;

    return 0;
}
