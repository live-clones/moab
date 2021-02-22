#include <cmath>
#include "moab/ReadUtilIface.hpp"
#include "IntxUtilsCSLAM.hpp"

using namespace moab;

#ifndef CORRTAGNAME
#define CORRTAGNAME "__correspondent"
#endif

// page 4 Nair Lauritzen paper
// param will be: (la1, te1), (la2, te2), b, c; hmax=1, r=1/2
double IntxUtilsCSLAM::quasi_smooth_field( double lam, double tet, double* params )
{
    double la1   = params[0];
    double te1   = params[1];
    double la2   = params[2];
    double te2   = params[3];
    double b     = params[4];
    double c     = params[5];
    double hmax  = params[6];  // 1;
    double r     = params[7];  // 0.5;
    double r1    = moab::IntxUtils::distance_on_sphere( lam, tet, la1, te1 );
    double r2    = moab::IntxUtils::distance_on_sphere( lam, tet, la2, te2 );
    double value = b;
    if( r1 < r ) { value += c * hmax / 2 * ( 1 + cos( M_PI * r1 / r ) ); }
    if( r2 < r ) { value += c * hmax / 2 * ( 1 + cos( M_PI * r2 / r ) ); }
    return value;
}

// page 4
// params are now x1, y1, ..., y2, z2 (6 params)
// plus h max and b0 (total of 8 params); radius is 1
double IntxUtilsCSLAM::smooth_field( double lam, double tet, double* params )
{
    moab::IntxUtils::SphereCoords sc;
    sc.R         = 1.;
    sc.lat       = tet;
    sc.lon       = lam;
    double hmax  = params[6];
    double b0    = params[7];
    CartVect xyz = moab::IntxUtils::spherical_to_cart( sc );
    CartVect c1( params );
    CartVect c2( params + 3 );
    double expo1 = -b0 * ( xyz - c1 ).length_squared();
    double expo2 = -b0 * ( xyz - c2 ).length_squared();
    return hmax * ( exp( expo1 ) + exp( expo2 ) );
}

// page 5
double IntxUtilsCSLAM::slotted_cylinder_field( double lam, double tet, double* params )
{
    double la1 = params[0];
    double te1 = params[1];
    double la2 = params[2];
    double te2 = params[3];
    double b   = params[4];
    double c   = params[5];
    // double hmax = params[6]; // 1;
    double r      = params[6];  // 0.5;
    double r1     = moab::IntxUtils::distance_on_sphere( lam, tet, la1, te1 );
    double r2     = moab::IntxUtils::distance_on_sphere( lam, tet, la2, te2 );
    double value  = b;
    double d1     = fabs( lam - la1 );
    double d2     = fabs( lam - la2 );
    double rp6    = r / 6;
    double rt5p12 = r * 5 / 12;

    if( r1 <= r && d1 >= rp6 ) value = c;
    if( r2 <= r && d2 >= rp6 ) value = c;
    if( r1 <= r && d1 < rp6 && tet - te1 < -rt5p12 ) value = c;
    if( r2 <= r && d2 < rp6 && tet - te2 > rt5p12 ) value = c;

    return value;
}

/*
 * based on paper A class of deformational flow test cases for linear transport problems on the
 * sphere longitude lambda [0.. 2*pi) and latitude theta (-pi/2, pi/2) lambda: -> lon (0, 2*pi)
 *  theta : -> lat (-pi/2, po/2)
 *  Radius of the sphere is 1 (if not, everything gets multiplied by 1)
 *
 *  cosine bell: center lambda0, theta0:
 */
void IntxUtilsCSLAM::departure_point_case1( CartVect& arrival_point, double t, double delta_t,
                                            CartVect& departure_point )
{
    // always assume radius is 1 here?
    moab::IntxUtils::SphereCoords sph = moab::IntxUtils::cart_to_spherical( arrival_point );
    double T                          = 5;    // duration of integration (5 units)
    double k                          = 2.4;  // flow parameter
    /*     radius needs to be within some range   */
    double sl2      = sin( sph.lon / 2 );
    double pit      = M_PI * t / T;
    double omega    = M_PI / T;
    double costetha = cos( sph.lat );
    // double u = k * sl2*sl2 * sin(2*sph.lat) * cos(pit);
    double v = k * sin( sph.lon ) * costetha * cos( pit );
    // double psi = k * sl2 * sl2 *costetha * costetha * cos(pit);
    double u_tilda = 2 * k * sl2 * sl2 * sin( sph.lat ) * cos( pit );

    // formula 35, page 8
    // this will approximate dep point using a Taylor series with up to second derivative
    // this will be O(delta_t^3) exact.
    double lon_dep =
        sph.lon - delta_t * u_tilda -
        delta_t * delta_t * k * sl2 *
            ( sl2 * sin( sph.lat ) * sin( pit ) * omega - u_tilda * sin( sph.lat ) * cos( pit ) * cos( sph.lon / 2 ) -
              v * sl2 * costetha * cos( pit ) );
    // formula 36, page 8 again
    double lat_dep = sph.lat - delta_t * v -
                     delta_t * delta_t / 4 * k *
                         ( sin( sph.lon ) * cos( sph.lat ) * sin( pit ) * omega -
                           u_tilda * cos( sph.lon ) * cos( sph.lat ) * cos( pit ) +
                           v * sin( sph.lon ) * sin( sph.lat ) * cos( pit ) );

    moab::IntxUtils::SphereCoords sph_dep;
    sph_dep.R   = 1.;  // radius
    sph_dep.lat = lat_dep;
    sph_dep.lon = lon_dep;

    departure_point = moab::IntxUtils::spherical_to_cart( sph_dep );
    return;
}

void IntxUtilsCSLAM::velocity_case1( CartVect& arrival_point, double t, CartVect& velo )
{
    // always assume radius is 1 here?
    moab::IntxUtils::SphereCoords sph = moab::IntxUtils::cart_to_spherical( arrival_point );
    double T                          = 5;    // duration of integration (5 units)
    double k                          = 2.4;  // flow parameter
    /*     radius needs to be within some range   */
    double sl2 = sin( sph.lon / 2 );
    double pit = M_PI * t / T;
    // double omega = M_PI/T;
    double coslat = cos( sph.lat );
    double sinlat = sin( sph.lat );
    double sinlon = sin( sph.lon );
    double coslon = cos( sph.lon );
    double u      = k * sl2 * sl2 * sin( 2 * sph.lat ) * cos( pit );
    double v      = k / 2 * sinlon * coslat * cos( pit );
    velo[0]       = -u * sinlon - v * sinlat * coslon;
    velo[1]       = u * coslon - v * sinlat * sinlon;
    velo[2]       = v * coslat;
}

ErrorCode IntxUtilsCSLAM::create_span_quads( Interface* mb, EntityHandle euler_set, int rank )
{
    // first get all edges adjacent to polygons
    Tag dpTag = 0;
    std::string tag_name( "DP" );
    ErrorCode rval = mb->tag_get_handle( tag_name.c_str(), 3, MB_TYPE_DOUBLE, dpTag, MB_TAG_DENSE );
    // if the tag does not exist, get out early
    if( rval != MB_SUCCESS ) return rval;
    Range polygons;
    rval = mb->get_entities_by_dimension( euler_set, 2, polygons );
    if( MB_SUCCESS != rval ) return rval;
    Range iniEdges;
    rval = mb->get_adjacencies( polygons, 1, false, iniEdges, Interface::UNION );
    if( MB_SUCCESS != rval ) return rval;
    // now create some if missing
    Range allEdges;
    rval = mb->get_adjacencies( polygons, 1, true, allEdges, Interface::UNION );
    if( MB_SUCCESS != rval ) return rval;
    // create the vertices at the DP points, and the quads after that
    Range verts;
    rval = mb->get_connectivity( polygons, verts );
    if( MB_SUCCESS != rval ) return rval;
    int num_verts = (int)verts.size();
    // now see the departure points; to what boxes should we send them?
    std::vector< double > dep_points( 3 * num_verts );
    rval = mb->tag_get_data( dpTag, verts, (void*)&dep_points[0] );
    if( MB_SUCCESS != rval ) return rval;

    // create vertices corresponding to dp locations
    ReadUtilIface* read_iface;
    rval = mb->query_interface( read_iface );
    if( MB_SUCCESS != rval ) return rval;
    std::vector< double* > coords;
    EntityHandle start_vert, start_elem, *connect;
    // create verts, num is 2(nquads+1) because they're in a 1d row; will initialize coords in loop
    // over quads later
    rval = read_iface->get_node_coords( 3, num_verts, 0, start_vert, coords );
    if( MB_SUCCESS != rval ) return rval;
    // fill it up
    // Cppcheck warning (false positive): variable coords is assigned a value that is never used
    for( int i = 0; i < num_verts; i++ )
    {
        // block from interleaved
        coords[0][i] = dep_points[3 * i];
        coords[1][i] = dep_points[3 * i + 1];
        coords[2][i] = dep_points[3 * i + 2];
    }
    // create quads; one quad for each edge
    rval = read_iface->get_element_connect( allEdges.size(), 4, MBQUAD, 0, start_elem, connect );
    if( MB_SUCCESS != rval ) return rval;

    const EntityHandle* edge_conn = NULL;
    int quad_index                = 0;
    EntityHandle firstVertHandle  = verts[0];  // assume vertices are contiguous...
    for( Range::iterator eit = allEdges.begin(); eit != allEdges.end(); ++eit, quad_index++ )
    {
        EntityHandle edge = *eit;
        int num_nodes;
        rval = mb->get_connectivity( edge, edge_conn, num_nodes );
        if( MB_SUCCESS != rval ) return rval;
        connect[quad_index * 4]     = edge_conn[0];
        connect[quad_index * 4 + 1] = edge_conn[1];

        // maybe some indexing in range?
        connect[quad_index * 4 + 2] = start_vert + edge_conn[1] - firstVertHandle;
        connect[quad_index * 4 + 3] = start_vert + edge_conn[0] - firstVertHandle;
    }

    Range quads( start_elem, start_elem + allEdges.size() - 1 );
    EntityHandle outSet;
    rval = mb->create_meshset( MESHSET_SET, outSet );
    if( MB_SUCCESS != rval ) return rval;
    mb->add_entities( outSet, quads );

    Tag colTag;
    rval = mb->tag_get_handle( "COLOR_ID", 1, MB_TYPE_INTEGER, colTag, MB_TAG_DENSE | MB_TAG_CREAT );
    if( MB_SUCCESS != rval ) return rval;
    int j = 1;
    for( Range::iterator itq = quads.begin(); itq != quads.end(); ++itq, j++ )
    {
        EntityHandle q = *itq;
        rval           = mb->tag_set_data( colTag, &q, 1, &j );
        if( MB_SUCCESS != rval ) return rval;
    }
    std::stringstream outf;
    outf << "SpanQuads" << rank << ".h5m";
    rval = mb->write_file( outf.str().c_str(), 0, 0, &outSet, 1 );
    if( MB_SUCCESS != rval ) return rval;
    EntityHandle outSet2;
    rval = mb->create_meshset( MESHSET_SET, outSet2 );
    if( MB_SUCCESS != rval ) return rval;

    Range quadEdges;
    rval = mb->get_adjacencies( quads, 1, true, quadEdges, Interface::UNION );
    if( MB_SUCCESS != rval ) return rval;
    mb->add_entities( outSet2, quadEdges );

    std::stringstream outf2;
    outf2 << "SpanEdges" << rank << ".h5m";
    rval = mb->write_file( outf2.str().c_str(), 0, 0, &outSet2, 1 );
    if( MB_SUCCESS != rval ) return rval;

    // maybe some clean up
    mb->delete_entities( &outSet, 1 );
    mb->delete_entities( &outSet2, 1 );
    mb->delete_entities( quads );
    Range new_edges = subtract( allEdges, iniEdges );
    mb->delete_entities( new_edges );
    new_edges = subtract( quadEdges, iniEdges );
    mb->delete_entities( new_edges );
    Range new_verts( start_vert, start_vert + num_verts );
    mb->delete_entities( new_verts );

    return MB_SUCCESS;
}

// this simply copies the one mesh set into another, and sets some correlation tags
// for easy mapping back and forth
ErrorCode IntxUtilsCSLAM::deep_copy_set( Interface* mb, EntityHandle source_set, EntityHandle dest_set )
{
    // create the handle tag for the corresponding element / vertex

    EntityHandle dum = 0;
    Tag corrTag      = 0;  // it will be created here
    ErrorCode rval   = mb->tag_get_handle( CORRTAGNAME, 1, MB_TYPE_HANDLE, corrTag, MB_TAG_DENSE | MB_TAG_CREAT, &dum );MB_CHK_ERR( rval );

    // give the same global id to new verts and cells created in the lagr(departure) mesh
    Tag gid = mb->globalId_tag();

    Range polys;
    rval = mb->get_entities_by_dimension( source_set, 2, polys );MB_CHK_ERR( rval );

    Range connecVerts;
    rval = mb->get_connectivity( polys, connecVerts );MB_CHK_ERR( rval );

    std::map< EntityHandle, EntityHandle > newNodes;
    for( Range::iterator vit = connecVerts.begin(); vit != connecVerts.end(); ++vit )
    {
        EntityHandle oldV = *vit;
        CartVect posi;
        rval = mb->get_coords( &oldV, 1, &( posi[0] ) );MB_CHK_ERR( rval );
        int global_id;
        rval = mb->tag_get_data( gid, &oldV, 1, &global_id );MB_CHK_ERR( rval );

        EntityHandle new_vert;
        // duplicate the position
        rval = mb->create_vertex( &( posi[0] ), new_vert );MB_CHK_ERR( rval );
        newNodes[oldV] = new_vert;

        // set also the correspondent tag :)
        rval = mb->tag_set_data( corrTag, &oldV, 1, &new_vert );MB_CHK_ERR( rval );
        // also the other side
        // need to check if we really need this; the new vertex will never need the old vertex
        // we have the global id which is the same
        rval = mb->tag_set_data( corrTag, &new_vert, 1, &oldV );MB_CHK_ERR( rval );
        // set the global id on the corresponding vertex the same as the initial vertex
        rval = mb->tag_set_data( gid, &new_vert, 1, &global_id );MB_CHK_ERR( rval );
    }

    for( Range::iterator it = polys.begin(); it != polys.end(); ++it )
    {
        EntityHandle q = *it;
        int nnodes;
        const EntityHandle* conn;
        rval = mb->get_connectivity( q, conn, nnodes );MB_CHK_ERR( rval );

        int global_id;
        rval = mb->tag_get_data( gid, &q, 1, &global_id );MB_CHK_ERR( rval );
        EntityType typeElem = mb->type_from_handle( q );
        std::vector< EntityHandle > new_conn( nnodes );
        for( int i = 0; i < nnodes; i++ )
        {
            EntityHandle v1 = conn[i];
            new_conn[i]     = newNodes[v1];
        }
        EntityHandle newElement;
        rval = mb->create_element( typeElem, &new_conn[0], nnodes, newElement );MB_CHK_ERR( rval );

        // set the corresponding tag; not sure we need this one, from old to new
        rval = mb->tag_set_data( corrTag, &q, 1, &newElement );MB_CHK_ERR( rval );
        rval = mb->tag_set_data( corrTag, &newElement, 1, &q );MB_CHK_ERR( rval );

        // set the global id
        rval = mb->tag_set_data( gid, &newElement, 1, &global_id );MB_CHK_ERR( rval );
        rval = mb->add_entities( dest_set, &newElement, 1 );MB_CHK_ERR( rval );
    }

    return MB_SUCCESS;
}
