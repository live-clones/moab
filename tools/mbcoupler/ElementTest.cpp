#include "TestUtil.hpp"
#include "ElemUtil.hpp"
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include <iostream>

#ifndef MESHDIR
#error Specify MESHDIR to compile test
#endif

using namespace moab;

void test_tet();
void test_hex();
void test_spectral_hex();
void test_spectral_quad();
void test_spherical_quad();
void test_linear_tri();
void test_spherical_tri();

int main()
{
    int rval = 0;
    rval += RUN_TEST( test_tet );
    rval += RUN_TEST( test_hex );
    rval += RUN_TEST( test_spectral_hex );
    rval += RUN_TEST( test_spectral_quad );
    rval += RUN_TEST( test_spherical_quad );
    rval += RUN_TEST( test_linear_tri );
    rval += RUN_TEST( test_spherical_tri );
    return rval;
}

void test_tet()
{
    moab::Element::LinearTet tet;
}  // test_tet()

void test_hex()
{
    double positions[] = { 236.80706050970281, -139.55422526228017, 193.27999999999997,
                           236.47511729348639, -141.33020962638582, 193.27999999999997,
                           237.8457938295229,  -142.57076074835663, 193.27999999999997,
                           239.12702305519684, -139.96608933577852, 193.27999999999997,
                           236.80841321361444, -139.55341321335499, 202.654,
                           236.47655014713746, -141.32980272396816, 202.654,
                           237.8477913707564,  -142.57047282187165, 202.654,
                           239.12865103844533, -139.96531051891105, 202.654 };
    CartVect x( 235.96518686964933, -142.43503000077749, 188.19999999999987 );
    std::vector< CartVect > vertices;
    for( int i = 0; i < 8; i++ )
        vertices.push_back( CartVect( positions + 3 * i ) );

    moab::Element::LinearHex hex( vertices );
    double tol( 0.0001 );
    if( hex.inside_box( x, tol ) )
    {
        CartVect nat_par = hex.ievaluate( x, 0.0001 );
        std::cout << nat_par << "\n";
    }

    double positions2[] = { 49.890500000000024, -20.376134375374882, 312.72000000000003,
                            52.015875000000044, -19.149048546996006, 312.72000000000003,
                            48.430375821458099, -18.548796774572125, 312.72000000000003,
                            47.717616239031223, -21.191360829777231, 312.72000000000003,
                            49.890500000000024, -20.376134375374882, 322.88,
                            52.015875000000044, -19.149048546996006, 322.88,
                            48.429930354643275, -18.52828610485021,  322.88,
                            47.720552036968819, -21.167591146685712, 322.88 };

    CartVect x2( 51.469000000000015, -20.145482942833631, 317.80000000000001 );

    vertices.clear();
    for( int i = 0; i < 8; i++ )
        vertices.push_back( CartVect( positions2 + 3 * i ) );
    moab::Element::LinearHex hex2( vertices );
    if( hex2.inside_box( x2, tol ) )
    {
        try
        {
            CartVect nat_par = hex.ievaluate( x, 0.0001 );
            std::cout << nat_par << "\n";
        }
        catch( Element::Map::EvaluationError )
        {
            // nothing
        }
    }

}  // test_hex()

void test_spectral_hex()
{
    // first load a model that has spectral elements
    moab::Core* mb       = new moab::Core();
    std::string meshFile = STRINGIFY( MESHDIR ) "/spectral.h5m";
    moab::ErrorCode rval = mb->load_mesh( meshFile.c_str() );
    if( moab::MB_SUCCESS != rval )
    {
        std::cout << "Problems reading file " << meshFile << "." << std::endl;
        delete mb;
        return;
    }

    // get the ent set with SEM_DIMS tag
    moab::Range spectral_sets;
    moab::Tag sem_tag;
    rval = mb->tag_get_handle( "SEM_DIMS", 3, moab::MB_TYPE_INTEGER, sem_tag );
    if( moab::MB_SUCCESS != rval )
    {
        std::cout << "can't find tag, no spectral set\n";
        delete mb;
        return;
    }
    rval = mb->get_entities_by_type_and_tag( 0, moab::MBENTITYSET, &sem_tag, NULL, 1, spectral_sets );
    if( moab::MB_SUCCESS != rval || spectral_sets.empty() )
    {
        std::cout << "can't get sem set\n";
        delete mb;
        return;
    }

    moab::Range ents;

    int sem_dims[3];
    moab::EntityHandle firstSemSet = spectral_sets[0];
    rval                           = mb->tag_get_data( sem_tag, &firstSemSet, 1, (void*)sem_dims );
    if( moab::MB_SUCCESS != rval )
    {
        delete mb;
        return;
    }

    rval = mb->get_entities_by_dimension( firstSemSet, 3, ents );
    if( moab::MB_SUCCESS != rval )
    {
        delete mb;
        return;
    }
    std::cout << "Found " << ents.size() << " " << 3 << "-dimensional entities:" << std::endl;

    if( sem_dims[0] != sem_dims[1] || sem_dims[0] != sem_dims[2] )
    {
        std::cout << " dimensions are different. bail out\n";
        delete mb;
        return;
    }

    // get the SEM_X ...tags
    moab::Tag xm1Tag, ym1Tag, zm1Tag;
    int ntot = sem_dims[0] * sem_dims[1] * sem_dims[2];
    rval     = mb->tag_get_handle( "SEM_X", ntot, moab::MB_TYPE_DOUBLE, xm1Tag );
    if( moab::MB_SUCCESS != rval )
    {
        std::cout << "can't get xm1tag \n";
        delete mb;
        return;
    }
    rval = mb->tag_get_handle( "SEM_Y", ntot, moab::MB_TYPE_DOUBLE, ym1Tag );
    if( moab::MB_SUCCESS != rval )
    {
        std::cout << "can't get ym1tag \n";
        delete mb;
        return;
    }
    rval = mb->tag_get_handle( "SEM_Z", ntot, moab::MB_TYPE_DOUBLE, zm1Tag );
    if( moab::MB_SUCCESS != rval )
    {
        std::cout << "can't get zm1tag \n";
        delete mb;
        return;
    }
    moab::Tag velTag;

    rval = mb->tag_get_handle( "VX", ntot, moab::MB_TYPE_DOUBLE, velTag );
    if( moab::MB_SUCCESS != rval )
    {
        std::cout << "can't get veltag \n";
        delete mb;
        return;
    }
    moab::Element::SpectralHex specHex( sem_dims[0] );

    // compute the data for some elements
    for( moab::Range::iterator rit = ents.begin(); rit != ents.end(); ++rit )
    {
        // get the tag pointers to the internal storage for xm1, to not copy the values
        moab::EntityHandle eh = *rit;
        const double* xval;
        const double* yval;
        const double* zval;
        rval = mb->tag_get_by_ptr( xm1Tag, &eh, 1, (const void**)&xval );
        if( moab::MB_SUCCESS != rval )
        {
            std::cout << "can't get xm1 values \n";
            delete mb;
            return;
        }
        rval = mb->tag_get_by_ptr( ym1Tag, &eh, 1, (const void**)&yval );
        if( moab::MB_SUCCESS != rval )
        {
            std::cout << "can't get ym1 values \n";
            delete mb;
            return;
        }
        rval = mb->tag_get_by_ptr( zm1Tag, &eh, 1, (const void**)&zval );
        if( moab::MB_SUCCESS != rval )
        {
            std::cout << "can't get zm1 values \n";
            delete mb;
            return;
        }
        if( rit == ents.begin() )
        {
            std::cout << " xm1 for first element: \n";
            for( int i = 0; i < ntot; i++ )
                std::cout << " " << xval[i];
            std::cout << "\n";
        }
        specHex.set_gl_points( (double*)xval, (double*)yval, (double*)zval );
        // first evaluate a point, then inverse it to see if we get the same thing
        moab::CartVect rst( 0.1, -0.1, 0.5 );
        moab::CartVect pos     = specHex.evaluate( rst );
        moab::CartVect inverse = specHex.ievaluate( pos, 0.0001 );
        std::cout << "difference" << rst - inverse << "\n";
        Matrix3 jac = specHex.jacobian( rst );
        std::cout << "jacobian: \n" << jac << " \n determinant: " << jac.determinant() << "\n";
        // evaluate vx at some point
        const double* vx;
        rval = mb->tag_get_by_ptr( velTag, &eh, 1, (const void**)&vx );
        if( moab::MB_SUCCESS != rval )
        {
            std::cout << "can't get vel values \n";
            delete mb;
            return;
        }
        double vel1 = specHex.evaluate_scalar_field( rst, vx );
        std::cout << "velocity: " << vel1 << "\n";
        // compute integral over vx:
        double integral = specHex.integrate_scalar_field( vx );
        std::cout << "integral over vx: " << integral << "\n";
    }
    std::cout << "success...\n";

    delete mb;
}

void test_spectral_quad()
{
    // first load a model that has spectral elements
    moab::Core* mb = new moab::Core();
    // use the grid on Sphere from mbcslam folder
    std::string meshFile = STRINGIFY( MESHDIR ) "/mbcslam/eulerHomme.vtk";
    moab::ErrorCode rval = mb->load_mesh( meshFile.c_str() );
    if( moab::MB_SUCCESS != rval )
    {
        std::cout << "Problems reading file " << meshFile << "." << std::endl;
        delete mb;
        return;
    }

    // for each element, compute the gl points and project them on sphere
    // the radius is the same as the one from intersection test on sphere
    // double R = 6. * sqrt(3.) / 2; // input

    moab::Range ents;

    rval = mb->get_entities_by_type( 0, moab::MBQUAD, ents );  // get all quads
    if( moab::MB_SUCCESS != rval )
    {
        delete mb;
        return;
    }
    std::cout << "Found " << ents.size() << " " << 2 << "-dimensional entities:" << std::endl;

    //
    int NP = 4;  // test this....
    moab::Element::SpectralQuad specQuad( NP );

    // compute the gl points for some elements
    for( moab::Range::iterator rit = ents.begin(); rit != ents.end(); ++rit )
    {

        const moab::EntityHandle* conn4 = NULL;
        int num_nodes                   = 0;
        rval                            = mb->get_connectivity( *rit, conn4, num_nodes );
        if( moab::MB_SUCCESS != rval )
        {
            std::cout << "can't get connectivity for quad \n";
            delete mb;
            return;
        }
        assert( num_nodes == 4 );

        std::vector< CartVect > verts( 4 );
        rval = mb->get_coords( conn4, num_nodes, &( verts[0][0] ) );
        if( moab::MB_SUCCESS != rval )
        {
            std::cout << "can't get coords for quad \n";
            delete mb;
            return;
        }

        specQuad.set_vertices( verts );
        specQuad.compute_gl_positions();
        // do something with the gl positions, project them on a sphere, and create another mesh?
        if( rit == ents.begin() )
        {
            std::cout << " gl points for first element: \n";
            int size;
            double* xyz[3];
            specQuad.get_gl_points( xyz[0], xyz[1], xyz[2], size );
            for( int i = 0; i < size; i++ )
                std::cout << xyz[0][i] << " " << xyz[1][i] << " " << xyz[2][i] << "\n";
        }

        // project them on a sphere, and create another mesh with it?
    }
    std::cout << "success...\n";

    delete mb;
}
void test_spherical_quad()
{
    // example from one coupler test, run like this
    // ./mbcoupler_test -meshes sphere_16p.h5m mpas_p8.h5m -itag vertex_field -meth 4 -outfile
    // dd.h5m method 4 is spherical
    double positions[] = { -0.88882388032987436, -0.069951956448441419, 0.45287838714646161,   -0.88226455385461389,
                           -0.13973697758043971, 0.4495362433757738,    -0.84497006020160348,  -0.13383011007602069,
                           0.51779831884618843,  -0.85072691325794214,  -0.066953660115039074, 0.52132612293631853 };
    CartVect x( -0.85408569769999998, -0.12391301439999999, 0.50515659540000002 );
    std::vector< CartVect > vertices;
    for( int i = 0; i < 4; i++ )
        vertices.push_back( CartVect( positions + 3 * i ) );

    moab::Element::SphericalQuad squad( vertices );
    double tol( 0.0001 );
    if( squad.inside_box( x, tol ) )
    {
        CartVect nat_par = squad.ievaluate( x, 0.000001 );
        std::cout << nat_par << "\n";
    }
    std::cout << "success...\n";
}

void test_linear_tri()
{
    double positions[] = { 0, 0, 0, 2, 0, 0, 0, 3, 0 };
    CartVect x( 1, 0.5, 0 );
    std::vector< CartVect > vertices;
    for( int i = 0; i < 3; i++ )
        vertices.push_back( CartVect( positions + 3 * i ) );

    moab::Element::LinearTri tri( vertices );
    double tol( 0.0001 );
    if( tri.inside_box( x, tol ) )
    {
        CartVect nat_par = tri.ievaluate( x );
        std::cout << x << " :" << nat_par << "\n";
    }

    x = CartVect( 0, 2, 0 );
    if( tri.inside_box( x, tol ) )
    {
        CartVect nat_par = tri.ievaluate( x );
        std::cout << x << " :" << nat_par << "\n";
    }

    x = CartVect( 1, 0, 0.5 );
    if( tri.inside_box( x, tol ) )
    {
        CartVect nat_par = tri.ievaluate( x );
        std::cout << x << " :" << nat_par << "\n";
    }

    double positions2[] = { -0.866026, -0.500001, 0., 0.866026, -0.500001, 0., 0.000000, 100.000000, 0. };
    x                   = CartVect( 0, 0, 0 );
    std::vector< CartVect > vertices2;
    for( int i = 0; i < 3; i++ )
        vertices2.push_back( CartVect( positions2 + 3 * i ) );

    moab::Element::LinearTri tri2( vertices2 );

    if( tri2.inside_box( x, tol ) )
    {
        CartVect nat_par = tri2.ievaluate( x );
        std::cout << x << " :" << nat_par << "\n";
    }

    std::cout << "vertices2 " << vertices2[0] << " " << vertices2[1] << " " << vertices2[2] << "\n";

    x = CartVect( -0.866026, -0.500001, 0. );
    std::cout << x << " :" << tri2.ievaluate( x ) << "\n";

    x = CartVect( +0.866026, -0.500001, 0. );
    std::cout << x << " :" << tri2.ievaluate( x ) << "\n";
    x = CartVect( 0.000000, 100.000000, 0. );
    std::cout << x << " :" << tri2.ievaluate( x ) << "\n";

    std::cout << "success...\n";
}

void test_spherical_tri()
{
    // example from one coupler test, run like this
    // ./mbcoupler_test -meshes  tri_fl_8p.h5m mpas_p8.h5m -itag vertex_field -meth 4  -outfile
    // oo.h5m -eps 1.e-9 method 4 is spherical
    double positions[] = { -0.86339258282987197, -0.17004443185241255,  0.47501383044112816,
                           -0.80777478326268271, -0.15172299908552511,  0.5696314870803928,
                           -0.8655618847392077,  -0.061613422011313854, 0.49699739427361828 };
    CartVect x( -0.85408569769999998, -0.12391301439999999, 0.50515659540000002 );
    std::vector< CartVect > vertices;
    for( int i = 0; i < 3; i++ )
        vertices.push_back( CartVect( positions + 3 * i ) );

    moab::Element::SphericalTri sphtri( vertices );
    double tol( 0.000001 );
    if( sphtri.inside_box( x, tol ) )
    {
        CartVect nat_par = sphtri.ievaluate( x, 0.000001 );
        std::cout << nat_par << "\n";
    }
    std::cout << "success...\n";
}
