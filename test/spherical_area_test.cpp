/*
 * spherical_area_test.cpp
 *
 *  Created on: Feb 1, 2013
 */
#include <iostream>
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/IntxMesh/IntxUtils.hpp"
#include "TestUtil.hpp"

using namespace moab;

int main( int /* argc*/, char** /* argv[]*/ )
{
    // check command line arg
    const char* filename_mesh = STRINGIFY( MESHDIR ) "/mbcslam/eulerHomme.vtk";

    // read input mesh in a set
    Core moab;
    Interface* mb = &moab;  // global
    EntityHandle sf;
    ErrorCode rval = mb->create_meshset( MESHSET_SET, sf );
    if( MB_SUCCESS != rval ) return 1;

    rval = mb->load_file( filename_mesh, &sf );
    if( MB_SUCCESS != rval ) return 1;

    double R = 6.;  // should be input
    // compare total area with 4*M_PI * R^2

    const double area_sphere = R * R * M_PI * 4.;
    std::cout << "total area of the sphere        :  " << area_sphere << "\n";

    {
        moab::IntxAreaUtils areaAdaptor( moab::IntxAreaUtils::Girard );  // use_lHuiller = true
        double area1 = areaAdaptor.area_on_sphere( mb, sf, R, 0);
        std::cout << "total area with Girard          :  " << area1
                  << " rel error:" << fabs( ( area1 - area_sphere ) / area_sphere ) << "\n";
    }

    {
        moab::IntxAreaUtils areaAdaptor( moab::IntxAreaUtils::lHuiller );
        double area2 = areaAdaptor.area_on_sphere( mb, sf, R, 0 );
        std::cout << "total area with l'Huiller       : " << area2
                  << " rel error:" << fabs( ( area2 - area_sphere ) / area_sphere ) << "\n";
    }

#ifdef MOAB_HAVE_TEMPESTREMAP
    {
        moab::IntxAreaUtils areaAdaptor( moab::IntxAreaUtils::GaussQuadrature );
        double area3 = areaAdaptor.area_on_sphere( mb, sf, R, 0 );
        std::cout << "total area with GaussQuadrature : " << area3
                  << " rel error:" << fabs( ( area3 - area_sphere ) / area_sphere ) << "\n";
    }
#endif

    return 0;
}
