/*
 * IntxUtils.hpp
 *
 *  Created on: Oct 3, 2012
 */

#ifndef INTXUTILS_HPP_
#define INTXUTILS_HPP_

#include "moab/CartVect.hpp"
#include "moab/Core.hpp"

namespace moab
{

class IntxUtils
{
public:

    static double dist2(double * a, double * b);
    static double area2D(double *a, double *b, double *c);
    static int borderPointsOfXinY2(double * X, int nX, double * Y, int nY, double * P, int * side, double epsilon_area);
    static int SortAndRemoveDoubles2(double * P, int & nP, double epsilon);
    // the marks will show what edges of blue intersect the red

    static ErrorCode EdgeIntersections2(double * blue, int nsBlue, double * red, int nsRed,
        int * markb, int * markr, double * points, int & nPoints);

    // special one, for intersection between rll (constant latitude)  and cs quads
    static ErrorCode EdgeIntxRllCs(double * blue, CartVect * bluec, int * blueEdgeType, int nsBlue, double * red, CartVect * redc,
        int nsRed, int * markb, int * markr, int plane, double Radius, double * points, int & nPoints);

    // vec utils related to gnomonic projection on a sphere
    // vec utils
    /*
    *
    * position on a sphere of radius R
    * if plane specified, use it; if not, return the plane, and the point in the plane
    * there are 6 planes, numbered 1 to 6
    * plane 1: x=R, plane 2: y=R, 3: x=-R, 4: y=-R, 5: z=-R, 6: z=R
    *
    * projection on the plane will preserve the orientation, such that a triangle, quad pointing
    * outside the sphere will have a positive orientation in the projection plane
    * this is similar logic to /cesm1_0_4/models/atm/cam/src/dynamics/homme/share/coordinate_systems_mod.F90
    *    method: function cart2face(cart3D) result(face_no)
    */
    static void decide_gnomonic_plane(const CartVect & pos, int & oPlane);

    // point on a sphere is projected on one of six planes, decided earlier

    static ErrorCode gnomonic_projection(const CartVect & pos, double R, int plane, double & c1, double & c2);

    // given the position on plane (one out of 6), find out the position on sphere
    static ErrorCode reverse_gnomonic_projection(const double & c1, const double & c2, double R, int plane,
        CartVect & pos);

    /*
    *   other methods to convert from spherical coord to cartesian, and back
    *   A spherical coordinate triple is (R, lon, lat)
    *   should we store it as a CartVect? probably not ...
    *   /cesm1_0_4/models/atm/cam/src/dynamics/homme/share/coordinate_systems_mod.F90
    *
        enforce three facts:
        ! ==========================================================
        ! enforce three facts:
        !
        ! 1) lon at poles is defined to be zero
        !
        ! 2) Grid points must be separated by about .01 Meter (on earth)
        !    from pole to be considered "not the pole".
        !
        ! 3) range of lon is { 0<= lon < 2*pi }
        !
        ! 4) range of lat is from -pi/2 to pi/2; -pi/2 or pi/2 are the poles, so there the lon is 0
        !
        ! ==========================================================
    */
    struct SphereCoords{
        double R, lon, lat;
    };

    static SphereCoords cart_to_spherical(CartVect &) ;

    static CartVect spherical_to_cart (SphereCoords &) ;

    /*
    *  create a mesh used mainly for visualization for now, with nodes corresponding to
    *   GL points, a so-called refined mesh, with NP nodes in each direction.
    *   input: a range of quads (coarse), and  a desired order (NP is the number of points), so it
    *   is order + 1
    *
    *   output: a set with refined elements; with proper input, it should be pretty
    *   similar to a Homme mesh read with ReadNC
    */
    //ErrorCode SpectralVisuMesh(Interface * mb, Range & input, int NP, EntityHandle & outputSet, double tolerance);

    /*
    * given an entity set, get all nodes and project them on a sphere with given radius
    */
    static ErrorCode ScaleToRadius(Interface * mb, EntityHandle set, double R);

    static double distance_on_great_circle(CartVect & p1, CartVect & p2);

    static void departure_point_case1(CartVect & arrival_point, double t, double delta_t, CartVect & departure_point);

    static void velocity_case1(CartVect & arrival_point, double t, CartVect & velo);

    // break the nonconvex quads into triangles; remove the quad from the set? yes.
    // maybe radius is not needed;
    //
    static ErrorCode enforce_convexity(Interface * mb, EntityHandle set, int rank = 0);

    // this could be larger than PI, because of orientation; useful for non-convex polygons
    static double oriented_spherical_angle(double * A, double * B, double * C);

    // looking at DP tag, create the spanning quads, write a file (with rank) and
    // then delete the new entities (vertices) and the set of quads
    static ErrorCode create_span_quads(Interface * mb, EntityHandle euler_set, int rank);

    // looking at quad connectivity, collapse to triangle if 2 nodes equal
    // then delete the old quad
    static ErrorCode fix_degenerate_quads(Interface * mb, EntityHandle set);

    // distance along a great circle on a sphere of radius 1
    static double distance_on_sphere(double la1, double te1, double la2, double te2);

    /* Analytical functions */

    // page 4 Nair Lauritzen paper
    // param will be: (la1, te1), (la2, te2), b, c; hmax=1, r=1/2
    static double quasi_smooth_field(double lam, double tet, double * params);

    // page 4
    static double smooth_field(double lam, double tet, double * params);

    // page 5
    static double slotted_cylinder_field(double lam, double tet, double * params);

    /* End Analytical functions */

    /*
    * given 2 arcs AB and CD, compute the unique intersection point, if it exists
    *  in between
    */
    static ErrorCode intersect_great_circle_arcs(double * A, double * B, double * C, double * D, double R,
        double * E);
    /*
    * given 2 arcs AB and CD, compute the intersection points, if it exists
    *  AB is a great circle arc
    *  CD is a constant latitude arc
    */
    static ErrorCode intersect_great_circle_arc_with_clat_arc(double * A, double * B, double * C, double * D, double R,
        double * E, int & np);

    //ErrorCode  set_edge_type_flag(Interface * mb, EntityHandle sf1);

    static int  borderPointsOfCSinRLL(CartVect * redc, double * red2dc, int nsRed, CartVect *bluec, int nsBlue, int * blueEdgeType,
        double * P, int * side, double epsil);

    // copy the euler mesh into a new set, lagr_set (or lagr set into a new euler set)
    // it will be used in 3rd method, when the positions of nodes are modified, no new nodes are
    //  created
    // it will also be used to
    static ErrorCode  deep_copy_set(Interface * mb, EntityHandle source, EntityHandle dest);

    // used only by homme
    static ErrorCode  deep_copy_set_with_quads(Interface * mb, EntityHandle source_set, EntityHandle dest_set);

};

class IntxAreaUtils
{
private:
    /* data */
    bool use_lHuiller;

public:

    IntxAreaUtils(bool p_bUselHuiller = true) : use_lHuiller(p_bUselHuiller)
    {}

    ~IntxAreaUtils()
    {}

    /*
    * utilities to compute area of a polygon on which all edges are arcs of great circles on a sphere
    */
    /*
    * this will compute the spherical angle ABC, when A, B, C are on a sphere of radius R
    *  the radius will not be needed, usually, just for verification the points are indeed on that sphere
    *  the center of the sphere is at origin (0,0,0)
    *  this angle can be used in Girard's theorem to compute the area of a spherical polygon
    */
    double spherical_angle(double * A, double * B, double * C, double Radius);

    double area_spherical_triangle(double *A, double *B, double *C, double Radius);

    double area_spherical_polygon (double * A, int N, double Radius);

    double compute_area_spherical_triangle(double * A, double * B, double * C, double Radius = 1.0);

    double area_spherical_polygon_lHuiller (double * A, int N, double Radius);

    double area_spherical_polygon_lHuiller_check_sign(double * A, int N, double Radius, int * sign);

    double area_on_sphere(Interface * mb, EntityHandle set, double R);

    double area_on_sphere_lHuiller(Interface * mb, EntityHandle set, double R);

    double area_spherical_element(Interface * mb, EntityHandle  elem, double R);

    ErrorCode positive_orientation(Interface * mb, EntityHandle set, double R);

private:

#ifdef MOAB_HAVE_TEMPESTREMAP
    double area_spherical_triangle_GQ(double * ptA, double * ptB, double * ptC, double Radius);
#endif

    double area_spherical_triangle_lHuiller(double * ptA, double * ptB, double * ptC, double Radius);

};

}
#endif /* INTXUTILS_HPP_ */
