/*
 * CSLAMUtils.hpp
 *
 */

#ifndef CSLAMUTILS_HPP_
#define CSLAMUTILS_HPP_

#include "moab/CartVect.hpp"
#include "moab/Core.hpp"
#include "moab/IntxMesh/IntxUtils.hpp"

class IntxUtilsCSLAM
{
  public:
    /* Analytical functions */

    // page 4 Nair Lauritzen paper
    // param will be: (la1, te1), (la2, te2), b, c; hmax=1, r=1/2
    static double quasi_smooth_field( double lam, double tet, double* params );

    // page 4
    static double smooth_field( double lam, double tet, double* params );

    // page 5
    static double slotted_cylinder_field( double lam, double tet, double* params );

    // More MBCSLAM specific API
    static void departure_point_case1( moab::CartVect& arrival_point, double t, double delta_t,
                                       moab::CartVect& departure_point );

    static void velocity_case1( moab::CartVect& arrival_point, double t, moab::CartVect& velo );

    // looking at DP tag, create the spanning quads, write a file (with rank) and
    // then delete the new entities (vertices) and the set of quads
    static moab::ErrorCode create_span_quads( moab::Interface* mb, moab::EntityHandle euler_set,
                                              int rank );

    // copy the euler mesh into a new set, lagr_set (or lagr set into a new euler set)
    // it will be used in 3rd method, when the positions of nodes are modified, no new nodes are
    //  created
    // it will also be used to
    static moab::ErrorCode deep_copy_set( moab::Interface* mb, moab::EntityHandle source,
                                          moab::EntityHandle dest );
};

#endif
