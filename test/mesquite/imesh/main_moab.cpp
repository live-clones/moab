/* *****************************************************************
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000
    with Sandia Corporation, the U.S. Government retains certain
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov,
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov

  ***************************************************************** */
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3
// -*-
//
//   SUMMARY:
//     USAGE:
//
// ORIG-DATE: 19-Feb-02 at 10:57:52
//  LAST-MOD: 10-Feb-04 at 22:44:58 by Thomas Leurent
//
//
// DESCRIPTION:
// ============
/*! \file main.cpp

describe main.cpp here

 */
// DESCRIP-END.
//

#include "Mesquite.hpp"
#include "MsqMOAB.hpp"
#include "MeshImpl.hpp"
#include "MsqError.hpp"
#include "InstructionQueue.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"
#include "PlanarDomain.hpp"
#include "MeshWriter.hpp"
#include "TestUtil.hpp"

// algorithms
#include "IdealWeightInverseMeanRatio.hpp"
#include "EdgeLengthQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "FeasibleNewton.hpp"
#include "ConjugateGradient.hpp"
#include "SmartLaplacianSmoother.hpp"

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

#include "iBase.h"

using namespace MBMesquite;

std::string default_file_name = TestDir + "unittest/mesquite/3D/vtk/large_box_hex_1000.vtk";

void usage()
{
    cout << "main [-N] [filename]" << endl;
    cout << "  -N : Use native representation instead of TSTT implementation\n";
    cout << "  If no file name is specified, will use \"" << default_file_name << '"' << endl;
    exit( 1 );
}

// Construct a MeshTSTT from the file
Mesh* get_imesh_mesh( const char* file_name );

// Construct a MeshImpl from the file
Mesh* get_native_mesh( const char* file_name );

// Run FeasibleNewton solver
int run_global_smoother( Mesh* mesh, MsqError& err );

// Run SmoothLaplacian solver
int run_local_smoother( Mesh* mesh, MsqError& err );

int main( int argc, char* argv[] )
{
    MBMesquite::MsqPrintError err( cout );

    // command line arguments
    const char* file_name = 0;
    bool use_native = false, opts_done = false;
    for( int arg = 1; arg < argc; ++arg )
    {
        if( !opts_done && argv[arg][0] == '-' )
        {
            if( !strcmp( argv[arg], "-N" ) )
                use_native = true;
            else if( !strcmp( argv[arg], "--" ) )
                opts_done = true;
            else
                usage();
        }
        else if( !file_name )
            file_name = argv[arg];
        else
            usage();
    }
    if( !file_name )
    {
        file_name = default_file_name.c_str();
        cout << "No file specified: using default: " << default_file_name << endl;
    }

    // Try running a global smoother on the mesh
    Mesh* mesh = use_native ? get_native_mesh( file_name ) : get_imesh_mesh( file_name );
    if( !mesh )
    {
        std::cerr << "Failed to load input file.  Aborting." << std::endl;
        return 1;
    }

    MeshWriter::write_vtk( mesh, "original.vtk", err );
    if( err ) return 1;
    cout << "Wrote \"original.vtk\"" << endl;
    run_global_smoother( mesh, err );
    if( err ) return 1;

    // Try running a local smoother on the mesh
    mesh = use_native ? get_native_mesh( file_name ) : get_imesh_mesh( file_name );
    if( !mesh )
    {
        std::cerr << "Failed to load input file.  Aborting." << std::endl;
        return 1;
    }

    run_local_smoother( mesh, err );
    if( err ) return 1;

    return 0;
}

int run_global_smoother( Mesh* mesh, MsqError& err )
{
    double OF_value = 0.0001;

    // creates an intruction queue
    InstructionQueue queue1;

    // creates a mean ratio quality metric ...
    IdealWeightInverseMeanRatio* mean_ratio = new IdealWeightInverseMeanRatio( err );
    if( err ) return 1;
    mean_ratio->set_averaging_method( QualityMetric::SUM, err );
    if( err ) return 1;

    // ... and builds an objective function with it
    LPtoPTemplate* obj_func = new LPtoPTemplate( mean_ratio, 1, err );
    if( err ) return 1;

    // creates the feas newt optimization procedures
    FeasibleNewton* pass1 = new FeasibleNewton( obj_func );
    pass1->use_global_patch();
    if( err ) return 1;

    QualityAssessor stop_qa( mean_ratio );

    // **************Set stopping criterion****************
    TerminationCriterion tc_inner;
    tc_inner.add_absolute_vertex_movement( OF_value );
    if( err ) return 1;
    TerminationCriterion tc_outer;
    tc_outer.add_iteration_limit( 1 );
    pass1->set_inner_termination_criterion( &tc_inner );
    pass1->set_outer_termination_criterion( &tc_outer );

    queue1.add_quality_assessor( &stop_qa, err );
    if( err ) return 1;

    // adds 1 pass of pass1 to mesh_set1
    queue1.set_master_quality_improver( pass1, err );
    if( err ) return 1;

    queue1.add_quality_assessor( &stop_qa, err );
    if( err ) return 1;

    // launches optimization on mesh_set
    queue1.run_instructions( mesh, err );
    if( err ) return 1;

    MeshWriter::write_vtk( mesh, "feasible-newton-result.vtk", err );
    if( err ) return 1;
    cout << "Wrote \"feasible-newton-result.vtk\"" << endl;

    // print_timing_diagnostics( cout );
    return 0;
}

int run_local_smoother( Mesh* mesh, MsqError& err )
{
    double OF_value = 0.0001;

    // creates an intruction queue
    InstructionQueue queue1;

    // creates a mean ratio quality metric ...
    IdealWeightInverseMeanRatio* mean_ratio = new IdealWeightInverseMeanRatio( err );
    if( err ) return 1;
    mean_ratio->set_averaging_method( QualityMetric::SUM, err );
    if( err ) return 1;

    // ... and builds an objective function with it
    LPtoPTemplate* obj_func = new LPtoPTemplate( mean_ratio, 1, err );
    if( err ) return 1;

    // creates the smart laplacian optimization procedures
    SmartLaplacianSmoother* pass1 = new SmartLaplacianSmoother( obj_func );

    QualityAssessor stop_qa( mean_ratio );

    // **************Set stopping criterion****************
    TerminationCriterion tc_inner;
    tc_inner.add_absolute_vertex_movement( OF_value );
    TerminationCriterion tc_outer;
    tc_outer.add_iteration_limit( 1 );
    pass1->set_inner_termination_criterion( &tc_inner );
    pass1->set_outer_termination_criterion( &tc_outer );

    queue1.add_quality_assessor( &stop_qa, err );
    if( err ) return 1;

    // adds 1 pass of pass1 to mesh_set
    queue1.set_master_quality_improver( pass1, err );
    if( err ) return 1;

    queue1.add_quality_assessor( &stop_qa, err );
    if( err ) return 1;

    // launches optimization on mesh_set
    queue1.run_instructions( mesh, err );
    if( err ) return 1;

    MeshWriter::write_vtk( mesh, "smart-laplacian-result.vtk", err );
    if( err ) return 1;
    cout << "Wrote \"smart-laplacian-result.vtk\"" << endl;

    // print_timing_diagnostics( cout );
    return 0;
}

Mesh* get_imesh_mesh( const char* file_name )
{
    moab::Core* mb = new( std::nothrow ) moab::Core;
    if( NULL == mb ) return 0;

    moab::ErrorCode rval;
    // This file is in the mesh files directory
    rval = mb->load_file( file_name );MB_CHK_SET_ERR_RET_VAL( rval, "Failed to read", 0 );

    moab::Tag fixed_tag;
    rval = mb->tag_get_handle( "fixed", fixed_tag );

    moab::EntityHandle root_set = 0;
    MsqError err;
    Mesh* result = new MBMesquite::MsqMOAB( mb, root_set, moab::MBHEX, err, &fixed_tag );
    if( MSQ_CHKERR( err ) )
    {
        delete result;
        cerr << err << endl;
        return 0;
    }

    return result;
}

Mesh* get_native_mesh( const char* file_name )
{
    MsqError err;
    MeshImpl* mesh = new MeshImpl;
    if( !mesh )
    {
        cerr << "Failed during MeshImpl construction.\n";
        exit( 2 );
    }
    mesh->read_vtk( file_name, err );
    if( err )
    {
        cerr << err << endl;
        exit( 3 );
    }

    return mesh;
}
