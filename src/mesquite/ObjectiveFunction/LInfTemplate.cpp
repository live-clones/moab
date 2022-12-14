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
/*!
  \file   LInfTemplate.cpp
  \brief

  This Objective Function is evaluated using an L-infinity norm.
  total=max (abs(x))
  \author Michael Brewer
  \date   2002-07-3
*/
#include <cmath>
#include "LInfTemplate.hpp"
#include "QualityMetric.hpp"
#include "MsqError.hpp"

namespace MBMesquite
{

LInfTemplate::LInfTemplate( QualityMetric* qualitymetric ) : ObjectiveFunctionTemplate( qualitymetric ) {}

// Michael:  need to clean up here
LInfTemplate::~LInfTemplate() {}

ObjectiveFunction* LInfTemplate::clone() const
{
    return new LInfTemplate( get_quality_metric() );
}

void LInfTemplate::clear() {}

bool LInfTemplate::evaluate( EvalType type, PatchData& pd, double& value_out, bool free, MsqError& err )
{
    if( type != ObjectiveFunction::CALCULATE )
    {
        MSQ_SETERR( err )
        ( "LInfTemplate does not support block coodinate descent algoritms", MsqError::INVALID_STATE );
        return false;
    }

    QualityMetric* qm = get_quality_metric();
    qm->get_evaluations( pd, qmHandles, free, err );
    MSQ_ERRFALSE( err );
    const double negate = qm->get_negate_flag();

    // calculate OF value for just the patch
    std::vector< size_t >::const_iterator i;
    double value;
    value_out = -HUGE_VAL;
    for( i = qmHandles.begin(); i != qmHandles.end(); ++i )
    {
        bool result = qm->evaluate( pd, *i, value, err );
        if( MSQ_CHKERR( err ) || !result ) return false;

        value = negate * fabs( value );
        if( value > value_out ) value_out = value;
    }

    return true;
}

}  // namespace MBMesquite
