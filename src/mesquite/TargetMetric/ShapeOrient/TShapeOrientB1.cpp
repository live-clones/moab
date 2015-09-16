/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

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
 
    (2006) kraftche@cae.wisc.edu
    (2010) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file TShapeOrientB1.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TShapeOrientB1.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"
#include "TMPCommon.hpp"

namespace MESQUITE_NS {

std::string TShapeOrientB1::get_name() const
  { return "TShapeOrientB1"; }

TShapeOrientB1::~TShapeOrientB1() {}

template <unsigned DIM> static inline
bool eval( const MsqMatrix<DIM,DIM>& T, double& result )
{
  const double tau = det(T);
  if (TMetric::invalid_determinant(tau)) { // barrier
    result = 0.0;
    return false;
  }
  result = 0.5/tau * (Frobenius( T ) - trace(T)/DimConst<DIM>::sqrt());
  return true;
}


template <unsigned DIM> static inline
bool grad( const MsqMatrix<DIM,DIM>& T, 
           double& result, 
           MsqMatrix<DIM,DIM>& deriv )
{
  const double norm = Frobenius(T);
  const double invroot = 1.0/DimConst<DIM>::sqrt();
  const double tau = det(T);
  if (TMetric::invalid_determinant(tau)) { // barrier
    result = 0.0;
    return false;
  }
  const double inv_tau = 1.0/tau;
  const double invnorm = 1.0/norm;
  
  result = 0.5*inv_tau*(norm - invroot * trace(T));

  deriv = invnorm * T;
  pluseq_scaled_I( deriv, -invroot );
  deriv *= 0.5;
  deriv -= result * transpose_adj(T);
  deriv *= inv_tau;
  return true;
}

template <unsigned DIM> static inline
bool hess( const MsqMatrix<DIM,DIM>& T, 
           double& result, 
           MsqMatrix<DIM,DIM>& deriv, 
           MsqMatrix<DIM,DIM>* second )
{
  const double norm = Frobenius(T);
  const double invroot = 1.0/DimConst<DIM>::sqrt();
  const double tau = det(T);
  if (TMetric::invalid_determinant(tau)) { // barrier
    result = 0.0;
    return false;
  }
  const double inv_tau = 1.0/tau;
  const double invnorm = 1.0/norm;
  
  const double f = norm - invroot * trace(T);
  result = 0.5 * inv_tau * f;

  const MsqMatrix<DIM,DIM> adjt = transpose_adj(T);
  deriv = invnorm * T;
  pluseq_scaled_I( deriv, -invroot );
  deriv *= 0.5;
  deriv -= result * adjt;
  deriv *= inv_tau;
  
  const double a = 0.5 * inv_tau * invnorm;
  set_scaled_outer_product( second, -a*invnorm*invnorm, T );
  pluseq_scaled_I( second, a );
  pluseq_scaled_outer_product( second, f*inv_tau*inv_tau*inv_tau, adjt );
  pluseq_scaled_2nd_deriv_of_det( second, -0.5*f*inv_tau*inv_tau, T );
  pluseq_scaled_sum_outer_product( second, -0.5*inv_tau*inv_tau*invnorm, T, adjt );
  pluseq_scaled_sum_outer_product_I( second, 0.5*inv_tau*inv_tau*invroot, adjt );
  return true;
}


TMP_T_TEMPL_IMPL_COMMON(TShapeOrientB1)


} // namespace Mesquite
