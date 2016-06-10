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
   
  ***************************************************************** */


/** \file ElementSolIndQM.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "ElementSolIndQM.hpp"
#include "MsqError.hpp"
#include "PatchData.hpp"
#include <limits>

namespace MBMesquite {

ElementSolIndQM::ElementSolIndQM(std::vector<double>& solution_indicator)
    : AveragingQM(QualityMetric::LINEAR), m_solution_indicator(solution_indicator)
    {}
    
ElementSolIndQM::~ElementSolIndQM() {}

std::string ElementSolIndQM::get_name() const
{
  return "ElementSolIndQM(UserIndicator)";
}


int ElementSolIndQM::get_negate_flag() const
{
  return 1;
}

bool ElementSolIndQM::evaluate( PatchData& pd, 
                              size_t handle, 
                              double& value, 
                              MsqError& err )
{
  bool return_flag=true;
  // double met_vals[MSQ_MAX_NUM_VERT_PER_ENT]; // 27 ?
  value=MSQ_MAX_CAP;
  
  const MsqMeshEntity* element = &pd.element_by_index(handle);
  const size_t* v_i = element->get_vertex_index_array();
  const MsqVertex* vertices = pd.get_vertex_array();

  if (handle < m_solution_indicator.size()) {
    Vector3D centroid;
    element->get_centroid(centroid, pd, err);
    const double CONST_PI=3.1415926535897932;
    // Specific to triangle
    // double area = 0.5*((vertices[v_i[1]]-vertices[v_i[0]])*(vertices[v_i[2]]-vertices[v_i[0]])).length();
    // value = area * (sin(CONST_PI*(centroid.x()+0.5))*sin(CONST_PI*(0.5+centroid.y())));
    // if (centroid.x()*centroid.x()+centroid.y()*centroid.y()-0.25 < 0)
    // if (!(centroid.x()>-0.25 && centroid.x() <= 0.25) && (centroid.y() > -0.25 && centroid.y() <= 0.25)) {
    if (!(centroid.y() > -1.25 && centroid.y() <= 1.25)) {
      // value = (1+100*sin(CONST_PI*(centroid.x()+0.5))*sin(CONST_PI*(0.5+centroid.y())));
      value = 100.0;
      // std::cout << "Element " << handle << " with centroid matched criteria" << std::endl;
    }
    else
      value = 1.0;

    //value = m_solution_indicator[handle];
  }
  else {
    return_flag=false;
  }
  
  // EntityTopology type = element->get_element_type();

  // Look up conditionNumberQualityMetric for more details
    
  return return_flag;
}
/*
bool ElementSolIndQM::evaluate_with_gradient( PatchData& pd, 
                                            size_t handle, 
                                            double& value, 
                                            std::vector<size_t>& indices,
                                            std::vector<Vector3D>& gradient,
                                            MsqError& err )
{
  ElemSampleQM* qm = get_quality_metric();
  mHandles.clear();
  qm->get_element_evaluations( pd, handle, mHandles, err ); MSQ_ERRFALSE(err);

  bool valid = true;
  double tmpval;
  bool tmpvalid;
  unsigned count = 0;
  std::vector<size_t>::iterator h, i, j;

  value = -std::numeric_limits<double>::maximum();
  for (h = mHandles.begin(); h != mHandles.end(); ++h) {
    mIndices.clear();
    mGrad.clear();
    tmpvalid = qm->evaluate_with_gradient( pd, *h, tmpval, mIndices, mGrad, err );
    MSQ_ERRZERO(err);
    if (!tmpvalid) {
      valid = false;
    }
      // new value greater than previous max value
    else if (tmpval - value > 1e-6) {
      indices = mIndices;
      gradient = mGrad;
      count = 1;
      value = tmpval;
    }
      // new value equal to previous max value
    else if (tmpval - value >= -1e-6) {
      ++count;
      for (i = mIndices.begin(); i != mIndices.end(); ++i) {
        j = std::find( indices.begin(), indices.end(), *i );
        if (j == indices.end()) {
          indices.push_back( *i );
          gradient.push_back( mGrad[i - mIndices.begin()] );
        }
        else {
          gradient[j - indices.begin()] += mGrad[i - mIndices.begin()];
        }
      }
    }
  }
  if (count > 1) {
    const double inv_count = 1.0 / count;
    for (std::vector<Vector3D>::iterator g = gradient.begin(); g != gradient.end(); ++g)
      *g *= inv_count;
  }
    
  return valid;
}
*/

} // namespace Mesquite
