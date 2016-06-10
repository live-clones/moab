
/** \file ElementSolIndQM.hpp
 *  \brief 
 *  \author Vijay Mahadevan
 */

#ifndef MSQ_ELEMENT_SOL_IND_QM_HPP
#define MSQ_ELEMENT_SOL_IN_QM_HPP

#include "Mesquite.hpp"
#include "ElementQM.hpp"
#include "AveragingQM.hpp"
#include <vector>

namespace MBMesquite {

class ElemSampleQM;

class ElementSolIndQM : public ElementQM, public AveragingQM
{
public:

  MESQUITE_EXPORT ElementSolIndQM(std::vector<double>& solution_indicator);
  
  MESQUITE_EXPORT virtual ~ElementSolIndQM();

  MESQUITE_EXPORT virtual std::string get_name() const;

  MESQUITE_EXPORT virtual int get_negate_flag() const;

  MESQUITE_EXPORT virtual
  bool evaluate( PatchData& pd, 
                 size_t handle, 
                 double& value, 
                 MsqError& err );

private:

  const std::vector<double>& m_solution_indicator;

};


} // namespace Mesquite

#endif
