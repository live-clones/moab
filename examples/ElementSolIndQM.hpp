/** \file ElementSolIndQM.hpp
 *  \brief 
 *  \author Vijay Mahadevan
 */

#ifndef MSQ_ELEMENT_SOL_IND_QM_HPP
#define MSQ_ELEMENT_SOL_IN_QM_HPP

#include "Mesquite.hpp"
#include "VertexQM.hpp"
#include "AveragingQM.hpp"

namespace MBMesquite
{
     /*! \class ElementSolIndQM
       \brief Computes the local size metric for a given vertex.
       
        ElementSolIndQM is a vertex based metric which computes
        the corner volume (or area) for the element corners attached
        to a given element.  Then these volumes (or areas) are averaged
        together.  The default averaging method is QualityMetric::RMS.
     */
   class ElementSolIndQM : public VertexQM, public AveragingQM
   {
  public:
        //Default constructor. 
      ElementSolIndQM(std::vector<double>& solution_indicator) : AveragingQM(SUM_SQUARED), m_solution_indicator(solution_indicator) {}

       // virtual destructor ensures use of polymorphism during destruction
     virtual ~ElementSolIndQM() {}
     
     virtual std::string get_name() const;
     
     virtual int get_negate_flag() const;
     
     virtual
     bool evaluate( PatchData& pd, 
                    size_t handle, 
                    double& value, 
                    MsqError& err );
     
     virtual
     bool evaluate_with_indices( PatchData& pd,
                                 size_t handle,
                                 double& value,
                                 std::vector<size_t>& indices,
                                 MsqError& err );

  private:

    const std::vector<double>& m_solution_indicator;

  };

} //namespace

#endif
