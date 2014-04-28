/*
 * NCWriteHOMME.hpp
 *
 *  nc write helper for HOMME type data (CAM)
 *  Created on: April 9, 2014
 *
 */

#ifndef NCWRITEHOMME_HPP_
#define NCWRITEHOMME_HPP_

#include "NCWriteHelper.hpp"

namespace moab {

class NCWriteHOMME: public UcdNCWriteHelper
{
public:
  NCWriteHOMME(WriteNC* writeNC, int fileId, const FileOptions& opts, EntityHandle fileSet)
: UcdNCWriteHelper(writeNC, fileId, opts, fileSet) {}

  virtual ~NCWriteHOMME();

private:
  //! Implementation of NCWriteHelper::collect_mesh_info()
  virtual ErrorCode collect_mesh_info();

  //! Collect data for specified variables
  virtual ErrorCode collect_variable_data(std::vector<std::string>& var_names);

  //! Implementation of NCWriteHelper::write_values()
  virtual ErrorCode write_values(std::vector<std::string>& var_names);
};

} // namespace moab

#endif
