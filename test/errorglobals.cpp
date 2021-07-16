/*
 * =====================================================================================
 *
 *       Filename:  errorglobals.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  07/14/2021 10:09:01
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Vijay S. Mahadevan (vijaysm), mahadevan@anl.gov
 *        Company:  Argonne National Lab
 *
 * =====================================================================================
 */

/*
Demonstrate min/max angle quality metric for multiple face types
*/

#include "moab/Core.hpp"
#include "moab/Util.hpp"
#include "moab/ErrorHandler.hpp"

#include <iostream>

using namespace moab;
using namespace std;

#ifndef MESH_DIR
#define MESH_DIR "."
#endif

string test_file_name = "test_poly.h5m";

int main(int argc, char **argv)
{
  if (argc > 1) {
    // User has input a mesh file
    test_file_name = argv[1];
  }

  // Instantiate & load a mesh from a file
  Core* mb = new (std::nothrow) Core;
  if (NULL == mb)
    return 1;

  // read mesh first time and list all entities - should be no error
  ErrorCode rval = mb->load_mesh(test_file_name.c_str()); MB_CHK_SET_ERR_RET_VAL(rval, "Trouble loading mesh", rval);
  mb->list_entities(0, -1);
  
  Range duments;
  rval = mb->get_entities_by_dimension(0, 0, duments);
  double x, y, z;
  // just doing the following generates an error, which changes the global moab::lastError from ErrorHandler.cpp
  Util::normal(mb, *duments.begin(), x, y, z);


  // destroy the instance, re-create it, read the same file again, and list again
  delete mb;
  mb = new (std::nothrow) Core;

  rval = mb->load_mesh(test_file_name.c_str()); MB_CHK_SET_ERR_RET_VAL(rval, "Trouble loading mesh", rval);
  // second time, listing entities prints error
  mb->list_entities(0, -1);

  return 0;
}

