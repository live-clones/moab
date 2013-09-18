/*
 * proj1.cpp
 *
 *  project on a sphere of radius R
 */

#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include <iostream>
#include <math.h>

#include "CslamUtils.hpp"
#include <assert.h>
using namespace moab;

double radius = 1.;// in m:  6371220.

int main(int argc, char **argv)
{

  std::string extra_read_opts;
  // read a file and project on a sphere

  if (argc < 3)
    return 1;

  int index = 1;
  char * input_mesh1 = argv[1];
  char * output = argv[2];
  while (index < argc)
  {
    if (!strcmp(argv[index], "-R")) // this is for geometry tolerance
    {
      radius = atof(argv[++index]);
    }

    index++;
  }

  Core moab;
  Interface & mb = moab;

  ErrorCode rval;

  rval = mb.load_mesh(input_mesh1);

  std::cout  << " -R " << radius << " input: " << input_mesh1 <<
      "  output: " << output << "\n";

  Range verts;
  rval = mb.get_entities_by_dimension(0, 0, verts);
  if (MB_SUCCESS != rval)
    return 1;

  double *x_ptr, *y_ptr, *z_ptr;
  int count;
  rval = mb.coords_iterate(verts.begin(), verts.end(), x_ptr, y_ptr, z_ptr, count);
  if (MB_SUCCESS != rval)
      return 1;
  assert(count == (int) verts.size()); // should end up with just one contiguous chunk of vertices

  for (int v = 0; v < count; v++) {
     //EntityHandle v = verts[v];
     CartVect pos( x_ptr[v], y_ptr[v] , z_ptr[v]);
     pos = pos/pos.length();
     pos = radius*pos;
     x_ptr[v] = pos[0];
     y_ptr[v] = pos[1];
     z_ptr[v] = pos[2];
  }

  Range edges;
  rval = mb.get_entities_by_dimension(0, 1, edges);
  if (MB_SUCCESS != rval)
    return 1;
  mb.delete_entities(edges);
  mb.write_file(output);

  // remove all edges


  return 0;
}
