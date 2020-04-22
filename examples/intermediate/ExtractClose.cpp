/** \brief This test shows how to extract mesh from a model, based on distance.
 *
 * MOAB's It is needed to extract from large mesh files cells close to some point, where we suspect errors
 *  It would be useful for 9Gb input file that we cannot visualize, but we have overlapped elements.
 */

#include <iostream>
#include <cstdlib>
#include <cstdio>

#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/CartVect.hpp"

using namespace moab;
using namespace std;

#ifdef MOAB_HAVE_HDF5
string test_file_name = string(MESH_DIR) + string("/64bricks_512hex_256part.h5m");
#endif
int main(int argc, char **argv)
{
  ProgOptions opts;

  // vertex 25 in filt out  (-0.784768, 0.594952, -0.173698)
  double x=-0.784768, y=0.594952, z=-0.173698;

  string inputFile = test_file_name;
  opts.addOpt<string>("inFile,i", "Specify the input file name string ", &inputFile);

  string outFile = "out.h5m";
  opts.addOpt<string>("outFile,o", "Specify the output file name string ", &outFile);

  opts.addOpt<double>(string("xpos,x"), string("x position"), &x);
  opts.addOpt<double>(string("ypos,y"), string("y position"), &y);
  opts.addOpt<double>(string("zpos,z"), string("z position"), &z);

  double distance =0.01;
  opts.addOpt<double>(string("distance,d"), string("distance "), &distance);

  opts.parseCommandLine(argc, argv);

  // Instantiate
  Core mb;

  // Load the file
  ErrorCode rval = mb.load_file(inputFile.c_str());MB_CHK_SET_ERR(rval, "Error loading file");

  // Get all 2d elements in the file
  Range elems;
  rval = mb.get_entities_by_dimension(0, 2, elems);MB_CHK_SET_ERR(rval, "Error getting 2d elements");

  // create a meshset with close elements
  EntityHandle outSet;

      // create meshset
  rval = mb.create_meshset(MESHSET_SET, outSet);MB_CHK_ERR(rval);

  CartVect point(x,y,z);
  Range closeByCells;
  for (Range::iterator it=elems.begin(); it!=elems.end(); it++)
  {
    EntityHandle cell=*it;
    CartVect center;
    rval = mb.get_coords(&cell, 1, &(center[0])); MB_CHK_SET_ERR(rval, "Can't get cell center coords");
    double dist = (center-point).length();
    if (dist <= distance)
    {
      closeByCells.insert(cell);
    }
  }

  rval = mb.add_entities(outSet, closeByCells); MB_CHK_SET_ERR(rval, "Can't add to entity set");

  if (closeByCells.empty())
  {
    std::cout << " no close cells to the point [" << x << ", " << y << ", " << z <<"] at distance less than " << distance << "\n";
  }
  else
  {
    rval = mb.write_file(outFile.c_str(), 0, 0, &outSet, 1); MB_CHK_SET_ERR(rval, "Can't write file");
  }
  return 0;
}
