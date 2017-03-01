/** @example extractContinent
 * Description: read a mesh on a sphere and boundaries of continents, and extract a specified continent\n
 *  boundaries of continents exist as 2 files, boundary points and loops
 *    extractContinent  <input.h5m> <loop_number>
 *  default values: dual.h5m  <2> : loop number 2 is for Americas.
 */


#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"
#include "math.h"
#include "moab/CartVect.hpp"

#include <iostream>
#include <fstream>

using namespace moab;
using namespace std;

string bd_name = string("boundary_points.dat");
string loops = string("SaveLoopCounts");
string input_file = string("dual.h5m");

double getLat(CartVect p)  {
  p.normalize();
  return asin(p[2]);
}
double getLon(CartVect p)  {
  p.normalize();
  double lon;

  lon = atan2(p[1],p[0]);

  if(lon < 0){
    return 2.0 * M_PI + lon;
  } else {
    return lon;
  }
}


bool interior_point(vector<double> & coords, int& startLoop, int & endLoop, double lat, double lon)
{
  // compute oriented angle between point and edges in the loop
  // if angle ~ +/-2Pi, it is in the interior
  double totalAngle = 0;
  for (int i=startLoop; i<=endLoop; i++)
  {
    double x1=coords[2*i], y1 = coords[2*i+1];
    int i2 = i+1;
    if (i==endLoop)
      i2 = startLoop;
    double x2 = coords[2*i2], y2= coords[2*i2+1];
    CartVect P1( x1-lat, y1-lon, 0.);
    CartVect P2 (x2-lat, y2-lon, 0.);
    CartVect cross = P1*P2;
    double angle1 = angle(P1, P2);
    if (cross[2]<0)
      angle1=-angle1;
    totalAngle += angle1;
  }
  if ( fabs(totalAngle+2*M_PI)<1.e-3 || fabs(totalAngle-2*M_PI)<1.e-3)
    return true;
  return false;
}

int main(int argc, char **argv)
{
  // Get MOAB instance
  Interface* mb = new (std::nothrow) Core;
  if (NULL == mb)
    return 1;


  if (argc >1) {
    // User has input file for mesh file
    input_file = argv[1];
  }
  int loop_index = 2;
  if (argc > 2) {
    // User has input file for loops for continents/islands
    loop_index = atoi(argv[2]);
  }
  bool NorthAmerica = true;
  if (2!=loop_index)
    NorthAmerica = false;

  cout << "input file : " << input_file << " loop index:" << loop_index <<"\n";
  std::vector<double> coords;
  ifstream bdFile(bd_name.c_str());
  if(!bdFile) {
    cout << endl << "Failed to open file " << bd_name << endl;;
    return 1;
  }
  coords.resize(0);
  while(!bdFile.eof())
  {
    double co;
    bdFile >> co;
    coords.push_back(co);
  }
  cout<<" number of boundary points:" << coords.size()/3 << "\n";
  /// get loops for boundaries
  vector<int> loopsindx;
  ifstream loopFile(loops.c_str());
  while(!loopFile.eof())
  {
    int indx;
    loopFile >> indx;
    loopsindx.push_back(indx);
  }
  cout << "number of loops: " << loopsindx.size()/2 << "\n";
  // convert loops to lat-lon loops
  // it will be easier to determine in-out relations for points in 2d plane
  // careful with the international time zone; ( + - Pi longitude )
  vector<double> coords2d;
  coords2d.resize(coords.size()/3*2);
  for (int i=0; i<(int)coords.size()/3; i++)
  {
    CartVect p(coords[3*i], coords[3*i+1], coords[3*i+2]);
    coords2d[2*i]=getLat(p);
    coords2d[2*i+1]=getLon(p);
  }

  ErrorCode rval = mb-> load_file(input_file.c_str());  MB_CHK_SET_ERR(rval, "Can't load file");
  // look at the center of element, and see if it is inside the loop

  Range cells;
  rval = mb->get_entities_by_dimension(0, 2, cells);MB_CHK_SET_ERR(rval, "Can't get cells");

  cout <<"number of cells: " << cells.size() << "\n";

  EntityHandle islandSet;
  rval = mb->create_meshset(MESHSET_SET, islandSet);MB_CHK_SET_ERR(rval, "Can't create island set");

  int startLoop=loops[2*loop_index-2];
  int endLoop = loops[2*loop_index-1];
  // north america
  if (NorthAmerica)
  {
    startLoop=1421;
    endLoop=2082;
  }
  vector<EntityHandle> interiorCells;
  for (Range::iterator cit = cells.begin(); cit!=cells.end(); cit++)
  {
    EntityHandle cell= *cit;
    // see if it is in the interior of the loop
    CartVect center;
    rval = mb->get_coords(&cell, 1, &(center[0])); MB_CHK_SET_ERR(rval, "Can't get cell center coords");
    double lat=getLat(center), lon = getLon(center);
    // if NA, use some boxes too, for lat/lon
    if (NorthAmerica && (lat < 0.15 || lon < M_PI || lon > 2*M_PI - 0.5) )
      continue;

    if ( interior_point(coords2d, startLoop, endLoop, lat, lon))
      interiorCells.push_back(cell);
  }

  rval = mb->add_entities(islandSet, &interiorCells[0], interiorCells.size()); MB_CHK_SET_ERR(rval, "Can't add entities to set");

  std::stringstream islandFile;
  // erase 4 chars from input_file (extension, usually)
  if (input_file.size() > 4)
    input_file.resize(input_file.size()-4);
  islandFile<<input_file <<"_Island_"<< loop_index<<".h5m";

  rval = mb-> write_file(islandFile.str().c_str(), 0, 0, &islandSet, 1); MB_CHK_SET_ERR(rval, "Can't write island file");
  return 0;
}