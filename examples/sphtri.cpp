/** @example sphtri       
 * Description: read a mesh on a sphere\n
 *
 *    sphtri  <end_points.dat> <triangle.dat> <bisected_points.dat>  <boundary_points.dat> <SaveLoopCounts>
 * (default values can run if users don't specify a mesh file)
 */


#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"

#include <iostream>
#include<fstream>


using namespace moab;
using namespace std;


string points_name = string("end_points.dat");
string triangles_name = string("triangles.dat");
string bis_name = string("bisected_points.dat");

int main(int argc, char **argv)
{
  // Get MOAB instance
  Interface* mb = new (std::nothrow) Core;
  if (NULL == mb)
    return 1;

  // Need option handling here for input filename
  if (argc > 1) {
    // User has input a mesh file for end points
    points_name = argv[1];
  }
  if (argc > 2) {
    // User has input a mesh file for triangle connectivity
    triangles_name = argv[2];
  }
  if (argc > 3) {
    // User has input a mesh file for bisected points
    bis_name = argv[3];
  }

  ifstream pointsFile(points_name.c_str());
  if(!pointsFile) {
    cout << endl << "Failed to open file " << points_name << endl;;
    return 1;
  }
  vector<double>  coords;
  while(!pointsFile.eof()) 
  {
    double co;
    pointsFile >> co;
    coords.push_back(co);
  }
  cout<<" number of points:" << coords.size()/3 << "\n";
  ifstream triFile(triangles_name.c_str());
  if(!triFile) {
    cout << endl << "Failed to open file " << triangles_name << endl;;
    return 1;
  }
  vector<int>  conns;
  while(!triFile.eof()) 
  {
    int co;
    triFile >> co;
    conns.push_back(co);
  }
  cout<<" number of triangles:" << conns.size()/3 << "\n";
  
  ReadUtilIface* readMeshIface;
  mb->query_interface(readMeshIface);
  Range verts;
  ErrorCode rval = mb->create_vertices(&coords[0], 
                                      coords.size()/3,
                                      verts ); MB_CHK_SET_ERR(rval, "do not create v");

  int numv=(int)coords.size()/3;
  EntityHandle* conn_array;
  EntityHandle handle = 0;
  rval = readMeshIface->get_element_connect( conns.size()/3 , 3, MBTRI,
                                              1,
                                              handle, conn_array); MB_CHK_SET_ERR(rval, "do not elems");
  
  for (size_t i=0; i<conns.size(); i++)
     conn_array[i] = conns[i]+1; // vertices for sure start at 1
   
  EntityHandle tri_set;
  rval = mb->create_meshset(MESHSET_SET, tri_set);MB_CHK_SET_ERR(rval, "Can't create mesh set");
  Range tris(handle, handle+ conns.size()/3-1);
  rval = mb->add_entities(tri_set, tris);MB_CHK_SET_ERR(rval, "Can't add triangles to set");
  
  if (argc>3)
  {
    ifstream bisFile(bis_name.c_str());
    if(!bisFile) {
      cout << endl << "Failed to open file " << bis_name << endl;;
      return 1;
    }
    coords.resize(0);
    while(!bisFile.eof())
    {
      double co;
      bisFile >> co;
      coords.push_back(co);
    }
    cout<<" number of bisected points:" << coords.size()/3 - numv << "\n";
    verts.clear();
    rval = mb->create_vertices(&coords[3*numv], // the first numv points are repeated, they are the end points
        coords.size()/3 - numv,
                                        verts ); MB_CHK_SET_ERR(rval, "do not create bis vertices");
    rval = mb->add_entities(tri_set, verts);MB_CHK_SET_ERR(rval, "Can't add bisection points to set");

  }
  
  mb-> write_file("tris.vtk", 0, "CREATE_ONE_NODE_CELLS;", &tri_set, 1);  // this option is for bisection points

  delete mb;

  return 0;
}
