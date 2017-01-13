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
string bd_name = string("boundary_points.dat");
string loops = string("SaveLoopCounts");

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
  if (argc > 4) {
    // User has input file for boundary describing continents
    bd_name = argv[4];
  }
  if (argc > 5) {
    // User has input file for loops for continents/islands
    loops = argv[5];
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

  EntityHandle* conn_array;
  EntityHandle handle = 0;
  rval = readMeshIface->get_element_connect( conns.size()/3 , 3, MBTRI,
                                              1,
                                              handle, conn_array); MB_CHK_SET_ERR(rval, "do not elems");
  
  for (size_t i=0; i<conns.size(); i++)
     conn_array[i] = conns[i]+1; // vertices for sure start at 1
   
  
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
  cout<<" number of bisected points:" << coords.size()/3 << "\n";
  verts.clear();
  rval = mb->create_vertices(&coords[0], 
                                      coords.size()/3,
                                      verts ); MB_CHK_SET_ERR(rval, "do not create bis vertices");
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
  
  verts.clear();
  rval = mb->create_vertices(&coords[0], 
                                      coords.size()/3,
                                      verts ); MB_CHK_SET_ERR(rval, "do not create boundary vertices");

  
  int num_edges=0;
  for (size_t i=0; i<loopsindx.size()/2; i++)
  {
    num_edges += ( loopsindx[2*i+1] - loopsindx[2*i] +1);
  }
  rval = readMeshIface->get_element_connect( num_edges , 2, MBEDGE,
                                              1,
                                              handle, conn_array); MB_CHK_SET_ERR(rval, "do not elems");
  int i1 = 0; // index in edge connectivity
  for (size_t i=0; i<loopsindx.size()/2; i++)
  {
    for (int j=loopsindx[2*i]; j<=loopsindx[2*i+1]; j++)
    {
      int  j2;
      j2 = j;
      if (j == loopsindx[2*i+1]) 
        j2 = loopsindx[2*i]-1; // close the loop
      conn_array[2*i1] = verts[j-1]; // vertices for sure start at 1
      conn_array[2*i1+1] = verts[j2]; // vertices for sure start at 1
      i1++;
    }
    // the last vertex is first one in the loop
  }
 
  mb-> write_file("out.vtk", 0, "CREATE_ONE_NODE_CELLS;");  // this is for bisection points
  delete mb;

  return 0;
}
