/** @example sphtri       
 * Description: read a mesh on a sphere\n
 *
 *    sphtri  <end_points.dat> <triangle.dat> <bisected_points.dat>  <boundary_points.dat> <SaveLoopCounts>
 * (default values can run if users don't specify a mesh file)
 */


#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"
#include "math.h"
#include "moab/CartVect.hpp"

#include <iostream>
#include<fstream>

using namespace moab;
using namespace std;


string points_name = string("end_points.dat");
string triangles_name = string("triangles.dat");

void circumcenter(const CartVect &A,const CartVect &B,const CartVect &C, CartVect &cent){/*{{{*/
  double a, b, c;
  double pbc, apc, abp;

  a = (B-C).length_squared();
  b = (C-A).length_squared();
  c = (A-B).length_squared();

  pbc = a*(-a + b + c);
  apc = b*( a - b + c);
  abp = c*( a + b - c);

  cent = (pbc*A + apc*B + abp*C)/(pbc + apc + abp);

}/*}}}*/

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
   
  EntityHandle tri_set;
  rval = mb->create_meshset(MESHSET_SET, tri_set);MB_CHK_SET_ERR(rval, "Can't create mesh set");
  Range tris(handle, handle+ conns.size()/3-1);
  rval = mb->add_entities(tri_set, tris);MB_CHK_SET_ERR(rval, "Can't add triangles to set");
  
  mb-> write_file("tris.vtk", 0, 0, &tri_set, 1);

  // dual of the triangular mesh is the MPAS mesh
  // first find out all edges in the triangulation
  Range edges;
  rval = mb->get_adjacencies(tris,
      /*const int to_dimension*/ 1,
      /* const bool create_if_missing*/ true,
      /*Range &adj_entities*/ edges,
      /*const int operation_type = Interface::INTERSECT*/ Interface::UNION); MB_CHK_SET_ERR(rval, "Can't create edges");
  rval = mb->add_entities(tri_set, edges);MB_CHK_SET_ERR(rval, "Can't add edges to set");

  mb-> write_file("tris.vtk", 0, 0, &tri_set, 1);
  // create a dual polygon for each vertex, and a dual edge for each edge
  // create first a dual vertex for each triangle, at circumcenter
  /*tag_get_handle( const char* name,
                                      int size,
                                      DataType type,
                                      Tag& tag_handle,
                                      unsigned flags = 0,
                                      const void* default_value = 0,
                                      bool* created = 0 )*/
  Tag dualTag;
  EntityHandle defa=0;
  rval = mb->tag_get_handle("__dual", 1, MB_TYPE_HANDLE, dualTag, MB_TAG_CREAT|MB_TAG_DENSE, &defa); MB_CHK_SET_ERR(rval, "Can't create def set");

  Range dualVertices;
  // first determine circumcenters, and normalize them
  for (Range::iterator tit=tris.begin(); tit!=tris.end(); tit++)
  {
    EntityHandle triangle=*tit;
    int nnodes;
    const EntityHandle * conn;
    rval = mb->get_connectivity(triangle, conn, nnodes); MB_CHK_SET_ERR(rval, "Can't get connectivity");
    if (nnodes!=3)
      MB_CHK_SET_ERR(MB_FAILURE, "Can't get triangle conn");
    CartVect posi[3];
    rval = mb->get_coords(conn, 3, &(posi[0][0]) );  MB_CHK_SET_ERR(rval, "Can't get positions ");
    CartVect center;
    circumcenter(posi[0], posi[1], posi[2], center);
    center.normalize(); // on the sphere always
    EntityHandle dual_vertex;
    double coordsv[3] ;
    for (int i=0; i<3; i++)
      coordsv[i] = center[i];
    rval = mb->create_vertex( coordsv, dual_vertex); MB_CHK_SET_ERR(rval, "Can't create dual vertex ");

    dualVertices.insert(dual_vertex);

    rval = mb->tag_set_data(dualTag, &triangle, 1, &dual_vertex); MB_CHK_SET_ERR(rval, "Can't set dual tag vertex ");
  }

  Range dualEdges;
  for (Range::iterator eit=edges.begin(); eit!=edges.end(); eit++)
  {
    EntityHandle edge=*eit;
    // get triangles adjacent
    vector<EntityHandle> tris2;
    rval = mb->get_adjacencies(&edge, 1, 2, false, tris2); MB_CHK_SET_ERR(rval, "Can't get adjacent triangles");
    if ((int)tris2.size() !=2 )
      MB_CHK_SET_ERR(MB_FAILURE, "non manifold sphere ");
    vector<EntityHandle> dualVerts;
    dualVerts.resize(2); // 2 triangles adjacent
    rval = mb->tag_get_data(dualTag, &tris2[0], 2, &dualVerts[0]);MB_CHK_SET_ERR(rval, "Can't get dual verts");

    // create now an edge with these 2 vertices
    EntityHandle dualEdge ;
    rval = mb->create_element(MBEDGE, &dualVerts[0], 2, dualEdge);MB_CHK_SET_ERR(rval, "Can't create dual edge");

    dualEdges.insert(dualEdge);
  }

  EntityHandle dual_set;
  rval = mb->create_meshset(MESHSET_SET, dual_set);MB_CHK_SET_ERR(rval, "Can't create dual set");

  rval = mb->add_entities(dual_set, dualEdges); MB_CHK_SET_ERR(rval, "Can't add edges to dual set");

  // add to the dual set the original verts
  rval = mb->add_entities(dual_set, verts); MB_CHK_SET_ERR(rval, "Can't add original verts to dual set");

  // create polygons too, from dual verts
  mb-> write_file("dual.vtk", 0, "CREATE_ONE_NODE_CELLS;", &dual_set, 1);

  delete mb;

  return 0;
}
