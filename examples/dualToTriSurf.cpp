/** @example dualForMPAS
 * Description: read a delaunay triangular mesh on a sphere and convert to MPAS\n
 *  assume density_points.dat exists too, for density
 *    dualForMPAS  <end_points.dat> <triangle.dat> <offset>
 * (default values can run if users don't specify a mesh file)
 */


#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"
#include "math.h"
#include "moab/CartVect.hpp"

#include <iostream>
#include <fstream>

using namespace moab;
using namespace std;

string input_file = string("tris.h5m");

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
    input_file = argv[1];
  }

  cout<< " execute:  " << argv[0] << " " << input_file << "\n";

  ErrorCode rval = mb-> load_file(input_file.c_str()); MB_CHK_SET_ERR(rval, "Can't read input mesh\n");

  Range tris, verts;
  rval = mb->get_entities_by_type(0, MBTRI, tris);

  rval = mb->get_entities_by_type(0, MBVERTEX, verts);
  cout<<" number of triangles:" << tris.size()<< " number of vertices:" << verts.size()<<"\n";

  Tag gid;
  rval = mb->tag_get_handle("GLOBAL_ID", 1, MB_TYPE_INTEGER, gid); MB_CHK_SET_ERR(rval, "can't get global id tag");
  Tag densityTag;
  double densval = 0;
  rval = mb->tag_get_handle("Density", 1, MB_TYPE_DOUBLE, densityTag,
      MB_TAG_CREAT | MB_TAG_DENSE, &densval); MB_CHK_SET_ERR(rval, "can't get density tag");

  EntityHandle tri_set;
  rval = mb->create_meshset(MESHSET_SET, tri_set);MB_CHK_SET_ERR(rval, "Can't create mesh set");

  rval = mb->add_entities(tri_set, tris);MB_CHK_SET_ERR(rval, "Can't add triangles to set");

  // dual of the triangular mesh is the MPAS mesh
  // first find out all edges in the triangulation
  Range edges;
  rval = mb->get_adjacencies(tris,
      /*const int to_dimension*/ 1,
      /* const bool create_if_missing*/ true,
      /*Range &adj_entities*/ edges,
      /*const int operation_type = Interface::INTERSECT*/ Interface::UNION); MB_CHK_SET_ERR(rval, "Can't create edges");
  rval = mb->add_entities(tri_set, edges);MB_CHK_SET_ERR(rval, "Can't add edges to set");

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
  int global_id;
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
    rval = mb->tag_set_data(dualTag, &dual_vertex, 1, &triangle); MB_CHK_SET_ERR(rval, "Can't set dual tag triangle ");

    // also set the gid on the dual vertex, which will be the same as the triangle it came on (dual)

    rval = mb->tag_get_data(gid, &triangle, 1, &global_id); MB_CHK_SET_ERR(rval, "Can't get gid on triangle ");
    rval = mb->tag_set_data(gid, &dual_vertex, 1, &global_id); MB_CHK_SET_ERR(rval, "Can't set gid on dual vertex ");
  }

  Range dualEdges;
  for (Range::iterator eit=edges.begin(); eit!=edges.end(); eit++)
  {
    EntityHandle edge=*eit;
    // get triangles adjacent
    vector<EntityHandle> tris2;
    rval = mb->get_adjacencies(&edge, 1, 2, false, tris2); MB_CHK_SET_ERR(rval, "Can't get adjacent triangles");
    if ((int)tris2.size() !=2 )
      continue; // do not care about boundary edges
    vector<EntityHandle> dualVerts;
    dualVerts.resize(2); // 2 triangles adjacent
    rval = mb->tag_get_data(dualTag, &tris2[0], 2, &dualVerts[0]);MB_CHK_SET_ERR(rval, "Can't get dual verts");

    // create now an edge with these 2 vertices
    EntityHandle dualEdge ;
    rval = mb->create_element(MBEDGE, &dualVerts[0], 2, dualEdge);MB_CHK_SET_ERR(rval, "Can't create dual edge");

    rval = mb->tag_set_data(dualTag, &edge, 1, &dualEdge); MB_CHK_SET_ERR(rval, "Can't set dual tag edge ");
    rval = mb->tag_set_data(dualTag, &dualEdge, 1, &edge); MB_CHK_SET_ERR(rval, "Can't set dual tag edge ");
    dualEdges.insert(dualEdge);
  }

  // for each point in original triangulation, create the dual polygon

  Range dualPolygons;

  for (Range::iterator vit=verts.begin(); vit!=verts.end(); vit++)
  {
    EntityHandle vertex=*vit;
    // the issue is how to arrange the dual vertices in order
    rval = mb->tag_get_data(gid, &vertex, 1, &global_id); MB_CHK_SET_ERR(rval, "Can't get gid on vertex ");

    // start with first edge?
    vector<EntityHandle> edges2;
    rval = mb->get_adjacencies(&vertex, 1, 1, false, edges2); MB_CHK_SET_ERR(rval, "Can't get adjacent edges");
    vector<EntityHandle> dualEdges2;
    if (edges2.size()<3)
      continue; // less than 3 edges adjacent to a vertex; probably on boundary, skip anyway
    dualEdges2.resize(edges2.size());
    rval = mb->tag_get_data(dualTag, &edges2[0], (int)edges2.size(), &dualEdges2[0]); MB_CHK_SET_ERR(rval, "Can't get dual edges");

    Range dualEdgesRange;
    for (int j=0; j<(int)dualEdges2.size(); j++)
    {
      if (0!=dualEdges2[j])
        dualEdgesRange.insert(dualEdges2[j]);
    }

    if (dualEdgesRange.size() == dualEdges2.size())
    {
      // if some edges on boundary, no dual, skip for now zeros
      Range dualVerts;
      rval = mb->get_connectivity(dualEdgesRange, dualVerts); MB_CHK_SET_ERR(rval, "Can't get connectivity");
      if (dualVerts.size()<2)
        continue;
      // now form a loop with vertices in dualVerts
      vector<EntityHandle> polyv;
      polyv.resize(dualVerts.size());
      // first 2 vertices are from dualEdges[0]; orientation will follow later

      const EntityHandle * conn2;
      int nnodes;
      rval = mb->get_connectivity(dualEdgesRange[0], conn2, nnodes); MB_CHK_SET_ERR(rval, "Can't get conn of first edge ");
      polyv[0] = conn2[0];
      polyv[1] = conn2[1];
      // next, find the next nodes in the loop
      int currentIndex=1;
      dualEdgesRange.erase(dualEdgesRange[0]);
      while (currentIndex<(int)edges2.size()-1 )
      {
        EntityHandle dualv=polyv[currentIndex];
        // current edges adjacent
        Range dualEdgesFront;
        rval = mb->get_adjacencies(&dualv, 1, 1, false, dualEdgesFront);MB_CHK_SET_ERR(rval, "Can't get adj dual edges");

        Range front =intersect(dualEdgesRange, dualEdgesFront);
        if (front.size() != 1)
          break;
        rval = mb->get_connectivity(front[0], conn2, nnodes); MB_CHK_SET_ERR(rval, "Can't get conn of front edge ");
        EntityHandle nextV = conn2[0];
        dualEdgesRange.erase(front[0]);
        if (dualv==nextV)
          nextV = conn2[1];
        polyv[++currentIndex] = nextV;
      }
      if (currentIndex==(int)edges2.size()-1)
      {
      // finally create a polygon, if everything is fine
        EntityHandle polyg;
        rval = mb->create_element(MBPOLYGON, &polyv[0], polyv.size(), polyg ); MB_CHK_SET_ERR(rval, "Can't create polygon ");
        dualPolygons.insert(polyg);
        // also set the gid on the dual polyg, which should be the same as original vertex
        rval = mb->tag_set_data(gid, &polyg, 1, &global_id); MB_CHK_SET_ERR(rval, "Can't set gid on dual polygon ");
      }
    }
    else
    {
      // for each dual edge, add a triangle with original vertex
      for (int i=0; i<(int)dualEdgesRange.size(); i++)
      {
        EntityHandle dualEdge=dualEdgesRange[i];
        int nnodes;
        const EntityHandle * conn ;
        rval = mb->get_connectivity(dualEdge, conn, nnodes); MB_CHK_SET_ERR(rval, "Can't get conn on dual edge ");
        if (nnodes!=2)
          continue;
        EntityHandle conn3[3];
        conn3[0] = conn[0]; conn3[1]= conn[1];
        conn3[2] = vertex;
        EntityHandle partTri;
        rval = mb->create_element(MBTRI, conn3, 3, partTri); MB_CHK_SET_ERR(rval, "Can't create triangle ");
        rval = mb->tag_set_data(gid, &partTri, 1, &global_id); MB_CHK_SET_ERR(rval, "Can't set gid on dual triangle (part) ");
        dualPolygons.insert(partTri); // it is some triangles, not polygons
      }
    }

  }


  EntityHandle dual_set;
  rval = mb->create_meshset(MESHSET_SET, dual_set);MB_CHK_SET_ERR(rval, "Can't create dual set");

  rval = mb->add_entities(dual_set, dualPolygons); MB_CHK_SET_ERR(rval, "Can't add polygons to dual set");
  rval = mb->add_entities(dual_set, dualEdges); MB_CHK_SET_ERR(rval, "Can't add edges to dual set");

  // add to the dual set the original verts
  rval = mb->add_entities(dual_set, verts); MB_CHK_SET_ERR(rval, "Can't add original verts to dual set");

  // create polygons too, from dual verts
  std::stringstream dualFile;
  // erase 4 chars from input_file (extension, usually)
  if (input_file.size() > 4)
    input_file.resize(input_file.size()-4);
  dualFile<<input_file <<"_dual.h5m";

  cout << " writing " << dualFile.str() << "\n";
  rval = mb-> write_file(dualFile.str().c_str(), 0, 0, &dual_set, 1); MB_CHK_SET_ERR(rval, "Can't write dual mesh\n");

  delete mb;
  return 0;
}
