/** @example ExtrudeHelix.cpp
 * Description: read a 3d mesh in form  layers of hexes from top and bottom layers) \n
 *
 * To run: ./ExtrudeHelix [meshfile] [outfile] [topz] [bottomz] [layers_top] [layers_bot] \n
 */


#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"
#include "MBTagConventions.hpp"
#include <iostream>

using namespace moab;
using namespace std;

#ifndef MESH_DIR
#define MESH_DIR "."
#endif

// Note: change the file name below to test a trivial "No such file or directory" error
string test_file_name = string(MESH_DIR) + string("block100.exo");
string output = string("blo100.exo");
int layersTop = 30;
int layersBot = 30;
double topz = 75.;
double bottomz = 25.;

int main(int argc, char **argv)
{
  // Get MOAB instance
  Interface* mb = new (std::nothrow) Core;
  if (NULL == mb)
    return 1;

  // Need option handling here for input filename
  if (argc > 1) {
    // User has input a mesh file
    test_file_name = argv[1];
  }
  if (argc>2)
    output = argv[2];

  if (argc>3)
    topz   = atof(argv[3]);

  if (argc>4)
    bottomz  = atof(argv[4]);

  if (argc>5)
    layersTop = atoi(argv[5]);

  if (argc>6)
    layersBot  = atoi(argv[6]);

  std::cout << "Run: " << argv[0] << " " << test_file_name << " " << output<<
     " " << topz << " " << bottomz <<  " " << layersTop << " " << layersBot << " " << "\n";
  // Load the mesh from exo file
  ErrorCode rval = mb->load_mesh(test_file_name.c_str());MB_CHK_ERR(rval);

  // Get the tag handle for this tag name; tag should already exist (it was created during file read)
  Tag mtag, ntag;
  rval = mb->tag_get_handle("NEUMANN_SET", 1, MB_TYPE_INTEGER, ntag);MB_CHK_ERR(rval);
  rval = mb->tag_get_handle("MATERIAL_SET", 1, MB_TYPE_INTEGER, mtag);MB_CHK_ERR(rval);
 // Get all the sets having that neumann tag (with any value for that tag)
  Range sets;
  rval = mb->get_entities_by_type_and_tag(0, MBENTITYSET, &ntag, NULL, 1, sets);MB_CHK_ERR(rval);
  Range msets;
  rval = mb->get_entities_by_type_and_tag(0, MBENTITYSET, &mtag, NULL, 1, msets);MB_CHK_ERR(rval);

  vector<int> vals;
  vals.resize(sets.size());
  rval = mb-> tag_get_data(ntag, sets, &vals[0]); MB_CHK_ERR(rval);
  //assert(sets.size()>=2);
  EntityHandle setTop=sets[0] , setBot=sets[1];
  if (vals[0] == 2)
  {
    setTop=sets[1];
    setBot=sets[0];
  }
  Range quadsTop;
  rval = mb->get_entities_by_handle(setTop, quadsTop, true); MB_CHK_ERR(rval);
  Range quadsBot;
  rval = mb->get_entities_by_handle(setBot, quadsBot, true); MB_CHK_ERR(rval);

  
  // Get verts entities, by type
  Range verts;
  rval = mb->get_connectivity(quadsTop, verts);MB_CHK_ERR(rval);

  // Get faces, by dimension, so we stay generic to entity type
  cout << "Number of vertices on top is " << verts.size() << endl;
  cout << "Number of quads on top  is " << quadsTop.size() << endl;

  std::vector<double> coords;
  int nvPerLayer= (int)verts.size();
  coords.resize(3*nvPerLayer);

  rval = mb->get_coords(verts, &coords[0]); MB_CHK_ERR(rval);
  std::vector<double> zdelta;
  zdelta.resize(nvPerLayer);
  for (int i=0; i< nvPerLayer; i++)
  {
    zdelta[i] = (topz - coords[3*i+2])/layersTop; // each vertex can be different 
  }
  // create first vertices
  Range * newVerts = new Range [layersTop+1];
  newVerts[0] = verts; // just for convenience
  for (int ii=0; ii<layersTop; ii++)
  {
    for (int i=0;i<nvPerLayer; i++)
      coords[3*i+2] += zdelta[i];

    rval = mb->create_vertices(&coords[0], nvPerLayer, newVerts[ii+1]); MB_CHK_ERR(rval);
  }
  // for each quad, we will create layers hexas
  int nhexas = quadsTop.size()*layersTop;
  ReadUtilIface *read_iface;
  rval = mb->query_interface(read_iface);MB_CHK_SET_ERR(rval, "Error in query_interface");

  EntityHandle start_elem, *connect;
       // Create quads
  rval = read_iface->get_element_connect(nhexas, 8, MBHEX , 0, start_elem, connect);MB_CHK_SET_ERR(rval, "Error in get_element_connect");
  int nquads = (int)quadsTop.size();

  int indexConn=0;
  for (int j=0; j<nquads; j++)
  {
    EntityHandle quad = quadsTop[j];

    const EntityHandle *conn4 = NULL;
    int num_nodes;
    rval = mb->get_connectivity( quad, conn4, num_nodes) ; MB_CHK_ERR(rval);
    if (4!=num_nodes)
      MB_CHK_ERR(MB_FAILURE);
    int inx[4];
    for(int k=0; k<4; k++)
       inx[k] = verts.index(conn4[k]);

    for  (int ii=0; ii<layersTop; ii++)
    {
      connect[indexConn++] = newVerts[ii][inx[0]];
      connect[indexConn++] = newVerts[ii][inx[1]];
      connect[indexConn++] = newVerts[ii][inx[2]];
      connect[indexConn++] = newVerts[ii][inx[3]];
      connect[indexConn++] = newVerts[ii+1][inx[0]];
      connect[indexConn++] = newVerts[ii+1][inx[1]];
      connect[indexConn++] = newVerts[ii+1][inx[2]];
      connect[indexConn++] = newVerts[ii+1][inx[3]];
    }
  }
  Range newHexs(start_elem, start_elem+nhexas-1);
  mb->add_entities(msets[0], newHexs);
  
  verts.clear();
  rval = mb->get_connectivity(quadsBot, verts);MB_CHK_ERR(rval);

  // Get faces, by dimension, so we stay generic to entity type
  cout << "Number of vertices on bottom is " << verts.size() << endl;
  cout << "Number of quads on bottom  is " << quadsBot.size() << endl;

  nvPerLayer= (int)verts.size();
  coords.resize(3*nvPerLayer);

  rval = mb->get_coords(verts, &coords[0]); MB_CHK_ERR(rval);
  //std::vector<double> zdelta;
  zdelta.resize(nvPerLayer);
  for (int i=0; i< nvPerLayer; i++)
  {
    zdelta[i] = (coords[3*i+2]- bottomz)/layersBot; // each vertex can be different 
  }
  // create first vertices
  delete [] newVerts;
  newVerts = new Range [layersBot+1];
  newVerts[0] = verts; // just for convenience
  for (int ii=0; ii<layersBot; ii++)
  {
    for (int i=0;i<nvPerLayer; i++)
      coords[3*i+2] -= zdelta[i];

    rval = mb->create_vertices(&coords[0], nvPerLayer, newVerts[ii+1]); MB_CHK_ERR(rval);
  }
  // for each quad, we will create layers hexas
  nhexas = quadsBot.size()*layersBot;

  rval = read_iface->get_element_connect(nhexas, 8, MBHEX , 0, start_elem, connect);MB_CHK_SET_ERR(rval, "Error in get_element_connect");
  nquads = (int)quadsBot.size();

  indexConn=0;
  for (int j=0; j<nquads; j++)
  {
    EntityHandle quad = quadsBot[j];

    const EntityHandle *conn4 = NULL;
    int num_nodes;
    rval = mb->get_connectivity( quad, conn4, num_nodes) ; MB_CHK_ERR(rval);
    if (4!=num_nodes)
      MB_CHK_ERR(MB_FAILURE);
    int inx[4];
    for(int k=0; k<4; k++)
       inx[k] = verts.index(conn4[k]);

    for  (int ii=0; ii<layersBot; ii++)
    {
      connect[indexConn++] = newVerts[ii][inx[0]];
      connect[indexConn++] = newVerts[ii][inx[1]];
      connect[indexConn++] = newVerts[ii][inx[2]];
      connect[indexConn++] = newVerts[ii][inx[3]];
      connect[indexConn++] = newVerts[ii+1][inx[0]];
      connect[indexConn++] = newVerts[ii+1][inx[1]];
      connect[indexConn++] = newVerts[ii+1][inx[2]];
      connect[indexConn++] = newVerts[ii+1][inx[3]];
    }
  }
  Range newHexsBot(start_elem, start_elem+nhexas-1);
  mb->add_entities(msets[0], newHexsBot);
  rval = mb->write_file(output.c_str(), 0, 0, msets);MB_CHK_ERR(rval);

  delete mb;

  return 0;
}
