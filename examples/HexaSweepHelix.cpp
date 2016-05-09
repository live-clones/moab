/*
 * create a helix mesh;
 *  input file : surface side mesh 
 * 
 *  HexaSweepHelix -i side1.exo 
 */
#include "moab/Core.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/ParallelComm.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/CartVect.hpp"
#include "MBTagConventions.hpp"

#include <time.h>
#include <iostream>
#include <vector>

using namespace moab;
using namespace std;

int main(int argc, char **argv)
{
  double pitch = 15.0; // -1.875 total pitch is from one end to the other ; what is the unit? inch!
                        // will consider the pitch for  a full revolution
                        // for 45 degrees it will be 15/8, -1.875
  double angle = 45.;
  int orglayers = 128; // for 45 degrees, there are 128 orglayers; first we will see if we can generate 
//  qq starting from the (skin) boundary edges. will start from y=~0 plane (xz), 45 degrees; see if we 
// can get the other end
  ProgOptions opts;

  string inputf =  "side1.exo";
  opts.addOpt<string>("input,i", "Specify the input file name string (default input.cub)", &inputf);
  opts.addOpt<double>(string("pitch,p"), string("Total pitch (default=-15)"), &pitch);
  opts.addOpt<double>(string("angle,a"),  string("Total angle (default=-45.)"), &angle);
  opts.addOpt<int>(string("layers,l"),  string("num layers (default=128)"), &orglayers);

  opts.parseCommandLine(argc, argv);

  double angleRadInc = angle/180.*M_PI/orglayers; // 
  double deltazInc =  pitch/orglayers* (angle/360.) ;
  int numLayersSupport = 5;
  double supportAngle = 1.; // in degrees
  double angleIncSupport = supportAngle * M_PI/180./numLayersSupport; // so total support angle will be 1 degree 
  double deltazIncSupport = pitch/2/M_PI*angleIncSupport; // 

  Core mb;

  ErrorCode rval = mb.load_file(inputf.c_str()); MB_CHK_ERR(rval);

  ReadUtilIface* iface;
  rval = mb.query_interface(iface);MB_CHK_SET_ERR(rval, "Can't get reader interface");

  // get first the edges, by skinning the quads
  Range orgQuads; 
  rval = mb.get_entities_by_dimension(0, 2, orgQuads); MB_CHK_SET_ERR(rval, "Can't get quads");
  cout <<" number of quads: " << orgQuads.size() << endl;

  Range verts;
  rval = mb.get_connectivity(orgQuads, verts); MB_CHK_SET_ERR(rval, "Can't get vertices on quads");
  cout <<" number of vertices: " << verts.size() << endl;

  vector<double> coords;
  coords.resize(verts.size()*3);
  rval = mb.get_coords(verts, &coords[0]); MB_CHK_SET_ERR(rval, "Can't get coords of vertices");
// separate vertices based on coordinate y; those with y ~= 0 will be our base vertices
  

  int layers = 5  +  6; // 5 layers for support 1, 6 for support 2

  int num_newNodes = (layers) * (int)verts.size();

  vector<double*> arrays;
  EntityHandle startv;
  rval = iface->get_node_coords(3, num_newNodes, 0, startv, arrays);MB_CHK_SET_ERR(rval, "Can't get node coords");

  std::vector<double>  inico(3*verts.size());
  rval = mb.get_coords(verts, &inico[0]); MB_CHK_ERR(rval);

  // for each layer, compute an angle
  int nv = (int)verts.size();
  for (int i=0; i<layers; i++)
  {
    double arad=angleIncSupport*(i+1); // for last layer (128?), it will be angle * 180/M_PI
    double deltaz = deltazIncSupport* (i+1); // for last layer it will be total pitch 
    if (i>=5)
    {
      // add 128
      arad = angleRadInc*(128) + angleIncSupport*(i);
      deltaz = deltazInc* (128) + deltazIncSupport* (i);
    }
    for (int j=0; j<nv; j++)
    {
      double x=inico[3*j];

      double z=inico[3*j+2]; // basically original y is 0!!!
      arrays[0][i*nv+j] = x*cos(arad);
      arrays[1][i*nv+j] = x*sin(arad); // rotate counter-clockwise here, so y becomes positive here for positive x (angle is now positive)
      arrays[2][i*nv+j] = z+deltaz;
    }
  }

  int num_hexas = (layers-1) * (int)orgQuads.size(); // will generate layers of hexa for each quad in xz plane

  EntityHandle starte; // Connectivity
  EntityHandle* conn;
  //int num_v_per_elem = 8;

  rval = iface->get_element_connect(num_hexas, 8, MBHEX, 0, starte, conn);MB_CHK_SET_ERR(rval, "Can't get element connectivity for new hexas");

  int ixh=0; // index in conn array
  for (Range::iterator eit=orgQuads.begin(); eit!=orgQuads.end(); eit++)
  {
    EntityHandle quad=*eit;
    const EntityHandle * conn4 =NULL;
    int nno;
    rval = mb.get_connectivity(quad, conn4, nno); MB_CHK_ERR(rval);
    if (nno!=4 )continue;// should error out actually 

    int indices[4];
    for (int k=0; k<4; k++)
      indices[k] = verts.index(conn4[k]);

    // first layer hexa will use connectivity vertices conn4; next will use the indices
    // 
    conn[ixh++]= conn4[0];
    conn[ixh++]= conn4[1];
    conn[ixh++]= conn4[2];
    conn[ixh++]= conn4[3];
    conn[ixh++]= startv  + indices[0];
    conn[ixh++]= startv  + indices[1];
    conn[ixh++]= startv  + indices[2];
    conn[ixh++]= startv  + indices[3];
  }
  
  for (int i=2; i<=5; i++) // first layers, 2, 3, 4 and 5
  {
    std::cout <<" Layer : " << i << endl;
    for (Range::iterator eit=orgQuads.begin(); eit!=orgQuads.end(); eit++)
    {
      EntityHandle quad=*eit;
      const EntityHandle * conn4 =NULL;
      int nno;
      rval = mb.get_connectivity(quad, conn4, nno); MB_CHK_ERR(rval);
      if (nno!=4 )continue;// should error out actually 

      int indices[4];
      for (int k=0; k<4; k++)
        indices[k] = verts.index(conn4[k]);
     
      conn[ixh++]= startv + indices[0]+(i-2)*nv;
      conn[ixh++]= startv + indices[1]+(i-2)*nv;
      conn[ixh++]= startv + indices[2]+(i-2)*nv;
      conn[ixh++]= startv + indices[3]+(i-2)*nv;
      conn[ixh++]= startv + indices[0]+(i-1)*nv;
      conn[ixh++]= startv + indices[1]+(i-1)*nv;
      conn[ixh++]= startv + indices[2]+(i-1)*nv;
      conn[ixh++]= startv + indices[3]+(i-1)*nv;
    }
  }

  Tag dir_tag;
  rval = mb.tag_get_handle(MATERIAL_SET_TAG_NAME, dir_tag); MB_CHK_SET_ERR(rval, "Can't get dirichlet tag");

  Tag ptag;
  int dum_id = -1;
  rval = mb.tag_get_handle("PARALLEL_PARTITION", 1, MB_TYPE_INTEGER,
                             ptag, MB_TAG_CREAT | MB_TAG_SPARSE, &dum_id);MB_CHK_SET_ERR(rval, "Can't get parallel partition tag");

  Range sets;
  EntityHandle  setSupp1,  setSupp2;

  Range newHexSupport1(starte, starte+(int)orgQuads.size()*5-1);

  rval = mb.create_meshset(MESHSET_SET, setSupp1);MB_CHK_ERR(rval);
  rval = mb.add_entities(setSupp1, newHexSupport1);
  sets.insert(setSupp1);

  int valDirTag = 100;
  rval = mb.tag_set_data(dir_tag, &setSupp1, 1, &valDirTag);MB_CHK_ERR(rval);

  dum_id = 10;
  rval = mb.tag_set_data(ptag, &setSupp1, 1, &dum_id);MB_CHK_ERR(rval);
  

  for (int i=6; i<=10; i++) // layers,  6, 7, .. 11 
  {
    std::cout <<" Layer : " << i << endl;
    for (Range::iterator eit=orgQuads.begin(); eit!=orgQuads.end(); eit++)
    {
      EntityHandle quad=*eit;
      const EntityHandle * conn4 =NULL;
      int nno;
      rval = mb.get_connectivity(quad, conn4, nno); MB_CHK_ERR(rval);
      if (nno!=4 )continue;// should error out actually 

      int indices[4];
      for (int k=0; k<4; k++)
        indices[k] = verts.index(conn4[k]);
     
      conn[ixh++]= startv + indices[0]+(i-1)*nv;
      conn[ixh++]= startv + indices[1]+(i-1)*nv;
      conn[ixh++]= startv + indices[2]+(i-1)*nv;
      conn[ixh++]= startv + indices[3]+(i-1)*nv;
      conn[ixh++]= startv + indices[0]+(i)*nv;
      conn[ixh++]= startv + indices[1]+(i)*nv;
      conn[ixh++]= startv + indices[2]+(i)*nv;
      conn[ixh++]= startv + indices[3]+(i)*nv;
    }
  }
  Range newHexSupport2( starte+(int)orgQuads.size()*5, starte+(int)orgQuads.size()*10 -1);

  rval = mb.create_meshset(MESHSET_SET, setSupp2);MB_CHK_ERR(rval);
  rval = mb.add_entities(setSupp2, newHexSupport2);
  sets.insert(setSupp2);

  valDirTag = 110;
  rval = mb.tag_set_data(dir_tag, &setSupp2, 1, &valDirTag);MB_CHK_ERR(rval);

  dum_id = 11;
  rval = mb.tag_set_data(ptag, &setSupp2, 1, &dum_id);MB_CHK_ERR(rval);
  

  rval = mb.write_file("hexes1.h5m", 0, 0, sets);MB_CHK_ERR(rval);
  return 0;
}
