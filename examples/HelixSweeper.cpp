/*
 * create a helix mesh;
 *  input file : surface mesh tube
 * 
 *  HelixSweep -i qq.exo -p 10 -a 360
 */
#include "moab/Core.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/ParallelComm.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/CartVect.hpp"
#include "moab/Skinner.hpp"
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

  string inputf =  "qq.exo";
  opts.addOpt<string>("input,i", "Specify the input file name string (default input.cub)", &inputf);
  opts.addOpt<double>(string("pitch,p"), string("Total pitch (default=-15)"), &pitch);
  opts.addOpt<double>(string("angle,a"),  string("Total angle (default=-45.)"), &angle);
  opts.addOpt<int>(string("layers,l"),  string("num layers (default=128)"), &orglayers);

  opts.parseCommandLine(argc, argv);

  double angleRadInc = angle/180.*M_PI/orglayers; // 
  double deltazInc =  pitch/orglayers* (angle/360.) ;
  Core mb;

  ErrorCode rval = mb.load_file(inputf.c_str()); MB_CHK_ERR(rval);

  ReadUtilIface* iface;
  rval = mb.query_interface(iface);MB_CHK_SET_ERR(rval, "Can't get reader interface");

  // get first the edges, by skinning the quads
  Range orgQuads; 
  rval = mb.get_entities_by_dimension(0, 2, orgQuads); MB_CHK_SET_ERR(rval, "Can't get quads");
  cout <<" number of quads: " << orgQuads.size() << endl;
  Skinner    skinner(&mb);
  Range edges; 
  rval = skinner.find_skin(0, orgQuads, 1, edges,true, true); MB_CHK_SET_ERR(rval, "Can't get skinned edges");
  cout <<" number of edges on skin: " << edges.size() << endl;

  Range verts;
  rval = mb.get_connectivity(edges, verts); MB_CHK_SET_ERR(rval, "Can't get vertices on skinned edges");
  cout <<" number of vertices on skin: " << verts.size() << endl;

  vector<double> coords;
  coords.resize(verts.size()*3);
  rval = mb.get_coords(verts, &coords[0]); MB_CHK_SET_ERR(rval, "Can't get coords of vertices");
// separate vertices based on coordinate y; those with y ~= 0 will be our base vertices
  
  Range xzverts, otherverts;
  Range::iterator vit = verts.begin();
  for (size_t i = 0; i< verts.size(); i++, vit++) 
  {
     EntityHandle v = *vit;
     if (fabs(coords[3*i+1]) < 0.1) // it means it is in the plane xz
        xzverts.insert(v);
     else
        otherverts.insert(v);
  }
  
  cout <<" number of vertices on xz plane: " << xzverts.size() << " range psize:" << xzverts.psize()<< endl;
  //xzverts.print(0); // show the range 
  cout <<" number of vertices on 45 degree plane: " << otherverts.size() << " range psize:" << otherverts.psize()<< endl;
  //otherverts.print(0); // show the range 
  Range xzedges;
   
  
  rval = mb.get_adjacencies(xzverts, 1, false, xzedges, Interface::UNION);
   cout <<" number of edges on xz plane: " << xzedges.size() << " psize: " << xzedges.psize() << endl;

  int layers = 3 + orglayers + 3 + 20; // 3 layers for support 1, 3 for support 2, 20 free 

  int num_newNodes = (layers) * (int)xzverts.size();

  vector<double*> arrays;
  EntityHandle startv;
  rval = iface->get_node_coords(3, num_newNodes, 0, startv, arrays);MB_CHK_SET_ERR(rval, "Can't get node coords");

  std::vector<double>  inico(3*xzverts.size());
  rval = mb.get_coords(xzverts, &inico[0]); MB_CHK_ERR(rval);

  // for each layer, compute an angle
  int nv = (int)xzverts.size();
  for (int i=0; i<layers; i++)
  {
    double arad=angleRadInc*(i+1); // for last layer (128?), it will be angle * 180/M_PI
    double deltaz = deltazInc* (i+1); // for last layer it will be total pitch 
    for (int j=0; j<nv; j++)
    {
      double x=inico[3*j];

      double z=inico[3*j+2]; // basically original y is 0!!!
      arrays[0][i*nv+j] = x*cos(arad);
      arrays[1][i*nv+j] = x*sin(arad); // rotate counter-clockwise here, so y becomes positive here for positive x (angle is now positive)
      arrays[2][i*nv+j] = z+deltaz;
    }
  }

  int num_quads = layers * (int)xzedges.size(); // will generate layers of quads for each edge in xz plane

  EntityHandle starte; // Connectivity
  EntityHandle* conn;
  //int num_v_per_elem = 4;

  rval = iface->get_element_connect(num_quads, 4, MBQUAD, 0, starte, conn);MB_CHK_SET_ERR(rval, "Can't get element connectivity for new quads");

  int ixh=0; // index in conn array
  for (Range::iterator eit=xzedges.begin(); eit!=xzedges.end(); eit++)
  {
    EntityHandle edge=*eit;
    const EntityHandle * conn2 =NULL;
    int nno;
    rval = mb.get_connectivity(edge, conn2, nno); MB_CHK_ERR(rval);
    if (nno!=2 )continue;// should error out actually 

    int indices[2];
    for (int k=0; k<2; k++)
      indices[k] = xzverts.index(conn2[k]);

    // first layer quad will use connectivity vertices conn2; next will use the indices
    // 
    conn[ixh++]= conn2[0];
    conn[ixh++]= conn2[1];
    conn[ixh++]= startv  + indices[1];
    conn[ixh++]= startv  + indices[0];
  }
  
  for (int i=2; i<=layers; i++)
  {
    std::cout <<" Layer : " << i << endl;
    for (Range::iterator eit=xzedges.begin(); eit!=xzedges.end(); eit++)
    {
      EntityHandle edge=*eit;
      const EntityHandle * conn2 =NULL;
      int nno;
      rval = mb.get_connectivity(edge, conn2, nno); MB_CHK_ERR(rval);
      if (nno!=2 )continue;// should error out actually 

      int indices[2];
      for (int k=0; k<2; k++)
        indices[k] = xzverts.index(conn2[k]);
     
      conn[ixh++]= startv + indices[0]+(i-2)*nv;
      conn[ixh++]= startv + indices[1]+(i-2)*nv;
      conn[ixh++]= startv + indices[1]+(i-1)*nv;
      conn[ixh++]= startv + indices[0]+(i-1)*nv;
    }
  }

  Tag dir_tag;
  rval = mb.tag_get_handle(MATERIAL_SET_TAG_NAME, dir_tag); MB_CHK_SET_ERR(rval, "Can't get dirichlet tag");

  Tag ptag;
  int dum_id = -1;
  rval = mb.tag_get_handle("PARALLEL_PARTITION", 1, MB_TYPE_INTEGER,
                             ptag, MB_TAG_CREAT | MB_TAG_SPARSE, &dum_id);MB_CHK_SET_ERR(rval, "Can't get parallel partition tag");

  Range sets;
  EntityHandle set1, setSupp1, set2, setSupp2, set3;
  rval = mb.create_meshset(MESHSET_SET, set1);MB_CHK_ERR(rval);
  rval = mb.add_entities(set1, orgQuads);
  sets.insert(set1);
  Range newQuadsSupport1(starte, starte+(int)xzedges.size()*3-1);

  rval = mb.create_meshset(MESHSET_SET, setSupp1);MB_CHK_ERR(rval);
  rval = mb.add_entities(setSupp1, newQuadsSupport1);
  sets.insert(setSupp1);

  int valDirTag = 10;
  rval = mb.tag_set_data(dir_tag, &set1, 1, &valDirTag);MB_CHK_ERR(rval);
  valDirTag = 20;
  rval = mb.tag_set_data(dir_tag, &setSupp1, 1, &valDirTag);MB_CHK_ERR(rval);
  dum_id = 0;
  rval = mb.tag_set_data(ptag, &set1, 1, &dum_id);MB_CHK_ERR(rval);
  dum_id = 1;
  rval = mb.tag_set_data(ptag, &setSupp1, 1, &dum_id);MB_CHK_ERR(rval);

  Range newQuads2(starte+(int)xzedges.size()*3, starte+(int)xzedges.size()*(3+128)-1);
  rval = mb.create_meshset(MESHSET_SET, set2);MB_CHK_ERR(rval);
  rval = mb.add_entities(set2, newQuads2);
  sets.insert(set2);
  valDirTag = 30;
  rval = mb.tag_set_data(dir_tag, &set2, 1, &valDirTag);MB_CHK_ERR(rval);
  dum_id = 2;
  rval = mb.tag_set_data(ptag, &set2, 1, &dum_id);MB_CHK_ERR(rval);

  
  Range newQuadsSupport2(starte+(int)xzedges.size()*(3+128), starte+(int)xzedges.size()*(3+128+3)-1);
  rval = mb.create_meshset(MESHSET_SET, setSupp2);MB_CHK_ERR(rval);
  rval = mb.add_entities(setSupp2, newQuadsSupport2);
  sets.insert(setSupp2);
  valDirTag = 40;
  rval = mb.tag_set_data(dir_tag, &setSupp2, 1, &valDirTag);MB_CHK_ERR(rval);
  dum_id = 3;
  rval = mb.tag_set_data(ptag, &setSupp2, 1, &dum_id);MB_CHK_ERR(rval);

  Range newQuads3(starte+(int)xzedges.size()*(3+128+3), starte+(int)xzedges.size()*(3+128+3+20)-1);
  rval = mb.create_meshset(MESHSET_SET, set3);MB_CHK_ERR(rval);
  rval = mb.add_entities(set3, newQuads3);
  sets.insert(set3);
  valDirTag = 50;
  rval = mb.tag_set_data(dir_tag, &set3, 1, &valDirTag);MB_CHK_ERR(rval);
  dum_id = 4;
  rval = mb.tag_set_data(ptag, &set3, 1, &dum_id);MB_CHK_ERR(rval);
  

  Range otheredges;
   
  
  rval = mb.get_adjacencies(otherverts, 1, false, otheredges, Interface::UNION);
   cout <<" number of edges at 45 degree plane: " << otheredges.size() << " psize: " << otheredges.psize() << endl;

  // int layers = 3 + orglayers + 3 + 20; // 3 layers for support 1, 3 for support 2, 20 free 

  num_newNodes = (layers) * (int)xzverts.size();

  // vector<double*> arrays;
  // EntityHandle startv; // a new set of vertices, in the other direction
  rval = iface->get_node_coords(3, num_newNodes, 0, startv, arrays);MB_CHK_SET_ERR(rval, "Can't get node coords");

  //std::vector<double>  inico(3*xzverts.size()); // the same number of vertices ....
  rval = mb.get_coords(otherverts, &inico[0]); MB_CHK_ERR(rval);

  // for each layer, compute an angle
  // int nv = (int)xzverts.size();
  for (int i=0; i<layers; i++)
  {
    double arad=  - angleRadInc*(i+1); // this angle is now negative (clockwise)
    double deltaz = deltazInc*(i+1); // for last layer it will be total pitch 
    for (int j=0; j<nv; j++)
    {
      double x = inico[3*j];
      double y = inico[3*j+1];
      double z = inico[3*j+2]; // basically original y is 0!!!
      arrays[0][i*nv+j] = x*cos(arad) - y*sin(arad);
      arrays[1][i*nv+j] = x*sin(arad) + y*cos(arad); // rotate clockwise, so y becomes positive here for positive x (angle si now negative)
      arrays[2][i*nv+j] = z - deltaz;
    }
  }

  num_quads = layers * (int)otheredges.size(); // will generate layers of quads for each edge in 45 degree plane

  //EntityHandle starte; // Connectivity
  //EntityHandle* conn;
  //int num_v_per_elem = 4;

  rval = iface->get_element_connect(num_quads, 4, MBQUAD, 0, starte, conn);MB_CHK_SET_ERR(rval, "Can't get element connectivity for new quads");

  /*int */ ixh=0; // index in conn array, reset to 0
  for (Range::iterator eit=otheredges.begin(); eit!=otheredges.end(); eit++)
  {
    EntityHandle edge=*eit;
    const EntityHandle * conn2 =NULL;
    int nno;
    rval = mb.get_connectivity(edge, conn2, nno); MB_CHK_ERR(rval);
    if (nno!=2 )continue;// should error out actually 

    int indices[2];
    for (int k=0; k<2; k++)
      indices[k] = otherverts.index(conn2[k]);

    // first layer quad will use connectivity vertices conn2; next will use the indices
    // 
    conn[ixh++]= conn2[0];
    conn[ixh++]= conn2[1];
    conn[ixh++]= startv  + indices[1];
    conn[ixh++]= startv  + indices[0];
  }
  
  for (int i=2; i<=layers; i++)
  {
    std::cout <<" Layer : " << i << endl;
    for (Range::iterator eit=otheredges.begin(); eit!=otheredges.end(); eit++)
    {
      EntityHandle edge=*eit;
      const EntityHandle * conn2 =NULL;
      int nno;
      rval = mb.get_connectivity(edge, conn2, nno); MB_CHK_ERR(rval);
      if (nno!=2 )continue;// should error out actually 

      int indices[2];
      for (int k=0; k<2; k++)
        indices[k] = otherverts.index(conn2[k]);
     
      conn[ixh++]= startv + indices[0]+(i-2)*nv;
      conn[ixh++]= startv + indices[1]+(i-2)*nv;
      conn[ixh++]= startv + indices[1]+(i-1)*nv;
      conn[ixh++]= startv + indices[0]+(i-1)*nv;
    }
  }

  EntityHandle setSupp3, set4, setSupp4, set5;
  
  Range newQuadsSupport3(starte, starte+(int)otheredges.size()*3 -1);
  rval = mb.create_meshset(MESHSET_SET, setSupp3);MB_CHK_ERR(rval);
  rval = mb.add_entities(setSupp3, newQuadsSupport3);
  sets.insert(setSupp3);
  valDirTag = 60;
  rval = mb.tag_set_data(dir_tag, &setSupp3, 1, &valDirTag);MB_CHK_ERR(rval);
  dum_id = 5;
  rval = mb.tag_set_data(ptag, &setSupp3, 1, &dum_id);MB_CHK_ERR(rval);

  Range newQuads4(starte+(int)otheredges.size()*3, starte+(int)otheredges.size()*(3+128)-1);
  rval = mb.create_meshset(MESHSET_SET, set4);MB_CHK_ERR(rval);
  rval = mb.add_entities(set4, newQuads4);
  sets.insert(set4);
  valDirTag = 70;
  rval = mb.tag_set_data(dir_tag, &set4, 1, &valDirTag);MB_CHK_ERR(rval);
  dum_id = 6;
  rval = mb.tag_set_data(ptag, &set4, 1, &dum_id);MB_CHK_ERR(rval);

  
  Range newQuadsSupport4(starte+(int)otheredges.size()*(3+128), starte+(int)otheredges.size()*(3+128+3)-1);
  rval = mb.create_meshset(MESHSET_SET, setSupp4);MB_CHK_ERR(rval);
  rval = mb.add_entities(setSupp4, newQuadsSupport4);
  sets.insert(setSupp4);
  valDirTag = 80;
  rval = mb.tag_set_data(dir_tag, &setSupp4, 1, &valDirTag);MB_CHK_ERR(rval);
  dum_id = 7;
  rval = mb.tag_set_data(ptag, &setSupp4, 1, &dum_id);MB_CHK_ERR(rval);

  Range newQuads5(starte+(int)otheredges.size()*(3+128+3), starte+(int)otheredges.size()*(3+128+3+20)-1);
  rval = mb.create_meshset(MESHSET_SET, set5);MB_CHK_ERR(rval);
  rval = mb.add_entities(set5, newQuads5);
  sets.insert(set5);
  valDirTag = 90;
  rval = mb.tag_set_data(dir_tag, &set5, 1, &valDirTag);MB_CHK_ERR(rval);
  dum_id = 8;
  rval = mb.tag_set_data(ptag, &set5, 1, &dum_id);MB_CHK_ERR(rval);

  rval = mb.write_file("out.h5m", 0, 0, sets);MB_CHK_ERR(rval);
  return 0;
}
