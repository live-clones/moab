/*
 * create a helix mesh;
 *  input file : surface side mesh 
 * 
 *  HexaSweepHelix -i start.exo 
 */
#include "moab/Core.hpp"
#include "moab/ProgOptions.hpp"

#include "moab/ReadUtilIface.hpp"
#include "moab/CartVect.hpp"
#include "MBTagConventions.hpp"

#include <time.h>
#include <iostream>
#include <vector>

using namespace moab;
using namespace std;

ErrorCode extrude_quads(Core & mb, ReadUtilIface * iface, Range & quads, int & layers,
  double angleIncSupport, double deltazIncSupport, Tag mtag, int valTag, EntityHandle & startv, 
     EntityHandle & starte, EntityHandle & outset)
{
  Range verts;
  ErrorCode rval = mb.get_connectivity(quads, verts); MB_CHK_SET_ERR(rval, "Can't get vertices on quads 50");
  cout <<" number of vertices: " << verts.size() << endl;

  vector<double> coords;
  coords.resize(verts.size()*3);
  rval = mb.get_coords(verts, &coords[0]); MB_CHK_SET_ERR(rval, "Can't get coords of vertices");

  int num_newNodes = (layers) * (int)verts.size();

  vector<double*> arrays;

  rval = iface  ->get_node_coords(3, num_newNodes, 0, startv, arrays);MB_CHK_SET_ERR(rval, "Can't get node coords");

  std::vector<double>  inico(3*verts.size());
  rval = mb.get_coords(verts, &inico[0]); MB_CHK_ERR(rval);

  // for each layer, compute an angle and a delta z
  int nv = (int)verts.size();
  for (int i=0; i<layers; i++)
  {
    double arad=angleIncSupport*(i+1); // for last layer (8?), it will be angle * 180/M_PI
    double deltaz = deltazIncSupport* (i+1); // for last layer it will be total pitch 
    
    for (int j=0; j<nv; j++)
    {
      double x=inico[3*j];
      double y=inico[3*j+1];
      double z=inico[3*j+2]; // basically 
      arrays[0][i*nv+j] = x*cos(arad) - y * sin(arad);
      arrays[1][i*nv+j] = x*sin(arad) + y*cos(arad); // rotate clockwise here, 
      arrays[2][i*nv+j] = z+deltaz;
    }
  }

  int num_hexas = (layers) * (int)quads.size(); // will generate layers of hexa for each quad in xz plane

  //EntityHandle starte; // Connectivity
  EntityHandle* conn;
  //int num_v_per_elem = 8;

  rval = iface->get_element_connect(num_hexas, 8, MBHEX, 0, starte, conn);MB_CHK_SET_ERR(rval, "Can't get element connectivity for new hexas");

  int ixh=0; // index in conn array
  for (Range::iterator eit=quads.begin(); eit!=quads.end(); eit++)
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
  
  for (int i=2; i<=layers; i++) // 
  {
    std::cout <<" Layer : " << i << endl;
    for (Range::iterator eit=quads.begin(); eit!=quads.end(); eit++)
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


  Range newHex(starte, starte+(int)quads.size()*layers-1);
  //EntityHandle setSupp1;
  rval = mb.create_meshset(MESHSET_SET, outset);MB_CHK_ERR(rval);
  rval = mb.add_entities(outset, newHex);MB_CHK_ERR(rval);
  rval = mb.tag_set_data(mtag, &outset, 1, &valTag); MB_CHK_ERR(rval);
  return MB_SUCCESS;
}

ErrorCode extrude_edges(Core & mb, ReadUtilIface * iface, EntityHandle eset, int & layers,
  double angleIncSupport, double deltazIncSupport, Tag mtag, int valTag, EntityHandle & outset)
{
  EntityHandle startv , starte;
  // now add quads 
  Range edges;
  ErrorCode rval = mb.get_entities_by_dimension(eset, 1, edges);MB_CHK_SET_ERR(rval, "Can't get edges from set");

  Range verts;
  rval = mb.get_connectivity(edges, verts); MB_CHK_SET_ERR(rval, "Can't get vertices on verts");
  cout <<"Extrude edges:  number of vertices: " << verts.size() << endl;

  vector<double> coords;
  coords.resize(verts.size()*3);
  rval = mb.get_coords(verts, &coords[0]); MB_CHK_SET_ERR(rval, "Can't get coords of vertices");

  int num_newNodes = (layers) * (int)verts.size();

  vector<double*> arrays;

  rval = iface  ->get_node_coords(3, num_newNodes, 0, startv, arrays);MB_CHK_SET_ERR(rval, "Can't get node coords");

  std::vector<double>  inico(3*verts.size());
  rval = mb.get_coords(verts, &inico[0]); MB_CHK_ERR(rval);

  // for each layer, compute an angle and a delta z
  int nv = (int)verts.size();
  for (int i=0; i<layers; i++)
  {
    double arad=angleIncSupport*(i+1); // for last layer (8?), it will be angle * 180/M_PI
    double deltaz = deltazIncSupport* (i+1); // for last layer it will be total pitch 
    
    for (int j=0; j<nv; j++)
    {
      double x=inico[3*j];
      double y=inico[3*j+1];
      double z=inico[3*j+2]; // basically 
      arrays[0][i*nv+j] = x*cos(arad) - y * sin(arad);
      arrays[1][i*nv+j] = x*sin(arad) + y*cos(arad); // rotate clockwise here, 
      arrays[2][i*nv+j] = z+deltaz;
    }
  }
  // we have the edges, generate quads
  int num_quads = (layers) * (int)edges.size(); // will generate layers of quads for each edge in set

  //EntityHandle starte; // Connectivity
  EntityHandle* conn;
  //int num_v_per_elem = 8;

  rval = iface->get_element_connect(num_quads, 4, MBQUAD, 0, starte, conn);MB_CHK_SET_ERR(rval, "Can't get element connectivity for new quads");

  int ixh=0; // index in conn array
  for (Range::iterator eit=edges.begin(); eit!=edges.end(); eit++)
  {
    EntityHandle edge=*eit;
    const EntityHandle * conn2 =NULL;
    int nno;
    rval = mb.get_connectivity(edge, conn2, nno); MB_CHK_ERR(rval);
    if (nno!=2 )continue;// should error out actually 

    int indices[2];
    for (int k=0; k<2; k++)
      indices[k] = verts.index(conn2[k]);

    // first layer hexa will use connectivity vertices conn4; next will use the indices
    // 
    conn[ixh++]= conn2[0];
    conn[ixh++]= conn2[1];
    conn[ixh++]= startv  + indices[1];
    conn[ixh++]= startv  + indices[0];
  }
  for (int i=2; i<=layers; i++) // 
  {
    std::cout <<" Layer quad: " << i << endl;
    for (Range::iterator eit=edges.begin(); eit!=edges.end(); eit++)
    {
      EntityHandle edge=*eit;
      const EntityHandle * conn2 =NULL;
      int nno;
      rval = mb.get_connectivity(edge, conn2, nno); MB_CHK_ERR(rval);
      if (nno!=2 )continue;// should error out actually 

      int indices[2];
      for (int k=0; k<2; k++)
        indices[k] = verts.index(conn2[k]);

	    // first layer hexa will use connectivity vertices conn4; next will use the indices
	    // 
      conn[ixh++]= startv + indices[0]+(i-2)*nv;
      conn[ixh++]= startv + indices[1]+(i-2)*nv;
      conn[ixh++]=  startv + indices[1]+(i-1)*nv;
      conn[ixh++]=  startv + indices[0]+(i-1)*nv;
    }

  }


  Range newQuads(starte, starte+(int)edges.size()*layers-1);
  //EntityHandle setSupp1;
  rval = mb.create_meshset(MESHSET_SET, outset);MB_CHK_ERR(rval);
  rval = mb.add_entities(outset, newQuads);MB_CHK_ERR(rval);
  int val = valTag;
  rval = mb.tag_set_data(mtag, &outset, 1, &val); MB_CHK_ERR(rval);
// end copy
  
  return MB_SUCCESS;
}


ErrorCode extrude_quads_flat(Core & mb, ReadUtilIface * iface, EntityHandle qset, 
  EntityHandle eset, int & layers,
  double distFlat, double tpitch, Tag mtag, int valTag,  EntityHandle & outset, EntityHandle & outquadset,
     EntityHandle & outEset)
{
  EntityHandle  startv, starte; // local here
  Range quads;
 
  ErrorCode rval = mb.get_entities_by_dimension(qset, 2, quads); MB_CHK_SET_ERR(rval, "Can't get quads");
  cout <<" number of quads in set  " << quads.size() << endl;
  Range verts;
  rval = mb.get_connectivity(quads, verts); MB_CHK_SET_ERR(rval, "Can't get vertices on quads ");
  cout <<" number of vertices: " << verts.size() << endl;

  vector<double> coords;
  coords.resize(verts.size()*3);
  rval = mb.get_coords(verts, &coords[0]); MB_CHK_SET_ERR(rval, "Can't get coords of vertices");

  int num_newNodes = (layers) * (int)verts.size();

  vector<double*> arrays;

  rval = iface  ->get_node_coords(3, num_newNodes, 0, startv, arrays);MB_CHK_SET_ERR(rval, "Can't get node coords");

  std::vector<double>  inico(3*verts.size());
  rval = mb.get_coords(verts, &inico[0]); MB_CHK_ERR(rval);

  // for each layer, compute an angle and a delta z
  int nv = (int)verts.size();
  for (int i=0; i<layers; i++)
  {
    double heightLayer=distFlat*(i+1)/layers;
    //double arad=angleIncSupport*(i+1); // for last layer (8?), it will be angle * 180/M_PI
    //double deltaz = deltazIncSupport* (i+1); // for last layer it will be total pitch 
    
    for (int j=0; j<nv; j++)
    {
      double x=inico[3*j];
      double y=inico[3*j+1];
      double z=inico[3*j+2]; 
      // first compute angle from x, y, along circumference 
      double arad = heightLayer/sqrt( x*x+y*y); // these are radians; for small angles, these are fine
      // now what is z? 
      double deltaz = tpitch/(2*M_PI)*arad; // total pitch is for 360 degrees = 2*M_PI
      arrays[0][i*nv+j] = x*cos(arad) - y * sin(arad);
      arrays[1][i*nv+j] = x*sin(arad) + y*cos(arad); // rotate clockwise here, 
      arrays[2][i*nv+j] = z+deltaz;
    }
  }

  int num_hexas = (layers) * (int)quads.size(); // will generate layers of hexa for each quad in xz plane

  //EntityHandle starte; // Connectivity
  EntityHandle* conn;
  //int num_v_per_elem = 8;

  rval = iface->get_element_connect(num_hexas, 8, MBHEX, 0, starte, conn);MB_CHK_SET_ERR(rval, "Can't get element connectivity for new hexas");

  int ixh=0; // index in conn array
  for (Range::iterator eit=quads.begin(); eit!=quads.end(); eit++)
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
  
  for (int i=2; i<=layers; i++) // 
  {
    std::cout <<" Layer : " << i << endl;
    for (Range::iterator eit=quads.begin(); eit!=quads.end(); eit++)
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


  Range newHex(starte, starte+(int)quads.size()*layers-1);
  //EntityHandle setSupp1;
  rval = mb.create_meshset(MESHSET_SET, outset);MB_CHK_ERR(rval);
  rval = mb.add_entities(outset, newHex);MB_CHK_ERR(rval);
  rval = mb.tag_set_data(mtag, &outset, 1, &valTag); MB_CHK_ERR(rval);

  // now add quads from edge set
  Range edges;
  rval = mb.get_entities_by_dimension(eset, 1, edges);MB_CHK_SET_ERR(rval, "Can't get edges from set");

  
  int num_quads = (layers) * (int)edges.size(); // will generate layers of quads for each edge in set

  //EntityHandle starte; // Connectivity
  // EntityHandle* conn;
  //int num_v_per_elem = 8;

  rval = iface->get_element_connect(num_quads, 4, MBQUAD, 0, starte, conn);MB_CHK_SET_ERR(rval, "Can't get element connectivity for new quads");

  ixh=0; // index in conn array
  for (Range::iterator eit=edges.begin(); eit!=edges.end(); eit++)
  {
    EntityHandle edge=*eit;
    const EntityHandle * conn2 =NULL;
    int nno;
    rval = mb.get_connectivity(edge, conn2, nno); MB_CHK_ERR(rval);
    if (nno!=2 )continue;// should error out actually 

    int indices[2];
    for (int k=0; k<2; k++)
      indices[k] = verts.index(conn2[k]);

    // first layer hexa will use connectivity vertices conn4; next will use the indices
    // 
    conn[ixh++]= conn2[0];
    conn[ixh++]= conn2[1];
    conn[ixh++]= startv  + indices[1];
    conn[ixh++]= startv  + indices[0];
  }
  for (int i=2; i<=layers; i++) // 
  {
    std::cout <<" Layer quad: " << i << endl;
    for (Range::iterator eit=edges.begin(); eit!=edges.end(); eit++)
    {
      EntityHandle edge=*eit;
      const EntityHandle * conn2 =NULL;
      int nno;
      rval = mb.get_connectivity(edge, conn2, nno); MB_CHK_ERR(rval);
      if (nno!=2 )continue;// should error out actually 

      int indices[2];
      for (int k=0; k<2; k++)
        indices[k] = verts.index(conn2[k]);

	    // first layer hexa will use connectivity vertices conn4; next will use the indices
	    // 
      conn[ixh++]= startv + indices[0]+(i-2)*nv;
      conn[ixh++]= startv + indices[1]+(i-2)*nv;
      conn[ixh++]=  startv + indices[1]+(i-1)*nv;
      conn[ixh++]=  startv + indices[0]+(i-1)*nv;
    }

  }


  Range newQuads(starte, starte+(int)edges.size()*layers-1);
  //EntityHandle setSupp1;
  rval = mb.create_meshset(MESHSET_SET, outquadset);MB_CHK_ERR(rval);
  rval = mb.add_entities(outquadset, newQuads);MB_CHK_ERR(rval);
  int val = valTag+10;
  rval = mb.tag_set_data(mtag, &outquadset, 1, &val); MB_CHK_ERR(rval);

  // add new edges, in a set  
  int num_edges = (int)edges.size(); // will generate 1 set of edges at end

  //EntityHandle starte; // Connectivity
  // EntityHandle* conn;
  //int num_v_per_elem = 8;

  rval = iface->get_element_connect(num_edges, 2, MBEDGE, 0, starte, conn);MB_CHK_SET_ERR(rval, "Can't get element connectivity for new edges");
  ixh=0; // index in conn array
  for (Range::iterator eit=edges.begin(); eit!=edges.end(); eit++)
  {
    EntityHandle edge=*eit;
    const EntityHandle * conn2 =NULL;
    int nno;
    rval = mb.get_connectivity(edge, conn2, nno); MB_CHK_ERR(rval);
    if (nno!=2 )continue;// should error out actually 

    int indices[2];
    for (int k=0; k<2; k++)
      indices[k] = verts.index(conn2[k]);

    // first layer hexa will use connectivity vertices conn4; next will use the indices
    // 
    conn[ixh++]= startv  + (layers-1)*nv + indices[1]; // reverse the order?
    conn[ixh++]= startv  + (layers-1)*nv + indices[0]; // reverse the order?
  }

  Range newEdges(starte, starte+(int)edges.size() - 1);

  rval = mb.create_meshset(MESHSET_SET, outEset);MB_CHK_ERR(rval);
  rval = mb.add_entities(outEset, newEdges);MB_CHK_ERR(rval);
  val = valTag+20;
  rval = mb.tag_set_data(mtag, &outEset, 1, &val); MB_CHK_ERR(rval);

  return MB_SUCCESS;
}


int main(int argc, char **argv)
{
  double pitch =28.9284 *0.0254 ; // full revolution, in meters
  double delta1 = 0.0367474744147 ; // angle for support between sections 2 and 3
  double delta2 = 0.0367474742259 ; // angle for support between section 3 and 4
  int orglayers = 100; // 
  int numLayersSupport = 8; // 
  ProgOptions opts;

  string inputf =  "start.exo";
  opts.addOpt<string>("input,i", "Specify the input file name string (default start.exo)", &inputf);
  opts.addOpt<double>(string("pitch,p"), string("Total pitch (default=0.73478136)"), &pitch);
  //opts.addOpt<double>(string("angle,a"),  string("Total angle (default=-45.)"), &angle);
  opts.addOpt<int>(string("supplayers,l"),  string("num layers (default=8)"), &numLayersSupport);

  opts.parseCommandLine(argc, argv);

  //double angleRadInc = angle/180.*M_PI/orglayers; // 
  //double deltazInc =  pitch/orglayers* (angle/360.) ;

  double supportAngle = delta1*180/M_PI; // in degrees
  double angleIncSupport = supportAngle * M_PI/180./numLayersSupport; // so total support angle will be about delta1
  double deltazIncSupport = pitch/2/M_PI*angleIncSupport; // 
  
  cout << " support angle: " << supportAngle << " degrees\n"; 
  cout << "angleIncSupport: " << angleIncSupport << " rad\n";
  cout << " deltazIncSupport: " << deltazIncSupport << " m\n";
  Core mb;

  ErrorCode rval = mb.load_file(inputf.c_str()); MB_CHK_ERR(rval);

  ReadUtilIface* iface;
  rval = mb.query_interface(iface);MB_CHK_SET_ERR(rval, "Can't get reader interface");

  Tag mtag; // 
  rval = mb.tag_get_handle(MATERIAL_SET_TAG_NAME, 1, MB_TYPE_INTEGER, mtag);MB_CHK_ERR(rval);
  Range sets;
  rval = mb.get_entities_by_type_and_tag(0, MBENTITYSET, &mtag, NULL, 1, sets);MB_CHK_ERR(rval);

  vector<int>  vals;
  vals.resize(sets.size());
  rval = mb.tag_get_data(mtag, sets, &vals[0]);
 
  cout<< " vals: \n" ; 
  for (int i=0; i<(int)vals.size(); i++)
     cout << "set " << mb.id_from_handle(sets[i]) << " " << vals[i] << "\n" ;
  
  EntityHandle set50 = sets[4]; // has material tag 50
  EntityHandle set40 = sets[3]; // edge set, around tubes (left side of set 10)
  
  

  EntityHandle setSupp1, setQuad1, setEdge1;
  double flatDist = 0.04463; // 1.75825*0.0254? first was with 0.04465
/*
ErrorCode extrude_quads_flat(Core & mb, ReadUtilIface * iface, EntityHandle qset, 
  EntityHandle eset, int & layers,
  double distFlat, double tpitch, Tag mtag, int valTag, int vedTag, EntityHandle & startv, 
     EntityHandle & starte, EntityHandle & outset, EntityHandle & outquadset,
     EntityHandle & outEset)
*/
  rval = extrude_quads_flat( mb, iface, set50, set40, numLayersSupport,
    flatDist, /*total pitch*/ pitch, mtag, /*val tag*/ 140, 
    setSupp1, setQuad1, setEdge1); MB_CHK_SET_ERR(rval, "Can't extrude quads");
  sets.insert(setSupp1);
  sets.insert(setQuad1);
  sets.insert(setEdge1);

  EntityHandle setE20 = sets[1];
  EntityHandle setQ30 = sets[2];
  EntityHandle setSupp2, setQ2, setE2;
  rval = extrude_quads_flat( mb, iface, setQ30, setE20, numLayersSupport,
    -flatDist, /*total pitch*/ pitch, mtag, /*val tag*/ 170, 
    setSupp2, setQ2, setE2); MB_CHK_SET_ERR(rval, "Can't extrude quads");
  sets.insert(setSupp2);
  sets.insert(setQ2);
  sets.insert(setE2);

  EntityHandle setE80 = sets[7];
  EntityHandle setQ90 = sets[8];
  EntityHandle setSupp3, setQ3, setE3;
  rval = extrude_quads_flat( mb, iface, setQ90, setE80, numLayersSupport,
    flatDist, /*total pitch*/ pitch, mtag, /*val tag*/ 200, 
    setSupp3, setQ3, setE3); MB_CHK_SET_ERR(rval, "Can't extrude quads");
  sets.insert(setSupp3);
  sets.insert(setQ3);
  sets.insert(setE3);

  EntityHandle setE120 = sets[11];
  EntityHandle setQ130 = sets[12];
  EntityHandle setSupp4, setQ4, setE4;
  rval = extrude_quads_flat( mb, iface, setQ130, setE120, numLayersSupport,
    -flatDist, /*total pitch*/ pitch, mtag, /*val tag*/ 230, 
    setSupp4, setQ4, setE4); MB_CHK_SET_ERR(rval, "Can't extrude quads");
  sets.insert(setSupp4);
  sets.insert(setQ4);
  sets.insert(setE4);

/*
  ErrorCode extrude_edges(Core & mb, ReadUtilIface * iface, EntityHandle eset, int & layers,
  double angleIncSupport, double deltazIncSupport, Tag mtag, int valTag, EntityHandle & outset)
*/
  EntityHandle setQ5; // from setE3
  int qlayers = 80;
  rval = extrude_edges(mb, iface, setE3, qlayers,
  angleIncSupport, deltazIncSupport, mtag, 260, setQ5);
  sets.insert(setQ5);

  EntityHandle setQ6; // from setE4

  rval = extrude_edges(mb, iface, setE4, qlayers,
  -angleIncSupport, -deltazIncSupport, mtag, 270, setQ6);
  sets.insert(setQ6);

  rval = mb.write_file("new2.exo");MB_CHK_ERR(rval);
 

  return 0;
}
