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
#include "moab/AdaptiveKDTree.hpp"
#include "moab/CartVect.hpp"

#include <time.h>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace moab;
using namespace std;

class Rounds 
{
 public:
  Rounds(Interface * mb, vector<EntityHandle> cir) : mbi(mb) {
   edges=cir;
   // compute center, then order vertices such that they start from top, towards smaller y
   Range verts;
   mbi->get_connectivity(&edges[0], edges.size(), verts); 
   vector<CartVect> coords;
   center = CartVect(0.);
   coords.resize(verts.size());
   mb->get_coords(verts, &(coords[0][0]));
   for (int i=0; i< (int)verts.size(); i++)
     center = center + coords[i];
   center = center/verts.size();
   //cout << center << "\n";
   // find max z vertex; then start from there in order
   double maxz = -2;
   double imaxz=-1;
   for (int j=0; j<(int)verts.size(); j++)
   {
     if (coords[j][2] > maxz)
     {
       maxz = coords[j][2] ;
       imaxz = j;
     }
   }
   EntityHandle startv = verts[imaxz];
   vertices.push_back(startv);
   vector<EntityHandle> adjs;
   mb->get_adjacencies(&startv,1, /* dim*/1, false, adjs);
   Range verts2 = verts;
   verts2.erase(startv);
   const EntityHandle * conn2 = NULL;
   int nnodes;
   mb->get_connectivity(adjs[0], conn2, nnodes); 
   EntityHandle nextv = conn2[0];
   if (startv == nextv) 
     nextv = conn2[1];
   int ix = verts.index(nextv);
   if ( coords[imaxz][1] > coords[ix][1] )
   { 
     // we are fine
     vertices.push_back(nextv); 
     verts2.erase(nextv);
     
   }
   else // try adj1
   {
     mb->get_connectivity(adjs[1], conn2, nnodes); 
     nextv = conn2[0];
     if (startv == nextv) 
       nextv = conn2[1];
     ix = verts.index(nextv);
     if ( coords[imaxz][1] > coords[ix][1] )
     { 
       // fine
       vertices.push_back(nextv); 
       verts2.erase(nextv);
     }
     else 
        cout << " big problem 2 \n"; 
   }
   while(!verts2.empty())
   {
     adjs.clear();
     mb->get_adjacencies(&nextv,1, /* dim*/1, false, adjs);
     Range vv;
     mb->get_connectivity(&adjs[0], adjs.size(), vv);
     for (int k=0; k<(int)vv.size(); k++)
     {
       EntityHandle nnv = vv[k];
       if ( verts2.index(nnv) >=0 )
       {
         nextv = nnv;
         vertices.push_back(nextv);
         verts2.erase(nextv);
       }
     } 
   }
   
  }
  
  Interface * mbi;
  CartVect center;
  vector<EntityHandle> edges;
  vector<EntityHandle> vertices;
  
};

ErrorCode  break_edgeset(Core & mb, EntityHandle eset,  vector< Rounds > & rounds)
{
  vector< vector<EntityHandle> >  circles;
  Range edges;
  ErrorCode rval = mb.get_entities_by_dimension(eset, 1, edges); MB_CHK_ERR(rval);
  Range verts;
  rval = mb.get_connectivity(edges, verts); MB_CHK_ERR(rval);
  cout<< " edges: " << edges.size() << " verts: " << verts.size() << "\n";
  // form loops of edges
  while (!edges.empty() )
  {
    vector<EntityHandle> circle;
    circle.push_back( edges.pop_front() ); // remove first from edges
    
    while(1)
    {
      EntityHandle edge = circle[circle.size()-1]; // current edge
      const EntityHandle * conn2 = NULL;
      int nnodes;
      rval = mb.get_connectivity(edge, conn2, nnodes);  MB_CHK_ERR(rval);
      vector<EntityHandle> adjs;
      rval = mb.get_adjacencies(&conn2[0],1, /* dim*/1, false, adjs);
      if (adjs.size() !=2 ) return MB_FAILURE;
      if ( edges.find(adjs[0]) !=edges.end() )
      {
        edges.erase(adjs[0]);
        circle.push_back(adjs[0]);
      }
      else if (edges.find(adjs[1]) !=edges.end() )
      {
        edges.erase(adjs[1]);
        circle.push_back(adjs[1]);
      }
      else
        break;
    }
    circles.push_back(circle);
  }
  cout<< "circles.size() " << circles.size() << "\n";
  for (int i=0; i< (int)circles.size(); i++)
  {
    rounds.push_back( Rounds(&mb, circles[i]) );
  }
  return MB_SUCCESS;
}

ErrorCode merge_two_vertices(Core & mb, EntityHandle keepv, EntityHandle rmv)
{
  std::vector<EntityHandle> adjs, all;
  ErrorCode rval = mb.get_adjacencies(&rmv, 1, 1, false, adjs); MB_CHK_ERR(rval);
  all = adjs;
  adjs.clear();
  mb.get_adjacencies(&rmv, 1, 2, false, adjs); MB_CHK_ERR(rval);
  all.insert( all.end(), adjs.begin(), adjs.end() );
  adjs.clear();
  mb.get_adjacencies(&rmv, 1, 3, false, adjs); MB_CHK_ERR(rval);
  all.insert( all.end(), adjs.begin(), adjs.end() );
  for (int i=0; i< (int)all.size(); i++)
  {
      // instead of adjs replace in connectivity
    vector<EntityHandle> connect;
   
    rval = mb.get_connectivity(&all[i], 1, connect);  MB_CHK_ERR(rval);
    replace(connect.begin(), connect.end(), rmv, keepv);
    rval = mb.set_connectivity(all[i], &connect[0], connect.size());  MB_CHK_ERR(rval);
   
  }
  mb.delete_entities(&rmv, 1);
  return rval;
}
int main(int argc, char **argv)
{
  ProgOptions opts;

  string inputf =  "new2.exo";
  opts.addOpt<string>("input,i", "Specify the input file name string (default start.exo)", &inputf);
  //opts.addOpt<double>(string("pitch,p"), string("Total pitch (default=0.73478136)"), &pitch);
  //opts.addOpt<double>(string("angle,a"),  string("Total angle (default=-45.)"), &angle);
  //opts.addOpt<int>(string("supplayers,l"),  string("num layers (default=8)"), &numLayersSupport);

  opts.parseCommandLine(argc, argv);
  Core mb;

  ErrorCode rval = mb.load_file(inputf.c_str()); MB_CHK_ERR(rval);

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
  
  EntityHandle set60 = sets[5]; // has material tag 60
  EntityHandle set160 = sets[15]; // edge set, around tubes ; this needs to collapse in set60
 
  vector< Rounds >  rounds;
  rval = break_edgeset(mb, set60, rounds); MB_CHK_ERR(rval);

  vector< Rounds >  rounds1;
  rval = break_edgeset(mb, set160, rounds1); MB_CHK_ERR(rval);

  for (int i=0; i<(int)rounds.size(); i++)
  {
    double minDist=1;
    int idx = -1;
    CartVect center=rounds[i].center;
    for(int j=0; j< (int)rounds1.size(); j++)
    {
      double dist = (center-rounds1[j].center).length();
      if (  dist< minDist )
      {
        minDist = dist;
        idx = j;
      }
    }
    cout << i << " center: " << center << " min dist " << minDist << " at idx: " <<idx<< "  esize: " << 
          rounds[i].edges.size() << " " << rounds1[idx].edges.size() << 
          " vsize: " << rounds[i].vertices.size() << " " << rounds1[idx].vertices.size() <<  "\n";
    if (rounds[i].edges.size() != rounds1[idx].edges.size())
       cout << " big problem\n";
    for ( int j=0; j< (int) rounds[i].vertices.size(); j++)
    {
      rval = merge_two_vertices(mb, rounds[i].vertices[j], rounds1[idx].vertices[j]);  MB_CHK_ERR(rval);
      
    }
  }
  EntityHandle set100 = sets[9]; // has material tag 100
  EntityHandle set190 = sets[18]; // edge set, around tubes ; this needs to collapse in set190
 
  vector< Rounds >  rounds2;
  rval = break_edgeset(mb, set100, rounds2); MB_CHK_ERR(rval);

  vector< Rounds >  rounds3;
  rval = break_edgeset(mb, set190, rounds3); MB_CHK_ERR(rval);

  for (int i=0; i<(int)rounds2.size(); i++)
  {
    double minDist=1;
    int idx = -1;
    CartVect center=rounds2[i].center;
    for(int j=0; j< (int)rounds3.size(); j++)
    {
      double dist = (center-rounds3[j].center).length();
      if ( dist< minDist )
      {
        minDist = dist;
        idx = j;
      }
    }
    cout << i << " center: " << center << " min dist " << minDist << " at idx: " <<idx<< "  isize: " << 
          rounds2[i].edges.size() << " " << rounds3[idx].edges.size() << 
          " vsize: " << rounds2[i].vertices.size() << " " << rounds3[idx].vertices.size() <<  "\n";
    if (rounds2[i].edges.size() != rounds3[idx].edges.size())
      cout << " big problem\n";
    for ( int j=0; j< (int) rounds2[i].vertices.size(); j++)
    {
      rval = merge_two_vertices(mb, rounds2[i].vertices[j], rounds3[idx].vertices[j]);  MB_CHK_ERR(rval);
      
    }
  }

  // delete edge sets
  Range edges;
  rval = mb.get_entities_by_dimension(0, 1, edges);MB_CHK_ERR(rval);
  rval = mb.delete_entities(edges);
  Range extraQuads;
  //EntityHandle qsets[4]={sets[2], sets[4], sets[8], sets[10]};
  rval =  mb.get_entities_by_dimension(sets[2], 2, extraQuads);MB_CHK_ERR(rval);
  rval = mb.delete_entities(extraQuads);
  extraQuads.clear();
  rval =  mb.get_entities_by_dimension(sets[4], 2, extraQuads);MB_CHK_ERR(rval);
  rval = mb.delete_entities(extraQuads);
  extraQuads.clear();
rval =  mb.get_entities_by_dimension(sets[8], 2, extraQuads);MB_CHK_ERR(rval);
  rval = mb.delete_entities(extraQuads);
  extraQuads.clear();
rval =  mb.get_entities_by_dimension(sets[12], 2, extraQuads);MB_CHK_ERR(rval);
  rval = mb.delete_entities(extraQuads);
  extraQuads.clear();
  rval = mb.write_file("new4.exo");MB_CHK_ERR(rval);
 

  return 0;
}
