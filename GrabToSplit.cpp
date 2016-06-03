/** @example SetsNTags.cpp
 * Description: Get the sets representing materials and Dirichlet/Neumann boundary conditions and list their contents.\n
 * This example shows how to get entity sets, and tags on those sets.
 *
 * To run: ./SetsNTags [meshfile]\n
 * (default values can run if users don't specify a mesh file)
 */

#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "MBTagConventions.hpp"
#include "math.h"

#include <iostream>

using namespace moab;
using namespace std;

#ifndef MESH_DIR
#define MESH_DIR "."
#endif

string test_file_name("newSetsScale2.exo");

ErrorCode     opposite(Interface * mb, EntityHandle hex, EntityHandle & opphex, 
 Range & levelVerts, Range & oppVerts)
{
  Range nodes;
  ErrorCode rval = mb->get_connectivity(&hex, 1, nodes); MB_CHK_ERR(rval);
  oppVerts = subtract(nodes, levelVerts);
  Range hexas;
  rval = mb->get_adjacencies(oppVerts, 3, false, hexas); // intersect, so only 2, in general
  hexas.erase(hex);
  if (hexas.size() == 1)
  {
    opphex = hexas[0];
    return MB_SUCCESS;
  }
  return MB_FAILURE; // end of the column probably
}
// Tag names for these conventional tags come from MBTagConventions.hpp
const char *tag_nms[] = {MATERIAL_SET_TAG_NAME, DIRICHLET_SET_TAG_NAME, NEUMANN_SET_TAG_NAME};

int main(int argc, char **argv)
{
  // Get the material set tag handle
  Tag mtag;
  ErrorCode rval;
  Range sets;

  // Get MOAB instance
  Interface* mb = new (std::nothrow) Core;
  if (NULL == mb)
    return 1;

  // Need option handling here for input filename
  if (argc > 1) {
    // User has input a mesh file
    test_file_name = argv[1];
  }

  // Load a file
  rval = mb->load_file(test_file_name.c_str());MB_CHK_ERR(rval);

  rval = mb->tag_get_handle(tag_nms[0], 1, MB_TYPE_INTEGER, mtag);MB_CHK_ERR(rval);
  rval = mb->get_entities_by_type_and_tag(0, MBENTITYSET, &mtag, NULL, 1, sets);MB_CHK_ERR(rval);
  vector<int>  vals;
  vals.resize(sets.size());
  rval = mb->tag_get_data(mtag, sets, &vals[0]);
 
  cout<< " vals and type \n" ; 
  for (int i=0; i<(int)vals.size(); i++)
  {
     Range ents;
     rval = mb->get_entities_by_handle(sets[i], ents); MB_CHK_ERR(rval);
     EntityType etype = mb->type_from_handle(ents[0]);
     cout << "set " << mb->id_from_handle(sets[i]) << " " << vals[i] << " etype:" << etype<< "\n" ;
  }
 
  Tag dtag; 
  rval = mb->tag_get_handle(tag_nms[1], 1, MB_TYPE_INTEGER, dtag);MB_CHK_ERR(rval);
  
  Range dsets; 
  rval = mb->get_entities_by_type_and_tag(0, MBENTITYSET, &dtag, NULL, 1, dsets);MB_CHK_ERR(rval);
  Range ents;
  EntityHandle set1=dsets[0] ; //  set  with vertices
  cout << endl;
  rval = mb->get_entities_by_handle(set1, ents);MB_CHK_ERR(rval);
  std::cout<<" vertices in dset: ents.size():" << ents.size() << "\n";
  vector<double> coordsOrg;
  coordsOrg.resize(ents.size()*3); MB_CHK_ERR(rval);
  rval = mb->get_coords(ents, &coordsOrg[0]); MB_CHK_ERR(rval);
  double topZ = coordsOrg[2];
  cout<<" topZ: " << topZ << "\n";
  Range newHexas;
  rval = mb->get_adjacencies(ents, 3, false, newHexas, Interface::UNION); MB_CHK_ERR(rval);
  Range frontVerts;
  rval = mb->get_connectivity(newHexas, frontVerts); MB_CHK_ERR(rval);
  vector<double> coords;
  frontVerts = subtract(frontVerts, ents);
  coords.resize(frontVerts.size()*3);
  rval = mb->get_coords(frontVerts, &coords[0]); MB_CHK_ERR(rval);
  Range newFront;
  for (int i=0; i<(int) frontVerts.size(); i++)
  {
    EntityHandle newVert= frontVerts[i];
    if ( fabs(coords[3*i+2]-topZ) <=0.00000001)
    {
       newFront.insert(newVert);
    }
  }
  cout << " new front: " << newFront.size() << "\n";
  Range topVerts= unite(ents, newFront);
  Range splitTopVerts;
  Range rightHexas;
  Range leftHexas;
  Range topHexas = newHexas; // accumulate new hexas
  for (int j=0; j<6*8; j++) // add 
  {
    newHexas.clear();
    rval = mb->get_adjacencies(newFront, 3, false, newHexas, Interface::UNION); MB_CHK_ERR(rval);
    newHexas = subtract(newHexas, topHexas);
    topHexas.merge(newHexas);
    frontVerts.clear();
    rval = mb->get_connectivity(newHexas, frontVerts); MB_CHK_ERR(rval);
    frontVerts = subtract(frontVerts, topVerts);
    coords.resize(frontVerts.size()*3);
    rval = mb->get_coords(frontVerts, &coords[0]); MB_CHK_ERR(rval);
    newFront.clear();
    for (int i=0; i<(int) frontVerts.size(); i++)
    {
      EntityHandle newVert= frontVerts[i];
      if ( fabs(coords[3*i+2]-topZ) <=0.00000001)
      {
         newFront.insert(newVert);
      }
    }
    topVerts.merge( newFront);
    if (j%8==6)
    {
      splitTopVerts.merge( newFront);
      rightHexas.merge( newHexas);
    }
    if(j%8==7)
      leftHexas.merge( newHexas);
    
  }
  EntityHandle newSet ;
  int valNew = 160;
  rval = mb->create_meshset(MESHSET_SET, newSet);  MB_CHK_ERR(rval);
  rval = mb-> tag_set_data(dtag, &newSet, 1 , &valNew);MB_CHK_ERR(rval);
  rval = mb->add_entities(newSet, splitTopVerts);

  Range splitHexas = unite(leftHexas, rightHexas);
  
  // now we have splitTopVerts, we need to advance towards bottom, form splitVerts; leave 
  // on left and right the hexas that will be separated
  // find first all verts connected to each top vertex on the same vertical
  Range splitVerts= splitTopVerts;
  coords.resize(splitTopVerts.size()*3);
  rval = mb->get_coords(splitTopVerts, &coords[0]); MB_CHK_ERR(rval);
  for (int j=0; j<(int)splitTopVerts.size(); j++)
  {
    // find all verts on the same vertical, advancing front
    double x = coords[3*j], y = coords[3*j+1];
    Range connHex;
    EntityHandle currentVertex= splitTopVerts[j];
    rval = mb->get_adjacencies(&currentVertex, 1, 3, false, connHex, Interface::UNION);MB_CHK_ERR(rval);
    Range newPotentialVerts;
    rval = mb->get_connectivity(connHex, newPotentialVerts);
    Range alreadySeen= intersect(newPotentialVerts, topVerts);
    newPotentialVerts = subtract(newPotentialVerts, alreadySeen);
    while (1)
    {
      vector<double> xyz; 
      xyz.resize(3*newPotentialVerts.size());
      rval = mb->get_coords(newPotentialVerts, &xyz[0]); MB_CHK_ERR(rval);
      int found = 0;
      for (int k=0; k<(int) newPotentialVerts.size(); k++)
      {
        if( (fabs(xyz[3*k] - x) < 0.000001) && 
            (fabs(xyz[3*k+1] - y) < 0.000001) )
        {
          found = 1;
          currentVertex=newPotentialVerts[k];
          break; // stop at first
        }
      } 

      if(!found)
        break;
      else
      {
        // continue advance
        splitVerts.insert(currentVertex);
        connHex.clear();
        rval = mb->get_adjacencies(&currentVertex, 1, 3, false, connHex);MB_CHK_ERR(rval);
        newPotentialVerts.clear();
        rval = mb->get_connectivity(connHex, newPotentialVerts);
        alreadySeen.insert(currentVertex);
        newPotentialVerts = subtract(newPotentialVerts, alreadySeen);
        splitHexas.merge(connHex);
      }
    }
    cout << " split vertex j: " << j << " x, y " << x << " " << y << " splitVerts.size() " << splitVerts.size() << " hexas: " << splitHexas.size() << "\n";
  }
 
  Range visitedHexas;
  Range queueRight=rightHexas;
  while( !queueRight.empty() )
  {
     EntityHandle currHex = queueRight.pop_back();
     visitedHexas.insert(currHex);
     Range verts;
     rval = mb->get_connectivity(&currHex, 1,  verts); MB_CHK_ERR(rval);
     Range connHex;
     rval = mb->get_adjacencies(verts, 3, false, connHex, Interface::UNION);
     connHex = subtract(connHex, visitedHexas);
     // visitedHexas.merge(connHex);
     // add to the queue only hexas that have 2 common verts with splitVerts
     for (int k=0; k<(int)connHex.size(); k++)
     {
       EntityHandle potRightHex = connHex[k];
       Range commonVerts;
       rval = mb->get_connectivity(&potRightHex, 1,  commonVerts); MB_CHK_ERR(rval);
       commonVerts = intersect (commonVerts, verts);
       if (commonVerts.size()!=4)
         continue;
       commonVerts = intersect(commonVerts, splitVerts);
       if (commonVerts.size() == 2)
       {
         queueRight.insert(potRightHex);
         rightHexas.insert(potRightHex);
       }
     }
  }
 
  valNew = 170;
  EntityHandle rightVset;
  rval = mb->create_meshset(MESHSET_SET, rightVset);  MB_CHK_ERR(rval);
  rval = mb-> tag_set_data(dtag, &rightVset, 1 , &valNew);MB_CHK_ERR(rval);
  Range rVerts;
  rval = mb->get_connectivity(rightHexas, rVerts); MB_CHK_ERR(rval);
  rval = mb->add_entities(rightVset, rVerts); MB_CHK_ERR(rval);

  leftHexas = subtract(splitHexas, rightHexas);
  cout << " leftHexas: " << leftHexas.size() << " rightHexas: " << rightHexas.size() << "\n";
#if 0
  rVerts.clear();
  EntityHandle leftVset;
  rval = mb->create_meshset(MESHSET_SET, leftVset);  MB_CHK_ERR(rval);
  valNew = 180;
  rval = mb-> tag_set_data(dtag, &leftVset, 1 , &valNew);MB_CHK_ERR(rval);
  rval = mb->get_connectivity(leftHexas, rVerts); MB_CHK_ERR(rval);
  rval = mb->add_entities(leftVset, rVerts); MB_CHK_ERR(rval);
  EntityHandle Vset;
  rval = mb->create_meshset(MESHSET_SET, Vset);  MB_CHK_ERR(rval);
  valNew = 190;
  rval = mb-> tag_set_data(dtag, &Vset, 1 , &valNew);MB_CHK_ERR(rval);
  rVerts.clear();
  rval = mb->get_connectivity(splitHexas, rVerts); MB_CHK_ERR(rval);
  rval = mb->add_entities(Vset, rVerts); MB_CHK_ERR(rval);
#endif
  // so now we have leftHexas, rightHexas, that need to be split at splitVerts
  // first duplicate the splitVerts; 
  // then all leftHexas will have new connectivity for splitNodes, efectively decoupling the stuff
  coords.resize(splitVerts.size()*3);
  rval = mb->get_coords(splitVerts, &coords[0]); MB_CHK_ERR(rval);
  Range newVerts;
  int numNewVerts = (int) splitVerts.size();
  rval = mb->create_vertices(&coords[0], numNewVerts, newVerts); MB_CHK_ERR(rval);
  for (Range::iterator eit=rightHexas.begin(); eit!=rightHexas.end(); eit++)
  {
    EntityHandle hex=*eit;
    vector<EntityHandle> vertices;
    rval = mb->get_connectivity(&hex, 1, vertices); MB_CHK_ERR(rval);
    for (int k=0; k<(int)vertices.size(); k++)
    {
      EntityHandle vertex=vertices[k];
      int index=splitVerts.index(vertex);
      if (index!=-1)
      {
        // replace a vertex:
        vertices[k] = newVerts[index];
      }
    }
    rval = mb->set_connectivity(hex, &vertices[0], 8); MB_CHK_ERR(rval);
  }
  
  
#if 0
  Range splitVerts= splitTopVerts;
  Range currSplit = splitTopVerts;
  Range levelVerts= topVerts;
  Range levelHexas = unite(rightHexas, leftHexas);
  Range levelRightHexas = rightHexas;
  Range levelLeftHexas = leftHexas;
  cout << " rightHexas and left Hexas: " << rightHexas.size() << " " << leftHexas.size() << "\n";
  
  for (int i=0; i< 300; i++) // probably only 2*30 + 5*40 ~= 60?
  {
    Range  rightVerts, leftVerts;
    Range rightNextHexas, leftNextHexas;
    int flag = 1; // keep going
    for (int k = 0; k< (int) levelRightHexas.size(); k++)
    {
      EntityHandle hex = levelRightHexas[k];
      EntityHandle opphex;
      Range oppVerts;
      rval = opposite(mb, hex, opphex, levelVerts, oppVerts); 
      if (MB_FAILURE == rval)
      {
        flag = 0;
        break; // we ar eprobably done 
      }
      rightVerts.merge(oppVerts);
      rightNextHexas.insert(opphex);
    } 
    for (int k = 0; k< (int) levelLeftHexas.size(); k++)
    {
      EntityHandle hex = levelLeftHexas[k];
      EntityHandle opphex;
      Range oppVerts;
      rval = opposite(mb, hex, opphex, levelVerts, oppVerts); 
      if (MB_FAILURE == rval)
      {
        flag = 0;
        break; // we ar eprobably done 
      }
      leftVerts.merge(oppVerts);
      leftNextHexas.insert(opphex);
    } 

    if (0==flag)
       break; // so we are done advancing along seams
    rightHexas.merge(rightNextHexas);
    leftHexas.merge(leftNextHexas);
    levelRightHexas = rightNextHexas;
    levelLeftHexas = leftNextHexas;
    
    levelVerts = unite(leftVerts, rightVerts); 
    currSplit = intersect(leftVerts, rightVerts);
    splitVerts.merge(currSplit);
    cout << " level i: " << i << " currSplit.size() = " <<  currSplit.size()  << "\n";
  }
 #endif 
  rval = mb->add_entities(newSet, splitVerts);

  mb->write_file("splits.exo");
  return 0;
}
