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

#include <iostream>

using namespace moab;
using namespace std;

#ifndef MESH_DIR
#define MESH_DIR "."
#endif

string test_file_name("extra.exo");

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
  
  Range ents;
  EntityHandle set1=sets[13] ; //  set with 140?
  cout << endl;
  rval = mb->get_entities_by_handle(set1, ents, true);MB_CHK_ERR(rval);
  Range extraEnts = ents;
  std::cout<<" ents.size():" << ents.size() << "\n";
  Range verts;
  rval = mb->get_connectivity(ents, verts); MB_CHK_ERR(rval);
  cout << " verts size: " << verts.size() << endl;
  Range newQuads;
  rval = mb->get_adjacencies(verts, 2, false, newQuads, Interface::UNION);
  Range frontQuads =subtract(newQuads, ents);
  //rval = mb->delete_entities(ents);  MB_CHK_ERR(rval);
  rval = mb->delete_entities(&set1, 1);  MB_CHK_ERR(rval);
  int i = 0;
  while (!frontQuads.empty() )
  {
    i++;
    cout <<"iteration: " << i <<  " frontQuads.size() : " << frontQuads.size() << " verts.size() " << verts.size() << "\n";
    verts.clear();
    rval = mb->get_connectivity(newQuads, verts); MB_CHK_ERR(rval);
    ents= newQuads;
    newQuads.empty();
    rval = mb->get_adjacencies(verts, 2, false, newQuads, Interface::UNION);
    frontQuads =subtract(newQuads, ents); 
  }
  newQuads = subtract(newQuads, extraEnts); // 
  rval = mb->delete_entities(extraEnts);  MB_CHK_ERR(rval);
  cout << " total new quads : " << newQuads.size() << "\n";
  for (i=0; i<9; i++)
  {
    ents.clear();
    rval = mb->get_entities_by_handle(sets[i], ents); MB_CHK_ERR(rval);
    int valnew = vals[i]+2;
    EntityHandle newSet ;
    rval = mb->create_meshset(MESHSET_SET, newSet);  MB_CHK_ERR(rval);
    rval = mb-> tag_set_data(mtag, &newSet, 1 , &valnew);MB_CHK_ERR(rval);
    Range inNewSet = intersect(ents, newQuads);
    rval = mb->add_entities(newSet, inNewSet);MB_CHK_ERR(rval);
    rval = mb->remove_entities(sets[i], inNewSet);
    cout << "set " << i << " val: " << vals[i] << " ents.size()=" << ents.size() << " inNewSet:" << inNewSet.size() << "\n";
  }
  
  mb->write_file("newSets.exo");
  return 0;
}
