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
#include "moab/MergeMesh.hpp"
#include "MBTagConventions.hpp"

#include <iostream>

using namespace moab;
using namespace std;

#ifndef MESH_DIR
#define MESH_DIR "."
#endif

string test_file_name("tubes.exo");
string test_file_name2("blo.exo");

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
  if (argc > 2) {
    // User has input a mesh file
    test_file_name2 = argv[2];
  }
  EntityHandle file_set1, file_set2;
  rval = mb->create_meshset(MESHSET_SET, file_set1);
  rval = mb->create_meshset(MESHSET_SET, file_set2);

  // Load a file
  rval = mb->load_file(test_file_name.c_str(), &file_set1);MB_CHK_ERR(rval);
  rval = mb->load_file(test_file_name2.c_str(), &file_set2);MB_CHK_ERR(rval);

  rval = mb->tag_get_handle(tag_nms[0], 1, MB_TYPE_INTEGER, mtag);MB_CHK_ERR(rval);
  rval = mb->get_entities_by_type_and_tag(0, MBENTITYSET, &mtag, NULL, 1, sets);MB_CHK_ERR(rval);
  vector<int>  vals;
  vals.resize(sets.size());
  rval = mb->tag_get_data(mtag, sets, &vals[0]);
 
  cout<< " vals: " ; 
  for (int i=0; i<(int)vals.size(); i++)
     cout << vals[i] << " " ;
  
  MergeMesh mm(mb);
  
  Range ents;
  EntityHandle set1=sets[1] , set2=sets[9]; // 20 with 100
  // match sets with some values  20 with 100; 40 to 110, 60 to 120, 80 to 130
  cout << endl;
  rval = mb->get_entities_by_handle(set1, ents, true);MB_CHK_ERR(rval);
  rval = mb->get_entities_by_handle(set2, ents, true);MB_CHK_ERR(rval);
  rval = mb->get_entities_by_handle(sets[3], ents, true);MB_CHK_ERR(rval);
  rval = mb->get_entities_by_handle(sets[5], ents, true);MB_CHK_ERR(rval);
  rval = mb->get_entities_by_handle(sets[7], ents, true);MB_CHK_ERR(rval);
  rval = mb->get_entities_by_handle(sets[10], ents, true);MB_CHK_ERR(rval);
  rval = mb->get_entities_by_handle(sets[11], ents, true);MB_CHK_ERR(rval);
  rval = mb->get_entities_by_handle(sets[12], ents, true);MB_CHK_ERR(rval);
  std::cout<<" ents.size():" << ents.size() << "\n";
  Range verts;
  EntityHandle setGoble;
  rval = mb->create_meshset(MESHSET_SET, setGoble);MB_CHK_ERR(rval);
  rval = mb->add_entities(setGoble, ents);
  rval = mb->get_connectivity(ents, verts); MB_CHK_ERR(rval);
  double merge_tol = 0.0001;
  cout << " start merge sets with values: " << vals[1] << " " << vals[9] << " \n";
  cout << " verts size: " << verts.size() << endl;
  mm.merge_all(setGoble, merge_tol);
  verts.clear();
  rval = mb->get_connectivity(ents, verts); MB_CHK_ERR(rval);
  cout << " after merge: verts size: " << verts.size() << endl;
  
  mb->clear_meshset(&setGoble, 1);

  ents.clear();

  mb->write_file("merged.exo");

  return 0;
}
