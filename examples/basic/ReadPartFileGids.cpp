/** @example ReadVarTag.cpp
 *
 * read variable sparse tag on a set
 */


#include "moab/Core.hpp"
#include <iostream>
#include <fstream>

using namespace moab;
using namespace std;

#ifndef MESH_DIR
#define MESH_DIR "."
#endif

// Note: change the file name below to test a trivial "No such file or directory" error
string test_file_name = string(MESH_DIR) + string("/3k-tri-sphere.vtk");
string part_file_name;
string gids_file;
int nparts;
int main(int argc, char **argv)
{
  // Get MOAB instance
  Interface* mb = new (std::nothrow) Core;
  if (NULL == mb)
    return 1;

  // Need option handling here for input filename
  if (argc > 5) {
    // User has input a mesh file
    test_file_name = argv[1];
    part_file_name = argv[2];
    nparts = atoi(argv[3]);
    gids_file = argv[4];
  }
  else
  {
    cerr << " usage is " << argv[0] << " <input file> <part file> <#parts> <gids_file> <output file> \n";
    exit(1);
  }
  ifstream inFile;
  inFile.open(part_file_name.c_str());
  if (!inFile) {
      cerr << "Unable to open file " << part_file_name << "\n";
      exit(1);   // call system to stop
  }
  // Load the mesh from file
  ErrorCode rval = mb->load_mesh(test_file_name.c_str());MB_CHK_ERR(rval);

  // Get sets entities, by type
  Range sets;
  rval = mb->get_entities_by_type(0, MBENTITYSET, sets);MB_CHK_ERR(rval);

  // Output the number of sets    
  cout << "Number of sets is " << sets.size() << endl;
  // remove the sets that have a PARALLEL_PARTITION tag
  
  Tag tag;
  rval = mb->tag_get_handle( "PARALLEL_PARTITION" , tag);MB_CHK_ERR(rval);

  int i = 0;
  int num_deleted_sets = 0;
  for (Range::iterator it=sets.begin(); it!=sets.end(); it++, i++)
  {
    EntityHandle eh = *it;
    // cout << " set :" << mb->id_from_handle(eh) <<"\n";
    int val=-1;
    rval = mb->tag_get_data(tag, &eh, 1, &val); 
    if (val != -1)
    { 
      num_deleted_sets++;
      rval = mb->delete_entities(&eh, 1); // delete the set, we will have a new partition soon
    }
  }
  if (num_deleted_sets)
    cout << "delete " << num_deleted_sets << " existing  partition sets, and create new ones \n";
  Range cells; // get them by dimension 2!
  rval = mb->get_entities_by_dimension(0, 2, cells); MB_CHK_ERR(rval);
  EntityHandle * psets = new EntityHandle[nparts];
  for (int i=0; i<nparts; i++)
  {
    rval = mb->create_meshset(MESHSET_SET, psets[i]);MB_CHK_ERR(rval);
    rval = mb->tag_set_data(tag, &(psets[i]),1, &i);MB_CHK_ERR(rval);
  }

  // get the gids of cells from the file
  ifstream gidFile;
  gidFile.open(gids_file.c_str());
  if (!gidFile) {
      cerr << "Unable to open file " << gids_file << "\n";
      exit(1);   // call system to stop
  }

  std::map<int, moab::EntityHandle> gidMap;
  moab::Tag gidTag;
  rval = mb->tag_get_handle("GLOBAL_ID", gidTag);MB_CHK_ERR(rval);
  std::vector<int> gids; gids.resize(cells.size());
  rval = mb->tag_get_data(gidTag, cells, &gids[0]); MB_CHK_ERR(rval);

  i=0;
  for (Range::iterator it = cells.begin(); it!=cells.end(); it++)
  {
    gidMap[gids[i++]] = *it;
  }
  for (int j=0; j<(int)cells.size();j++)
  {
    int part;
    int gid;

    inFile >> part;
    gidFile >> gid;
    EntityHandle eh = gidMap[gid];
    rval = mb->add_entities(psets[part], &eh, 1); MB_CHK_ERR(rval);
  }
  mb->write_file(argv[5]);

  delete mb;

  return 0;
}
