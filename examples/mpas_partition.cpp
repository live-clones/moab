#include "moab/Core.hpp"
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace moab;
using namespace std;

// Note: change the file name below to test a trivial "No such file or directory" error
string test_file_name = string("x4.163842.grid.nc");
string partition_file_name=  string("x4.163842.graph.info.part.48");

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

  // Load the mesh from mpas  file
  ErrorCode rval = mb->load_mesh(test_file_name.c_str());MB_CHK_ERR(rval);
  // how many faces in the file ?
  Range faces;
  rval = mb->get_entities_by_dimension(0, 2, faces); MB_CHK_ERR(rval);
  cout << " MPAS model has " << faces.size() << " polygons\n";
 
  if (argc > 2) {
    partition_file_name= argv[2];
  }

  // read the partition file
  ifstream partfile;
  partfile.open(partition_file_name.c_str());
  string line;
  std::vector<int> parts;
  parts.resize(faces.size());
  int i=0;
  if (partfile.is_open())
  {
    while ( getline (partfile,line) )
    {
      //cout << line << '\n';
      parts[i++] = atoi(line.c_str());
      if (i>(int)faces.size())
      {
         cout << " too many lines \n. bail out \n";
         return 1;
      }
    }
    partfile.close();
  }
  vector<int>::iterator pmax = max_element(parts.begin(), parts.end());
  vector<int>::iterator pmin = min_element(parts.begin(), parts.end());
  cout << " partitions range: " << *pmin << " " << *pmax << "\n";
  Tag part_set_tag;
  int dum_id = -1;
  rval = mb->tag_get_handle("PARALLEL_PARTITION", 1, MB_TYPE_INTEGER,
                                  part_set_tag, MB_TAG_SPARSE|MB_TAG_CREAT, &dum_id);  MB_CHK_ERR(rval);
  
    // get any sets already with this tag, and clear them
  // remove the parallel partition sets if they exist
  Range tagged_sets;
  rval = mb->get_entities_by_type_and_tag(0, MBENTITYSET, &part_set_tag, NULL, 1,
                                                tagged_sets, Interface::UNION);  MB_CHK_ERR(rval);
  if (!tagged_sets.empty()) {
    rval = mb->clear_meshset(tagged_sets); MB_CHK_ERR(rval);
    rval = mb->tag_delete_data(part_set_tag, tagged_sets); MB_CHK_ERR(rval);
  }
  Tag gid;
  rval = mb->tag_get_handle("GLOBAL_ID", gid);
  int num_sets =  *pmax + 1;
  if (*pmin!=0) 
  {
     cout << " problem reading parts; min is not 0 \n";
     return 1;
  }
  for (i = 0; i < num_sets; i++) {
    EntityHandle new_set;
    rval = mb->create_meshset(MESHSET_SET, new_set); MB_CHK_ERR(rval); 
    tagged_sets.insert(new_set);
  }
  int *dum_ids = new int[num_sets];
    for (i = 0; i < num_sets; i++) dum_ids[i] = i;
  
  rval = mb->tag_set_data(part_set_tag, tagged_sets, dum_ids); MB_CHK_ERR(rval);
  delete [] dum_ids; 
  
  std::vector<int> gids;
  int num_faces = (int)faces.size();
  gids.resize(num_faces);
  rval = mb->tag_get_data(gid, faces, &gids[0]); MB_CHK_ERR(rval);
  
  for (int j=0; j<num_faces; j++)
  {
    int eid = gids[j];
    EntityHandle eh = faces[j];
    int partition = parts[eid-1];
    if (partition <0 || partition >= num_sets)
    {
       cout << " wrong partition \n";
       return 1;
    }
    rval = mb->add_entities(tagged_sets[partition], &eh, 1); MB_CHK_ERR(rval);
  }
  
  // write the new file:
  string out_file("mpas_part.h5m");
  if (argc>3)
     out_file = argv[3];
  
  rval = mb-> write_file(out_file.c_str()); MB_CHK_ERR(rval);
  
  return 0;
}
