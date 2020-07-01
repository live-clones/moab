/** @example QuadTriConv.cpp \n
 * \brief Merge vertices in 2d cell, if they are repeated \n
 * <b>To run</b>: QuadTriConv  input_file output_file \n
 *
 * In this example, a mesh that has vertices repeated due to south or north pole, for 
 *  example in an rll mesh, need to be converted to triangles, preserving GLOBAL_ID tag
 */

#include <iostream>
#include <vector>
//#include <string>

// Include header for MOAB instance and tag conventions for
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"

using namespace moab;
using namespace std;


int main(int argc, char **argv)
{
  // Get MOAB instance
  Interface* mb = new (std::nothrow) Core;
  if (NULL == mb)
    return 1;

  std::string filename,  outfile;
  outfile = string("out.h5m");
  if (argc == 1)
    return 0;
  if (argc > 1)
    filename = string(argv[1]);
  if (argc > 2)
    outfile = string(argv[2]);

  // This file is in the mesh files directory
  ErrorCode rval = mb->load_file(filename.c_str());MB_CHK_SET_ERR(rval, "Failed to read");

  // get all cells of dimension 2; 
  Range cells;
  rval = mb->get_entities_by_dimension(0, 2, cells);MB_CHK_SET_ERR(rval, "Failed to get cells");

  cout << " number of cells : " << cells.size() << "\n";

  Tag gid; // global id tag
  rval = mb->tag_get_handle("GLOBAL_ID", gid);MB_CHK_SET_ERR(rval, "Failed to get Global ID tag");

  Range modifiedCells; // will be deleted at the end; keep the gid

  for (Range::iterator cit=cells.begin(); cit != cells.end(); cit++)
  {
    EntityHandle cell = *cit;
    const EntityHandle * connec = NULL;
    int num_verts = 0;
    rval = mb->get_connectivity(cell, connec, num_verts); MB_CHK_SET_ERR(rval, "Failed to get connectivity");

    vector<EntityHandle> newConnec;
    newConnec.push_back(connec[0]); // at least one vertex
    int index = 0;
    int new_size = 1;
    while (index<num_verts-2)
    {
      int next_index = (index+1);
      if (connec[next_index] != newConnec[new_size-1])
      {
        newConnec.push_back(connec[next_index]);
        new_size++;
      }
      index++;
    }
    // add the last one only if different from previous and first node
    if( ( connec[num_verts-1]!=connec[num_verts-2])  &&
        ( connec[num_verts-1]!= connec[0]) )
    {
      newConnec.push_back(connec[num_verts-1]);
      new_size++;
    }
    if (new_size < num_verts)
    {
      //cout << "new cell from " << cell << " has only " << new_size << " vertices \n";
      modifiedCells.insert(cell);
      // create a new cell with type triangle, quad or polygon
      EntityType type = MBTRI;
      if (new_size==3)
        type = MBTRI;
      else if (new_size == 4)
        type = MBQUAD;
      else if (new_size > 4)
        type = MBPOLYGON;

      // create new cell
      EntityHandle newCell;
      rval = mb->create_element(type, &newConnec[0], new_size, newCell); MB_CHK_SET_ERR(rval, "Failed to create new cell");
      // set the old id to the new element
      int id;
      rval = mb->tag_get_data(gid, &cell, 1, &id);  MB_CHK_SET_ERR(rval, "Failed to get global id");
      rval = mb->tag_set_data(gid, &newCell, 1, &id);  MB_CHK_SET_ERR(rval, "Failed to set global id on new cell");
    }
  }


  mb->delete_entities(modifiedCells);
  // in case global ids are not set for vertices, set them, in order;
  // they are needed for migrating later on
  // 
  Range verts;
  rval = mb->get_entities_by_dimension(0, 0, verts);MB_CHK_SET_ERR(rval, "Failed to get vertices");

  vector<int> gids;
  gids.resize(verts.size());
  rval = mb->tag_get_data(gid, verts, &gids[0]);  MB_CHK_SET_ERR(rval, "Failed to get gids on vertices");
  if (gids[0]<=0)
  {
    // gids were never set, assign them
    for (size_t i=1; i<=verts.size(); i++)
      gids[i-1] = (int)i;
    rval = mb->tag_set_data(gid, verts, &gids[0]);  MB_CHK_SET_ERR(rval, "Failed to set gids on vertices");
  }
  rval = mb->write_file(outfile.c_str()); MB_CHK_SET_ERR(rval, "Failed to write new file");

  cout << " wrote file " << outfile << " with " << modifiedCells.size() << " modified cells\n";

  delete mb;
  return 0;
}
