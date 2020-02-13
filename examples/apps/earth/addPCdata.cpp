/** @example addPCdata.cpp  Add point cloud data
 * this tool will take an existing h5m fine atm mesh file and add data from an h5m  type file with point cloud mesh
 * will support mainly showing the data with Visit
 *
 * example of usage:
 * ./mbaddpcdata -i wholeFineATM.h5m -s wholeLND_proj01.h5m -o atm_a2l.h5m -v a2lTbot_proj
 *
 * it should work also for pentagon file style data
 * ./mbaddpcdata -i <MOABsource>/MeshFiles/unittest/penta3d.h5m -s wholeLND_proj01.h5m -o atm_a2l.h5m -v a2lTbot_proj -p 1
 *
 * Basically, will output a new h5m file (atm_a2l.h5m), which has an extra tag, corresponding to the variable
 *   a2lTbot_proj from the file wholeLND_proj01.h5m; matching is based on the global ids between what we think is the order
 *   on the original file (wholeFineATM.h5m) and the order of surfdata_ne11np4_simyr1850_c160614.nc
 *
 *  file  wholeFineATM.h5m is obtained from a coupled run in e3sm, with the ne 11, np 4,
 *
 */
#include "moab/ProgOptions.hpp"
#include "moab/Core.hpp"

#include <math.h>
#include <sstream>

using namespace moab;

int main(int argc, char* argv[])
{

  ProgOptions opts;

  std::string inputfile, outfile("out.h5m"), sourcefile, variable_name;
  int dual_mesh = 0;
  opts.addOpt<std::string>("input,i", "input mesh filename", &inputfile);
  opts.addOpt<std::string>("source,s", "h5m file aligned with the mesh input file", &sourcefile);
  opts.addOpt<std::string>("output,o", "output mesh filename", &outfile);

  opts.addOpt<std::string>("var,v", "variable to extract and add to output file", &variable_name);
  opts.addOpt<int> ("pentagon,p", "switch for dual mesh ", &dual_mesh);

  opts.parseCommandLine(argc, argv);

  std::cout << "input file: " << inputfile << "\n";
  std::cout << "source file: " << sourcefile << "\n";
  std::cout << "variable to copy: " << variable_name << "\n";
  std::cout << "output file: " << outfile << "\n";


  if (inputfile.empty())
  {
    opts.printHelp();
    return 0;
  }
  ErrorCode rval;
  Core *mb = new Core();

  rval = mb->load_file(inputfile.c_str()); MB_CHK_SET_ERR(rval, "can't load input file");

  Core *mb2 = new Core();
  rval = mb2->load_file(sourcefile.c_str()); MB_CHK_SET_ERR(rval, "can't load source file");

  Tag sourceTag;
  rval = mb2-> tag_get_handle(variable_name.c_str(), sourceTag); MB_CHK_SET_ERR(rval, "can't get tag from file");

  int sizeTag =0;
  rval = mb2->tag_get_length(sourceTag, sizeTag);  MB_CHK_SET_ERR(rval, "can't get size of tag ");

  DataType type;
  rval = mb2->tag_get_data_type(sourceTag, type);  MB_CHK_SET_ERR(rval, "can't get type of tag ");

  int sizeInBytes = 0;
  rval = mb2->tag_get_bytes(sourceTag, sizeInBytes);  MB_CHK_SET_ERR(rval, "can't get size in bytes of tag ");

  Tag newTag;
  /*
   *  ErrorCode tag_get_handle( const char* name,
                                    int size,
                                    DataType type,
                                    Tag& tag_handle,
                                    unsigned flags = 0,
                                    const void* default_value = 0 )
   */
  void * defVal;
  int defInt = -100000;
  double defDouble = -1.e10;
  if (type==MB_TYPE_INTEGER)
  {
    defVal= (void*)(&defInt);
  }
  else if (type == MB_TYPE_DOUBLE)
  {
    defVal= (void*)(&defDouble);
  }

  rval = mb-> tag_get_handle(variable_name.c_str(), sizeTag, type, newTag, MB_TAG_DENSE|MB_TAG_CREAT, defVal);MB_CHK_SET_ERR(rval, "can't create new tag ");

  // get vertices on ini mesh; get global id on ini mesh
  // get global id on source mesh
  Tag gid;
  Tag gid2;
  rval = mb->tag_get_handle("GLOBAL_ID", gid);MB_CHK_SET_ERR(rval, "can't get GLOBAL_ID tag on ini mesh ");
  rval = mb2->tag_get_handle("GLOBAL_ID", gid2);MB_CHK_SET_ERR(rval, "can't get GLOBAL_ID tag on source mesh ");

  // get vertices on ini mesh; build
  Range iniVerts;
  rval = mb->get_entities_by_dimension(0, 0, iniVerts);MB_CHK_SET_ERR(rval, "can't get verts on initial mesh ");
  if (0!=dual_mesh)
  {
    // verts will be polygons
    rval = mb->get_entities_by_dimension(0, 2, iniVerts);MB_CHK_SET_ERR(rval, "can't get polygons on initial mesh ");
  }
  std::vector<int> gids;
  gids.resize(iniVerts.size());
  rval = mb->tag_get_data(gid, iniVerts, &(gids[0]));MB_CHK_SET_ERR(rval, "can't get gid on initial verts ");
  // build now the map
  std::map<int, EntityHandle> fromGidToEh;
  int i=0;
  for (Range::iterator vit=iniVerts.begin(); vit!=iniVerts.end(); vit++, i++)
  {
    fromGidToEh[gids[i]]=*vit;
  }
  // now get the source verts, and tags, and set it on new mesh

  char * valTag = new char[sizeInBytes];

  std::cout << " size of tag in bytes:" << sizeInBytes << "\n";
  Range sourceVerts;
  rval = mb2->get_entities_by_dimension(0, 0, sourceVerts);MB_CHK_SET_ERR(rval, "can't get verts on source mesh ");
  for (Range::iterator sit=sourceVerts.begin(); sit!=sourceVerts.end(); sit++)
  {
    int globalId =0;
    EntityHandle sourceHandle = *sit;
    rval = mb2->tag_get_data(gid2, &sourceHandle, 1, &globalId); MB_CHK_SET_ERR(rval, "can't get id on source mesh ");
    // iniVert could be a polygon, actually
    EntityHandle iniVert = fromGidToEh[globalId];
    // find the value of source tag
    rval = mb2->tag_get_data(sourceTag, &sourceHandle, 1, (void*)valTag); MB_CHK_SET_ERR(rval, "can't get value on source tag ");

    rval = mb->tag_set_data(newTag, &iniVert, 1, (void*)valTag); MB_CHK_SET_ERR(rval, "can't set value on initial mesh, new tag ");

  }

  // save file
  rval = mb->write_file(outfile.c_str()); MB_CHK_SET_ERR(rval, "can't write file");

  delete [] valTag;
  delete mb;
  delete mb2;

  return 0;
}

