/** @example DeformMeshRemap.cpp
 * Description: Account for mesh deformation of a solid due to structural mechanics\n
 * In this example there are two meshes, a "master" and "slave" mesh.  In the master mesh,
 * the solid material is deformed, to mimic what happens when a solid heats up and deforms.
 * The fluid mesh is smoothed to account for those deformations, and tags on the fluid are
 * remapped to those new positions.  Then mesh positions and state variables are transferred
 * to the slave mesh, mimicing another mesh used by some other physics.
 *
 * To run: ./DeformMeshRemap [<master_meshfile> <slave_meshfile>]\n
 * (default values can run if users don't specify the mesh files)
 */

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/LloydSmoother.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/BoundBox.hpp"
#include "moab/SpatialLocator.hpp"
#include "MBTagConventions.hpp"
#include "DataCoupler.hpp"

#ifdef USE_MPI
#  include "moab/ParallelComm.hpp"
#endif

#include <iostream>
#include <set>
#include <sstream>
#include <assert.h>

using namespace moab;
using namespace std;

#ifndef MESH_DIR
#define MESH_DIR "."
#endif

ErrorCode read_file(string &fname, EntityHandle &seth, 
                    Range &solids, Range &solid_elems, Range &fluids, Range &fluid_elems);
void deform_func(const BoundBox &bbox, double *xold, double *xnew);
ErrorCode deform_master(Range &fluid_elems, Range &solid_elems, Tag &xnew);
ErrorCode smooth_master(int dim, Tag xnew, EntityHandle &master, Range &fluids);
ErrorCode write_to_coords(Range &elems, Tag tagh);

const int SOLID_SETNO = 100, FLUID_SETNO = 200;

Interface *mb;
#define RR(a) if (MB_SUCCESS != rval) {cout << a << endl; return MB_FAILURE;}

const bool debug = true;

class DeformMeshRemap 
{
public:

    //! enumerator for solid/fluid, master/slave
  enum {MASTER=0, SLAVE, SOLID, FLUID};
  
    //! constructor
    //! if master is NULL, the MOAB part is run in serial; 
    //! if slave is NULL but the master isn't, the slave is copied from the master
    //! Create communicators using moab::ParallelComm::get_pcomm
  DeformMeshRemap(Interface *impl, ParallelComm *master = NULL, ParallelComm *slave = NULL);
  
    //! destructor
  ~DeformMeshRemap();

    //! execute the deformed mesh process
  ErrorCode execute();
  
    //! add a set number
  ErrorCode add_set_no(int fluid_or_solid, int set_no);
  
    //! remove a set number
  ErrorCode remove_set_no(int fluid_or_solid, int set_no);
  
    //! get the set numbers
  ErrorCode get_set_nos(int fluid_or_solid, std::set<int> &set_nos) const;

    //! get the xNew tag handle
  inline Tag x_new() const {return xNew;}

    //! get the tag name
  std::string x_new_name() const {return xNewName;}
  
    //! set the tag name
  void x_new_name(const std::string &name) {xNewName = name;}

    //! get/set the file name
  std::string get_file_name(int m_or_s) const;
  
    //! get/set the file name
  void set_file_name(int m_or_s, const std::string &name);
  
private:
    //! apply a known deformation to the solid elements, putting the results in the xNew tag; also
    //! write current coordinates to the xNew tag for fluid elements
  ErrorCode deform_master(Range &fluid_elems, Range &solid_elems, const char *tag_name = NULL);

    //! read a file and establish proper ranges
  ErrorCode read_file(int m_or_s, string &fname, EntityHandle &seth);

    //! write the input tag to the coordinates for the vertices in the input elems
  ErrorCode write_to_coords(Range &elems, Tag tagh);

    //! write the tag to the vertices, then save to the specified file
  ErrorCode write_and_save(Range &ents, EntityHandle seth, Tag tagh, const char *filename);

    //! moab interface
  Interface *mbImpl;

#ifdef USE_MPI
    //! ParallelComm for master, slave meshes
  ParallelComm *pcMaster, *pcSlave;
#endif
  
    //! material set numbers for fluid materials
  std::set<int> fluidSetNos;

    //! material set numbers for solid materials
  std::set<int> solidSetNos;

    //! sets defining master/slave meshes
  EntityHandle masterSet, slaveSet;

    //! sets in master/slave meshes
  Range fluidSets[2], solidSets[2];
  
    //! elements in master/slave meshes
  Range fluidElems[2], solidElems[2];
  
    //! filenames for master/slave meshes
  std::string masterFileName, slaveFileName;

    //! tag used for new positions
  Tag xNew;
  
    //! tag name used for new positions
  std::string xNewName;
};

  //! add a set number
inline ErrorCode DeformMeshRemap::add_set_no(int f_or_s, int set_no) 
{
  std::set<int> *this_set;
  switch (f_or_s) {
    case FLUID:
        this_set = &fluidSetNos; break;
    case SOLID:
        this_set = &solidSetNos; break;
    default:
        assert(false && "f_or_s should be FLUID or SOLID.");
        return MB_FAILURE;
  }

  this_set->insert(set_no);
  
  return MB_SUCCESS;
}
  
  //! remove a set number
inline ErrorCode DeformMeshRemap::remove_set_no(int f_or_s, int set_no) 
{
  std::set<int> *this_set;
  switch (f_or_s) {
    case FLUID:
        this_set = &fluidSetNos; break;
    case SOLID:
        this_set = &solidSetNos; break;
    default:
        assert(false && "f_or_s should be FLUID or SOLID.");
        return MB_FAILURE;
  }
  std::set<int>::iterator sit = this_set->find(set_no);
  if (sit != this_set->end()) {
    this_set->erase(*sit);
    return MB_SUCCESS;
  }

  return MB_FAILURE;
}
  
  //! get the set numbers
inline ErrorCode DeformMeshRemap::get_set_nos(int f_or_s, std::set<int> &set_nos) const
{
  const std::set<int> *this_set;
  switch (f_or_s) {
    case FLUID:
        this_set = &fluidSetNos; break;
    case SOLID:
        this_set = &solidSetNos; break;
    default:
        assert(false && "f_or_s should be FLUID or SOLID.");
        return MB_FAILURE;
  }

  set_nos = *this_set;
  
  return MB_SUCCESS;
}

ErrorCode DeformMeshRemap::execute() 
{
    // read master/slave files and get fluid/solid material sets
  ErrorCode rval = read_file(MASTER, masterFileName, masterSet);
  if (MB_SUCCESS != rval) return rval;
  
  rval = read_file(SLAVE, slaveFileName, slaveSet);
  if (MB_SUCCESS != rval) return rval;

  Range src_elems = solidElems[MASTER];
  src_elems.merge(fluidElems[MASTER]);
    // locate slave vertices in master, orig coords; do this with a data coupler, so you can
    // later interpolate
  Range tgt_verts, tmp_range = solidElems[SLAVE];
  tmp_range.merge(fluidElems[SLAVE]);
  rval = mbImpl->get_adjacencies(tmp_range, 0, false, tgt_verts, Interface::UNION);
  RR("Failed to get target verts.");
  

    // initialize data coupler on source elements
  DataCoupler dc_master(mbImpl, NULL, src_elems, 0);
  
    // locate slave vertices, caching results in dc
  rval = dc_master.locate_points(tgt_verts); RR("Point location of tgt verts failed.");
  int num_located = dc_master.spatial_locator()->local_num_located();
  if (num_located != (int)tgt_verts.size()) {
    rval = MB_FAILURE;
    std::cout << "Only " << num_located << " out of " << tgt_verts.size() << " target points successfully located." << std::endl;
    return rval;
  }

    // deform the master's solid mesh, put results in a new tag
  rval = deform_master(fluidElems[MASTER], solidElems[MASTER], "xnew"); RR("");

  { // to isolate the lloyd smoother & delete when done

      // smooth the master mesh
    LloydSmoother ll(mbImpl, NULL, fluidElems[MASTER], xNew);
    rval = ll.perform_smooth();
    RR("Failed in lloyd smoothing.");
    cout << "Lloyd smoothing required " << ll.num_its() << " iterations." << endl;
  }
  
    // map new locations to slave
    // interpolate xNew to slave points
  rval = dc_master.interpolate((int)DataCoupler::VOLUME, "xnew"); RR("Failed to interpolate target solution.");

    // transfer xNew to coords, for master and slave
  rval = write_to_coords(fluidElems[MASTER], xNew); RR("Failed writing tag to master fluid verts.");
  rval = write_to_coords(tgt_verts, xNew); RR("Failed writing tag to slave verts.");

  if (debug) {
    std::string str;
#ifdef USE_MPI
    if (pcMaster && pcMaster->size() > 1) 
      str = "PARALLEL=WRITE_PART";
#endif
    rval = mbImpl->write_file("smoothed_master.vtk", NULL, str.c_str(), &masterSet, 1);
#ifdef USE_MPI
    str.clear();
    if (pcSlave && pcSlave->size() > 1) 
      str = "PARALLEL=WRITE_PART";
#endif
    rval = mbImpl->write_file("slave_interp.vtk", NULL, str.c_str(), &slaveSet, 1);
  }

  return MB_SUCCESS;
}

std::string DeformMeshRemap::get_file_name(int m_or_s) const
{
  switch (m_or_s) {
    case MASTER:
        return masterFileName;
    case SLAVE:
        return slaveFileName;
    default:
        assert(false && "m_or_s should be MASTER or SLAVE.");
        return std::string();
  }
}
  
void DeformMeshRemap::set_file_name(int m_or_s, const std::string &name) 
{
  switch (m_or_s) {
    case MASTER:
        masterFileName = name; break;
    case SLAVE:
        slaveFileName = name; break;
    default:
        assert(false && "m_or_s should be MASTER or SLAVE.");
  }
}

DeformMeshRemap::DeformMeshRemap(Interface *impl, ParallelComm *master, ParallelComm *slave)  
        : mbImpl(impl), pcMaster(master), pcSlave(slave), masterSet(0), slaveSet(0), xNew(0), xNewName("xnew") 
{
  if (!pcSlave && pcMaster)
    pcSlave = pcMaster;
}
  
DeformMeshRemap::~DeformMeshRemap() 
{
    // delete the tag
  mbImpl->tag_delete(xNew);
}

int main(int argc, char **argv) {

  ErrorCode rval;

  ProgOptions po("Deformed mesh options");
  po.addOpt<std::string> ("master,m", "Specify the master meshfile name" );
  po.addOpt<std::string> ("slave,s", "Specify the slave meshfile name" );
  po.parseCommandLine(argc, argv);
  std::string foo;
  string masterf, slavef;
  if(!po.getOpt("master", &masterf))
    masterf = string(MESH_DIR) + string("/rodquad.g");
  if(!po.getOpt("slave", &slavef))
    slavef = string(MESH_DIR) + string("/rodtri.g");

  mb = new Core();

  DeformMeshRemap *dfr;
#ifdef USE_MPI
  ParallelComm *pc = new ParallelComm(mb, MPI_COMM_WORLD);
  dfr = new DeformMeshRemap(mb, pc);
#else  
  dfr = new DeformMeshRemap(mb);
#endif
  dfr->set_file_name(DeformMeshRemap::MASTER, masterf);
  dfr->set_file_name(DeformMeshRemap::SLAVE, slavef);
  rval = dfr->add_set_no(DeformMeshRemap::SOLID, SOLID_SETNO); RR("Failed to add solid set no.");
  rval = dfr->add_set_no(DeformMeshRemap::FLUID, FLUID_SETNO); RR("Failed to add fluid set no.");
  
  rval = dfr->execute();
  
  delete dfr;
  delete mb;
  
  return rval;
}

ErrorCode DeformMeshRemap::write_and_save(Range &ents, EntityHandle seth, Tag tagh, const char *filename) 
{
  ErrorCode rval = write_to_coords(ents, tagh); RR("");
  rval = mbImpl->write_file(filename, NULL, NULL, &seth, 1); RR("");
  return rval;
}
  
ErrorCode DeformMeshRemap::write_to_coords(Range &elems, Tag tagh) 
{
    // write the tag to coordinates
  Range verts;
  ErrorCode rval = mbImpl->get_adjacencies(elems, 0, false, verts, Interface::UNION);
  RR("Failed to get adj vertices.");
  std::vector<double> coords(3*verts.size());
  rval = mbImpl->tag_get_data(tagh, verts, &coords[0]);
  RR("Failed to get tag data.");
  rval = mbImpl->set_coords(verts, &coords[0]);
  RR("Failed to set coordinates.");
  return MB_SUCCESS;
}

void deform_func(const BoundBox &bbox, double *xold, double *xnew) 
{
/*  Deformation function based on max delx and dely at top of rod
    const double RODWIDTH = 0.2, RODHEIGHT = 0.5;
    // function: origin is at middle base of rod, and is .5 high
    // top of rod is (0,.55) on left and (.2,.6) on right
  double delx = 0.5*RODWIDTH;
  
  double xfrac = (xold[0] + .5*RODWIDTH)/RODWIDTH, yfrac = xold[1]/RODHEIGHT;
  xnew[0] = xold[0] + yfrac * delx;
  xnew[1] = xold[1] + yfrac * (1.0 + xfrac) * 0.05;
*/

/* Deformation function based on fraction of bounding box dimension in each direction */
  double frac = 0.01; // taken from approximate relative deformation from LLNL Diablo of XX09 assys
  CartVect *xo = reinterpret_cast<CartVect*>(xold), *xn = reinterpret_cast<CartVect*>(xnew);
  CartVect disp = frac * (*xo - bbox.bMin);
  *xn = *xo + disp;
}
  
ErrorCode DeformMeshRemap::deform_master(Range &fluid_elems, Range &solid_elems, const char *tag_name) 
{
    // deform elements with an analytic function

    // create the tag
  ErrorCode rval = mbImpl->tag_get_handle((tag_name ? tag_name : ""), 3, MB_TYPE_DOUBLE, xNew, MB_TAG_CREAT|MB_TAG_DENSE);
  RR("Failed to create xnew tag.");
  
    // get all the vertices and coords in the fluid, set xnew to them
  Range verts;
  rval = mbImpl->get_adjacencies(fluid_elems, 0, false, verts, Interface::UNION);
  RR("Failed to get vertices.");
  std::vector<double> coords(3*verts.size(), 0.0);
  rval = mbImpl->get_coords(verts, &coords[0]);
  RR("Failed to get vertex coords.");
  rval = mbImpl->tag_set_data(xNew, verts, &coords[0]);
  RR("Failed to set xnew tag on fluid verts.");

    // get the bounding box of the solid mesh
  BoundBox bbox;
  bbox.update(*mbImpl, solid_elems);
  
    // get all the vertices and coords in the solid
  verts.clear();
  rval = mbImpl->get_adjacencies(solid_elems, 0, false, verts, Interface::UNION);
  RR("Failed to get vertices.");
  coords.resize(3*verts.size(), 0.0);
  rval = mbImpl->get_coords(verts, &coords[0]);
  RR("Failed to get vertex coords.");
  unsigned int num_verts = verts.size();
  for (unsigned int i = 0; i < num_verts; i++)
    deform_func(bbox, &coords[3*i], &coords[3*i]);
    
    // set the new tag to those coords
  rval = mbImpl->tag_set_data(xNew, verts, &coords[0]);
  RR("Failed to set tag data.");
  
  return MB_SUCCESS;
}

ErrorCode DeformMeshRemap::read_file(int m_or_s, string &fname, EntityHandle &seth)
{
    // create meshset
  ErrorCode rval = mbImpl->create_meshset(0, seth);
  RR("Couldn't create master/slave set.");
  ostringstream ostr;
#ifdef USE_MPI
  ParallelComm *pc = (m_or_s == MASTER ? pcMaster : pcSlave);
  if (pc && pc->size() > 1) {
    if (debug) ostr << "DEBUG_IO=1;CPUTIME;";
    ostr << "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;"
         << "PARALLEL_GHOSTS=2.0.1;PARALLEL_COMM=" << pc->get_id();
  }
#endif  
  rval = mbImpl->load_file(fname.c_str(), &seth, ostr.str().c_str());
  RR("Couldn't load master/slave mesh.");

    // get material sets for solid/fluid
  Tag tagh;
  rval = mbImpl->tag_get_handle(MATERIAL_SET_TAG_NAME, tagh); RR("Couldn't get material set tag name.");
  for (std::set<int>::iterator sit = solidSetNos.begin(); sit != solidSetNos.end(); sit++) {
    Range sets;
    int set_no = *sit;
    const void *setno_ptr = &set_no;
    rval = mbImpl->get_entities_by_type_and_tag(seth, MBENTITYSET, &tagh, &setno_ptr, 1, sets);
    if (sets.empty()) rval = MB_FAILURE;
    RR("Couldn't get any solid sets.");
    solidSets[m_or_s].merge(sets);
  }

    // get solid entities, and dimension
  Range tmp_range;
  for (Range::iterator rit = solidSets[m_or_s].begin(); rit != solidSets[m_or_s].end(); rit++) {
    rval = mbImpl->get_entities_by_handle(*rit, tmp_range, true);
    RR("Failed to get entities in solid.");
  }
  int dim = mbImpl->dimension_from_handle(*tmp_range.rbegin());
  assert(dim > 0 && dim < 4);
  
  solidElems[m_or_s] = tmp_range.subset_by_dimension(dim);

  for (std::set<int>::iterator sit = fluidSetNos.begin(); sit != fluidSetNos.end(); sit++) {
    Range sets;
    int set_no = *sit;
    const void *setno_ptr = &set_no;
    rval = mbImpl->get_entities_by_type_and_tag(seth, MBENTITYSET, &tagh, &setno_ptr, 1, sets);
    if (sets.empty()) rval = MB_FAILURE;
    RR("Couldn't get any fluid sets.");
    fluidSets[m_or_s].merge(sets);
  }

    // get fluid entities, and dimension
  tmp_range.clear();
  for (Range::iterator rit = fluidSets[m_or_s].begin(); rit != fluidSets[m_or_s].end(); rit++) {
    rval = mbImpl->get_entities_by_handle(*rit, tmp_range, true);
    RR("Failed to get entities in fluid.");
  }
  
  fluidElems[m_or_s] = tmp_range.subset_by_dimension(dim);
  
  return MB_SUCCESS;
}
