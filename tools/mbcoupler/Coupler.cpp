#include "Coupler.hpp"
#include "moab/ParallelComm.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "ElemUtil.hpp"
#include "moab/CN.hpp"
#include "iMesh_extensions.h"
#include "moab/gs.hpp"
#include "moab/TupleList.hpp"
#include "moab/Error.hpp"
#include "iostream"
#include <stdio.h>
#include <algorithm>
#include <sstream>

#include "assert.h"

#define ERROR(a) {if (iBase_SUCCESS != err) std::cerr << a << std::endl;}
#define ERRORR(a,b) {if (iBase_SUCCESS != err) {std::cerr << a << std::endl; return b;}}
#define ERRORMPI(a,b) {if (MPI_SUCCESS != err) {std::cerr << a << std::endl; return b;}}

#define MASTER_PROC 0

namespace moab {

bool debug = false;
int pack_tuples(TupleList* tl, void **ptr);
void unpack_tuples(void *ptr, TupleList** tlp);

Coupler::Coupler(Interface *impl,
                     ParallelComm *pc,
                     Range &local_elems,
                     int coupler_id,
                     bool init_tree,
                     int max_ent_dim)
    : mbImpl(impl), myPc(pc), myId(coupler_id), numIts(3), max_dim(max_ent_dim)
{
  assert(NULL != impl && (pc || !local_elems.empty()));

    // keep track of the local points, at least for now
  myRange = local_elems;

    // now initialize the tree
  if (init_tree) initialize_tree();

    // initialize tuple lists to indicate not initialized
  mappedPts = NULL;
  targetPts = NULL;
  _spectralSource = _spectralTarget = NULL;

  ErrorCode rval = impl->query_interface(mError);
  if (MB_SUCCESS != rval) throw(rval);
}

  /* destructor
   */
Coupler::~Coupler()
{
  // this will clear the cache
  delete (moab::Element::SpectralHex*)_spectralSource;
  delete (moab::Element::SpectralHex*)_spectralTarget;
}


ErrorCode Coupler::initialize_tree()
{
  
  Range local_ents;

    //get entities on the local part
  ErrorCode result = MB_SUCCESS;
  if (myPc) result = myPc->get_part_entities(local_ents, max_dim);
  else local_ents = myRange;
  if (MB_SUCCESS != result || local_ents.empty()) {
    std::cout << "Problems getting source entities" << std::endl;
    return result;
  }

    // build the tree for local processor
  int max_per_leaf = 6;
  for (int i = 0; i < numIts; i++) {
    std::ostringstream str;
    str << "PLANE_SET=0;"
        << "MAX_PER_LEAF=" << max_per_leaf << ";";
    FileOptions opts(str.str().c_str());
    myTree = new AdaptiveKDTree(mbImpl);
    result = myTree->build_tree(local_ents, &localRoot, &opts);
    if (MB_SUCCESS != result) {
      std::cout << "Problems building tree";
      if (numIts != i) {
        delete myTree;
        max_per_leaf *= 2;
        std::cout << "; increasing elements/leaf to " 
                  << max_per_leaf << std::endl;;
      }
      else {
        std::cout << "; exiting" << std::endl;
        return result;
      }
    }
    else
      break; // get out of tree building
  }

    // get the bounding box for local tree
  if (myPc)
    allBoxes.resize(6*myPc->proc_config().proc_size());
  else allBoxes.resize(6);
  unsigned int my_rank = (myPc ? myPc->proc_config().proc_rank() : 0);
  BoundBox box;
  result = myTree->get_bounding_box(box, &localRoot);
  if (MB_SUCCESS != result) return result;
  box.bMin.get(&allBoxes[6*my_rank]);
  box.bMax.get(&allBoxes[6*my_rank+3]);
  
    // now communicate to get all boxes
    // use "in place" option
  if (myPc) {
    int mpi_err = MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                                &allBoxes[0], 6, MPI_DOUBLE, 
                                myPc->proc_config().proc_comm());
    if (MPI_SUCCESS != mpi_err) return MB_FAILURE;
  }

  /*  std::ostringstream blah;
  for(int i=0; i<allBoxes.size(); i++)
  blah << allBoxes[i] << " ";
  std::cout<<blah.str()<<"\n";*/


#ifndef NDEBUG
  double min[3] = {0,0,0}, max[3] = {0,0,0};
  unsigned int dep;
  myTree->get_info(localRoot, min, max, dep);
  std::cout << "Proc " << my_rank << ": box min/max, tree depth = ("
            << min[0] << "," << min[1] << "," << min[2] << "), ("
            << max[0] << "," << max[1] << "," << max[2] << "), "
            << dep << std::endl;
#endif  

  return result;
}

ErrorCode Coupler::initialize_spectral_elements(EntityHandle rootSource, EntityHandle rootTarget,
    bool & specSou, bool & specTar)
{
  /*void * _spectralSource;
    void * _spectralTarget;*/

  moab::Range spectral_sets;
  moab::Tag  sem_tag;
  int sem_dims[3];
  ErrorCode rval = mbImpl->tag_get_handle("SEM_DIMS", 3, moab::MB_TYPE_INTEGER, sem_tag);
  if (moab::MB_SUCCESS != rval)
  {
    std::cout<< "can't find tag, no spectral set\n";
    return MB_SUCCESS; // nothing to do, no spectral elements
  }
  rval = mbImpl->get_entities_by_type_and_tag(rootSource, moab::MBENTITYSET, &sem_tag, NULL, 1, spectral_sets);
  if (moab::MB_SUCCESS != rval || spectral_sets.empty())
  {
    std::cout<< "can't get sem set on source\n";
  }
  else
  {
    moab::EntityHandle firstSemSet=spectral_sets[0];
    rval = mbImpl->tag_get_data(sem_tag, &firstSemSet, 1, (void*)sem_dims);
    if (moab::MB_SUCCESS != rval)
      return MB_FAILURE;

    if (sem_dims[0]!=sem_dims[1] || sem_dims[0] != sem_dims[2])
    {
      std::cout << " dimensions are different. bail out\n";
      return MB_FAILURE;
    }
    // repeat for target sets
    spectral_sets.empty();
    // now initialize a source spectral element !
    _spectralSource = new moab::Element::SpectralHex(sem_dims[0]);
    specSou = true;
  }
  // repeat for target source
  rval = mbImpl->get_entities_by_type_and_tag(rootTarget, moab::MBENTITYSET, &sem_tag, NULL, 1, spectral_sets);
  if (moab::MB_SUCCESS != rval || spectral_sets.empty())
  {
    std::cout<< "can't get sem set on target\n";
  }
  else
  {

    moab::EntityHandle firstSemSet=spectral_sets[0];
    rval = mbImpl->tag_get_data(sem_tag, &firstSemSet, 1, (void*)sem_dims);
    if (moab::MB_SUCCESS != rval)
      return MB_FAILURE;

    if (sem_dims[0]!=sem_dims[1] || sem_dims[0] != sem_dims[2])
    {
      std::cout << " dimensions are different. bail out\n";
      return MB_FAILURE;
    }
    // repeat for target sets
    spectral_sets.empty();
    // now initialize a target spectral element !
    _spectralTarget = new moab::Element::SpectralHex(sem_dims[0]);
    specTar = true;
  }

  _ntot = sem_dims[0]*sem_dims[1]*sem_dims[2];
  rval = mbImpl->tag_get_handle("SEM_X", _ntot, moab::MB_TYPE_DOUBLE, _xm1Tag);
  if (moab::MB_SUCCESS != rval)
  {
     std::cout << "can't get xm1tag \n";
     return MB_FAILURE;
  }
  rval = mbImpl->tag_get_handle("SEM_Y", _ntot, moab::MB_TYPE_DOUBLE, _ym1Tag);
  if (moab::MB_SUCCESS != rval)
  {
     std::cout << "can't get ym1tag \n";
     return MB_FAILURE;
  }
  rval = mbImpl->tag_get_handle("SEM_Z", _ntot, moab::MB_TYPE_DOUBLE, _zm1Tag);
  if (moab::MB_SUCCESS != rval)
  {
     std::cout << "can't get zm1tag \n";
     return MB_FAILURE;
  }
  return MB_SUCCESS;
}

ErrorCode Coupler::locate_points(Range &targ_ents,
                                 double rel_eps, 
                                 double abs_eps,
                                 TupleList *tl,
                                 bool store_local)
{
    // get locations
  std::vector<double> locs(3*targ_ents.size());
  Range verts = targ_ents.subset_by_type(MBVERTEX);
  ErrorCode rval = mbImpl->get_coords(verts, &locs[0]);
  if (MB_SUCCESS != rval) return rval;
    // now get other ents; reuse verts
  unsigned int num_verts = verts.size();
  verts = subtract(targ_ents, verts);
    // compute centroids
  std::vector<EntityHandle> dum_conn(CN::MAX_NODES_PER_ELEMENT);
  std::vector<double> dum_pos(CN::MAX_NODES_PER_ELEMENT);
  const EntityHandle *conn;
  int num_conn;
  double *coords = &locs[num_verts];
    // do this here instead of a function to allow reuse of dum_pos and dum_conn
  for (Range::const_iterator rit = verts.begin(); rit != verts.end(); rit++) {
    rval = mbImpl->get_connectivity(*rit, conn, num_conn, false, &dum_conn);
    if (MB_SUCCESS != rval) return rval;
    rval = mbImpl->get_coords(conn, num_conn, &dum_pos[0]);
    if (MB_SUCCESS != rval) return rval;
    coords[0] = coords[1] = coords[2] = 0.0;
    for (int i = 0; i < num_conn; i++) {
      coords[0] += dum_pos[3*i];
      coords[1] += dum_pos[3*i+1];
      coords[2] += dum_pos[3*i+2];
    }
    coords[0] /= num_conn;
    coords[1] /= num_conn;
    coords[2] /= num_conn;
    coords += 3;
  }

  if (store_local) targetEnts = targ_ents;
  
  return locate_points(&locs[0], targ_ents.size(), rel_eps, abs_eps, tl, store_local);
}

ErrorCode Coupler::locate_points(double *xyz, int num_points,
                                 double rel_eps, 
                                 double abs_eps,
                                 TupleList *tl,
                                 bool store_local)
{
  assert(tl || store_local);

    // target_pts: TL(to_proc, tgt_index, x, y, z): tuples sent to source mesh procs representing pts to be located
    // source_pts: TL(from_proc, tgt_index, src_index): results of source mesh proc point location, ready to send
    //             back to tgt procs; src_index of -1 indicates point not located (arguably not useful...)
  TupleList target_pts;
  target_pts.initialize(2, 0, 0, 3, num_points);
  target_pts.enableWriteAccess();

  TupleList source_pts;
  mappedPts = new TupleList(0, 0, 1, 3, target_pts.get_max());
  mappedPts->enableWriteAccess();

  source_pts.initialize(3, 0, 0, 0, target_pts.get_max()); 
  source_pts.enableWriteAccess();

  mappedPts->set_n(0);
  source_pts.set_n(0);
  ErrorCode result;

    // for each point, find box(es) containing the point,
    // appending results to tuple_list;
    // keep local points separately, in local_pts, which has pairs
    // of <local_index, mapped_index>, where mapped_index is the index
    // into the mappedPts tuple list

  unsigned int my_rank = (myPc ? myPc->proc_config().proc_rank() : 0);
  bool point_located;
  
  for (int i = 0; i < 3*num_points; i+=3) 
  {
    for (unsigned int j = 0; j < (myPc ? myPc->proc_config().proc_size() : 0); j++)
    {
        // test if point is in proc's box
      if ( (allBoxes[6*j] <= xyz[i]+abs_eps) && ( xyz[i] <= allBoxes[6*j+3]+abs_eps) &&
           (allBoxes[6*j+1] <= xyz[i+1]+abs_eps) && (xyz[i+1] <= allBoxes[6*j+4]+abs_eps) &&
           (allBoxes[6*j+2] <= xyz[i+2]+abs_eps) && (xyz[i+2] <= allBoxes[6*j+5]+abs_eps))
      {
          // if in this proc's box, will send to proc to test further
          // check size, grow if we're at max
        if (target_pts.get_n() == target_pts.get_max()) 
          target_pts.resize(std::max(10.0, 1.5*target_pts.get_max()));
  
        target_pts.vi_wr[2*target_pts.get_n()] = j;
        target_pts.vi_wr[2*target_pts.get_n()+1] = i/3;

        target_pts.vr_wr[3*target_pts.get_n()] = xyz[i];
        target_pts.vr_wr[3*target_pts.get_n()+1] = xyz[i+1];
        target_pts.vr_wr[3*target_pts.get_n()+2] = xyz[i+2];
        target_pts.inc_n();
      }
    }
  }

  int num_to_me = 0;
  for (unsigned int i = 0; i < target_pts.get_n(); i++)
    if (target_pts.vi_rd[2*i] == (int)my_rank) num_to_me++;
#ifndef NDEBUG
  printf("rank: %d local points: %d, nb sent target pts: %d mappedPts: %d num to me: %d \n",
         my_rank, num_points, target_pts.get_n(), mappedPts->get_n(), num_to_me);
#endif
  // perform scatter/gather, to gather points to source mesh procs
  if (myPc) {
    (myPc->proc_config().crystal_router())->gs_transfer(1, target_pts, 0);

    num_to_me = 0;
    for (unsigned int i = 0; i < target_pts.get_n(); i++)
      if (target_pts.vi_rd[2*i] == (int)my_rank) num_to_me++;
#ifndef NDEBUG
    printf("rank: %d after first gs nb received_pts: %d; num_from_me = %d\n",
           my_rank, target_pts.get_n(), num_to_me);
#endif
      // after scatter/gather:
      // target_pts.set_n( # points local proc has to map );
      // target_pts.vi_wr[2*i] = proc sending point i
      // target_pts.vi_wr[2*i+1] = index of point i on sending proc
      // target_pts.vr_wr[3*i..3*i+2] = xyz of point i
      //
      // Mapping builds the tuple list:
      // source_pts.set_n (target_pts.get_n() )
      // source_pts.vi_wr[3*i] = target_pts.vi_wr[2*i] = sending proc
      // source_pts.vi_wr[3*i+1] = index of point i on sending proc
      // source_pts.vi_wr[3*i+2] = index of mapped point (-1 if not mapped)
      //
      // Also, mapping builds local tuple_list mappedPts:
      // mappedPts->set_n( # mapped points );
      // mappedPts->vul_wr[i] = local handle of mapped entity
      // mappedPts->vr_wr[3*i..3*i+2] = natural coordinates in mapped entity

      // test target points against my elements
    for (unsigned i = 0; i < target_pts.get_n(); i++) 
    {
      result = test_local_box(target_pts.vr_wr+3*i, 
                              target_pts.vi_rd[2*i], target_pts.vi_rd[2*i+1], i, 
                              point_located, rel_eps, abs_eps, &source_pts);
      if (MB_SUCCESS != result) return result;
    }

      // no longer need target_pts
    target_pts.reset();
#ifndef NDEBUG
    printf("rank: %d nb sent source pts: %d, mappedPts now: %d\n",
           my_rank, source_pts.get_n(),  mappedPts->get_n());
#endif
      // send target points back to target procs
    (myPc->proc_config().crystal_router())->gs_transfer(1, source_pts, 0);

#ifndef NDEBUG
    printf("rank: %d nb received source pts: %d\n",
           my_rank, source_pts.get_n());
#endif
  }

  
    // store proc/index tuples in targetPts, and/or pass back to application;
    // the tuple this gets stored to looks like:
    // tl.set_n( # mapped points );
    // tl.vi_wr[3*i] = remote proc mapping point
    // tl.vi_wr[3*i+1] = local index of mapped point
    // tl.vi_wr[3*i+2] = remote index of mapped point
    //
    // Local index is mapped into either myRange, holding the handles of
    // local mapped entities, or myXyz, holding locations of mapped pts

    // store information about located points
  TupleList *tl_tmp;
  if (!store_local) 
    tl_tmp = tl;
  else {
    targetPts = new TupleList();
    tl_tmp = targetPts;
  }

  tl_tmp->initialize(3, 0, 0, 0, num_points);
  tl_tmp->set_n(num_points); // automatically sets tl to write_enabled
    // initialize so we know afterwards how many pts weren't located
  std::fill(tl_tmp->vi_wr, tl_tmp->vi_wr+3*num_points, -1);

  unsigned int local_pts = 0;
  for (unsigned int i = 0; i < source_pts.get_n(); i++) {
    if (-1 != source_pts.vi_rd[3*i+2]) { //why bother sending message saying "i don't have the point" if it gets discarded?
      int tgt_index = 3*source_pts.vi_rd[3*i+1];
      // prefer always entities that are local, from the source_pts
      // if a local entity was already found to contain the target point, skip
      // tl_tmp->vi_wr[tgt_index] is -1 initially, but it could already be set with
      // a remote processor
      if (tl_tmp->vi_wr[tgt_index] != (int)my_rank)
      {
        tl_tmp->vi_wr[tgt_index]   = source_pts.vi_rd[3*i];
        tl_tmp->vi_wr[tgt_index+1] = source_pts.vi_rd[3*i+1];
        tl_tmp->vi_wr[tgt_index+2] = source_pts.vi_rd[3*i+2];
      }
    }
  }

    // count missing points
  unsigned int missing_pts = 0;
  for (int i = 0; i < num_points; i++) 
  {
    if (tl_tmp->vi_rd[3*i+1] == -1)
      missing_pts++;
    else
      if (tl_tmp->vi_rd[3*i]==(int)my_rank)
        local_pts++;
  }
#ifndef NDEBUG
  printf("rank: %d point location: wanted %d got %u locally, %d remote, missing %d\n",
         my_rank, num_points, local_pts,  num_points-missing_pts-local_pts, missing_pts);
#endif
  assert(0==missing_pts); //will likely break on curved geometries
  
    // no longer need source_pts
  source_pts.reset();

    // copy into tl if passed in and storing locally
  if (tl && store_local) {
    tl = new TupleList(3, 0, 0, 0, num_points);
    tl->enableWriteAccess();
    memcpy(tl->vi_wr, tl_tmp->vi_rd, 3*tl_tmp->get_n()*sizeof(int));
    tl->set_n( tl_tmp->get_n() );
    tl->disableWriteAccess();
  }

  tl_tmp->disableWriteAccess();

    // done
  return MB_SUCCESS;
}

ErrorCode Coupler::test_local_box(double *xyz, 
                                      int from_proc, int remote_index, int /*index*/,
                                      bool &point_located,
                                  double rel_eps, double abs_eps,
                                  TupleList *tl)
{
  std::vector<EntityHandle> entities;
  std::vector<CartVect> nat_coords;
  bool canWrite;
  if (tl) {
    canWrite = tl->get_writeEnabled();
    if(!canWrite) tl->enableWriteAccess();
  }

  if (rel_eps && !abs_eps) {
      // relative epsilon given, translate to absolute epsilon using box dimensions
    BoundBox box;
    myTree->get_bounding_box(box, &localRoot);
    abs_eps = rel_eps * box.diagonal_length();
  }
  
  ErrorCode result = nat_param(xyz, entities, nat_coords, abs_eps);
  if (MB_SUCCESS != result) return result;

    // if we didn't find any ents and we're looking locally, nothing more to do
  if (entities.empty())
  {
    if (tl->get_n() == tl->get_max())
      tl->resize(std::max(10.0, 1.5*tl->get_max()));

    tl->vi_wr[3 * tl->get_n()] = from_proc;
    tl->vi_wr[3 * tl->get_n() + 1] = remote_index;
    tl->vi_wr[3 * tl->get_n() + 2] = -1;
    tl->inc_n();

    point_located = false;
    return MB_SUCCESS;
  }

    // grow if we know we'll exceed size
  if (mappedPts->get_n()+entities.size() >= mappedPts->get_max())
    mappedPts->resize(std::max(10.0, 1.5*mappedPts->get_max()));

  std::vector<EntityHandle>::iterator eit = entities.begin();
  std::vector<CartVect>::iterator ncit = nat_coords.begin();

  mappedPts->enableWriteAccess();
  for (; eit != entities.end(); eit++, ncit++) {
      // store in tuple mappedPts
    mappedPts->vr_wr[3*mappedPts->get_n()] = (*ncit)[0];
    mappedPts->vr_wr[3*mappedPts->get_n()+1] = (*ncit)[1];
    mappedPts->vr_wr[3*mappedPts->get_n()+2] = (*ncit)[2];
    mappedPts->vul_wr[mappedPts->get_n()] = *eit;
    mappedPts->inc_n();

      // also store local point, mapped point indices
    if (tl->get_n() == tl->get_max()) 
      tl->resize(std::max(10.0, 1.5*tl->get_max()));

      // store in tuple source_pts
    tl->vi_wr[3*tl->get_n()] = from_proc;
    tl->vi_wr[3*tl->get_n()+1] = remote_index;
    tl->vi_wr[3*tl->get_n()+2] = mappedPts->get_n()-1;
    tl->inc_n();
  }

  point_located = true;
  
  if(tl && !canWrite) tl->disableWriteAccess();

  return MB_SUCCESS;
}

ErrorCode Coupler::interpolate(Coupler::Method method,
                                      const std::string &interp_tag,
                                      double *interp_vals,
                                      TupleList *tl,
                                      bool normalize)
{
  Tag tag;
  ErrorCode result ;
  if (_spectralSource)
    result = mbImpl->tag_get_handle(interp_tag.c_str(), _ntot, MB_TYPE_DOUBLE, tag);
  else
    result = mbImpl->tag_get_handle(interp_tag.c_str(), 1, MB_TYPE_DOUBLE, tag);
  if (MB_SUCCESS != result) {
    std::ostringstream str;
    str << "Failed to get handle for interpolation tag \"" << interp_tag << "\"";
    mError->set_last_error(str.str());
    return result;
  }
  return interpolate(method, tag, interp_vals, tl, normalize);
}

ErrorCode Coupler::interpolate(Coupler::Method *methods,
                               Tag *tags,
                               int *points_per_method,
                               int num_methods,
                               double *interp_vals,
                               TupleList *tl,
                               bool /*normalize*/)
{
  
    //if (!((LINEAR_FE == method) || (CONSTANT == method)))
    // return MB_FAILURE;

  TupleList *tl_tmp = (tl ? tl : targetPts);
    // remote pts first
  
  ErrorCode result = MB_SUCCESS;

  unsigned int pts_total = 0;
  for (int i = 0; i < num_methods; i++) pts_total += points_per_method[i];

    // if tl was passed in non-NULL, just have those points, otherwise have targetPts plus
    // locally mapped pts
  if (pts_total != tl_tmp->get_n())
    return MB_FAILURE;

  TupleList tinterp;
  tinterp.initialize(5, 0, 0, 1, tl_tmp->get_n());
  int t = 0;
  tinterp.enableWriteAccess();
  for (int i = 0; i < num_methods; i++) {
    for (int j = 0; j < points_per_method[i]; j++) {
      tinterp.vi_wr[5*t] = tl_tmp->vi_rd[3*t];
      tinterp.vi_wr[5*t+1] = tl_tmp->vi_rd[3*t+1];
      tinterp.vi_wr[5*t+2] = tl_tmp->vi_rd[3*t+2];
      tinterp.vi_wr[5*t+3] = methods[i];
      tinterp.vi_wr[5*t+4] = i;
      tinterp.vr_wr[t] = 0.0;
      tinterp.inc_n();
      t++;
    }
  }

    // scatter/gather interpolation points
  if (myPc) {
    (myPc->proc_config().crystal_router())->gs_transfer(1, tinterp, 0);

      // perform interpolation on local source mesh; put results into
      // tinterp.vr_wr[i]

    mappedPts->enableWriteAccess();
    for (unsigned int i = 0; i < tinterp.get_n(); i++) {
      int mindex = tinterp.vi_rd[5*i+2];
      Method method = (Method)tinterp.vi_rd[5*i+3];
      Tag tag = tags[tinterp.vi_rd[5*i+4]];

      result = MB_FAILURE;
      if(LINEAR_FE == method || QUADRATIC_FE == method){
        result = interp_field(mappedPts->vul_rd[mindex],
                              CartVect(mappedPts->vr_wr+3*mindex), 
                              tag, tinterp.vr_wr[i]);
      }else if (CONSTANT == method){
        result = constant_interp(mappedPts->vul_rd[mindex],
                                 tag, tinterp.vr_wr[i]);
      }

      if (MB_SUCCESS != result) return result;
    }
  
      // scatter/gather interpolation data
    myPc->proc_config().crystal_router()->gs_transfer(1, tinterp, 0);
  }

    // copy the interpolated field as a unit
  for (unsigned int i = 0; i < tinterp.get_n(); i++)
    interp_vals[tinterp.vi_rd[5*i+1]] = tinterp.vr_rd[i];
  
    // done
  return MB_SUCCESS;
}

ErrorCode Coupler::nat_param(double xyz[3], 
                                 std::vector<EntityHandle> &entities, 
                             std::vector<CartVect> &nat_coords,
                             double epsilon) 
{
  AdaptiveKDTreeIter treeiter;
  ErrorCode result = myTree->get_tree_iterator(localRoot, treeiter); 
  if (MB_SUCCESS != result) {
    std::cout << "Problems getting iterator" << std::endl;
    return result;
  }

  EntityHandle closest_leaf;
  if (epsilon) {
    std::vector<double> dists;
    std::vector<EntityHandle> leaves;
    result = myTree->distance_search(xyz, epsilon, leaves, 1.0e-10, 1.0e-6, &dists, NULL, &localRoot);
    if (leaves.empty()) 
      // not found returns success here, with empty list, just like case with no epsilon
      return MB_SUCCESS;
      // get closest leaf
    double min_dist = *dists.begin();
    closest_leaf = *leaves.begin();
    std::vector<EntityHandle>::iterator vit = leaves.begin()+1;
    std::vector<double>::iterator dit = dists.begin()+1;
    for (; vit != leaves.end() && min_dist; vit++, dit++) {
      if (*dit < min_dist) {
        min_dist = *dit;
        closest_leaf = *vit;
      }
    }
  }
  else {
    result = myTree->point_search(xyz, treeiter, 1.0e-10, 1.0e-6, NULL, &localRoot);
    if(MB_ENTITY_NOT_FOUND==result) //point is outside of myTree's bounding box
      return MB_SUCCESS; 
    else if (MB_SUCCESS != result) {
      std::cout << "Problems getting leaf \n";
      return result;
    }
    closest_leaf = treeiter.handle();
  }

    // find natural coordinates of point in element(s) in that leaf
  CartVect tmp_nat_coords; 
  Range range_leaf;
  result = mbImpl->get_entities_by_dimension(closest_leaf, max_dim, range_leaf, false);
  if(result != MB_SUCCESS) std::cout << "Problem getting leaf in a range" << std::endl;

  // loop over the range_leaf
  for(Range::iterator iter = range_leaf.begin(); iter != range_leaf.end(); iter++)
  {
    //test to find out in which entity the point is
    //get the EntityType and create the appropriate Element::Map subtype
    // if spectral, do not need coordinates, just the GL points
    EntityType etype = mbImpl->type_from_handle(*iter);
    if (NULL!= this->_spectralSource && etype==MBHEX)
    {
      EntityHandle eh = *iter;
      const double * xval;
      const double * yval;
      const double * zval;
      ErrorCode rval = mbImpl-> tag_get_by_ptr(_xm1Tag, &eh, 1,(const void **) &xval );
      if (moab::MB_SUCCESS != rval)
      {
        std::cout << "can't get xm1 values \n";
        return MB_FAILURE;
      }
      rval = mbImpl-> tag_get_by_ptr(_ym1Tag, &eh, 1, (const void **)&yval );
      if (moab::MB_SUCCESS != rval)
      {
        std::cout << "can't get ym1 values \n";
        return MB_FAILURE;
      }
      rval = mbImpl-> tag_get_by_ptr(_zm1Tag, &eh, 1, (const void **)&zval );
      if (moab::MB_SUCCESS != rval)
      {
        std::cout << "can't get zm1 values \n";
        return MB_FAILURE;
      }
      Element::SpectralHex * spcHex = ( Element::SpectralHex * ) _spectralSource;

      spcHex->set_gl_points((double*)xval, (double*)yval, (double*)zval);
      tmp_nat_coords = spcHex->ievaluate(CartVect(xyz));
      bool inside = spcHex->inside_nat_space(CartVect(xyz), epsilon);
      if (!inside) {
        std::cout << "point "<< xyz[0] << " " << xyz[1] << " " << xyz[2] <<
            " is not converging inside hex " << mbImpl->id_from_handle(eh) << "\n";
        continue; // it is possible that the point is outside, so it will not converge
      }
    }
    else
    {
      const EntityHandle *connect;
      int num_connect;

        //get connectivity
      result = mbImpl->get_connectivity(*iter, connect, num_connect, true);

        //get coordinates of the vertices
      std::vector<CartVect> coords_vert(num_connect);
      result = mbImpl->get_coords(connect, num_connect, &(coords_vert[0][0]));
      if (MB_SUCCESS != result) {
        std::cout << "Problems getting coordinates of vertices\n";
        return result;
      }

      if (etype == MBHEX) {
        if (8==num_connect)
        {
          Element::LinearHex hexmap(coords_vert);
          tmp_nat_coords = hexmap.ievaluate(CartVect(xyz), epsilon);
          bool inside = hexmap.inside_nat_space(tmp_nat_coords, epsilon);
          if (!inside)
            continue;
        }
        else if (27==num_connect)
        {
          Element::QuadraticHex hexmap(coords_vert);
          tmp_nat_coords = hexmap.ievaluate(CartVect(xyz), epsilon);
          bool inside = hexmap.inside_nat_space(tmp_nat_coords, epsilon);
          if (!inside)
            continue;
        }
        else // TODO this case not treated yet, no interpolation
          continue;
      }
      else if (etype == MBTET){
        Element::LinearTet tetmap(coords_vert);
        tmp_nat_coords = tetmap.ievaluate(CartVect(xyz));
        bool inside = tetmap.inside_nat_space(tmp_nat_coords, epsilon);
        if (!inside)
          continue;
      }
      else if (etype == MBQUAD){
        Element::LinearQuad quadmap(coords_vert);
        try {
          tmp_nat_coords = quadmap.ievaluate(CartVect(xyz), epsilon);
          bool inside = quadmap.inside_nat_space(tmp_nat_coords, epsilon);
          if (!inside) continue;
        }
        catch (Element::Map::EvaluationError) {
          continue;
        }
        if (!quadmap.inside_nat_space(tmp_nat_coords, epsilon))
          continue;
      }
      else {
        std::cout << "Entity not Hex or Tet or Quad" << std::endl;
        continue;
      }
    }
      //if we get here then we've found the coordinates.
      //save them and the entity and return success.
    entities.push_back(*iter);
    nat_coords.push_back(tmp_nat_coords);
    return MB_SUCCESS;
  }

  //didn't find any elements containing the point
  return MB_SUCCESS;
}

ErrorCode Coupler::interp_field(EntityHandle elem,
				CartVect nat_coord, 
				Tag tag,
				double &field)
{
  if (_spectralSource)
  {
    // get tag values at the GL points for some field (Tag)
    const double * vx;
    ErrorCode rval = mbImpl-> tag_get_by_ptr(tag, &elem, 1, (const void **) &vx);
    if (moab::MB_SUCCESS != rval)
    {
      std::cout << "can't get field values for the tag \n";
      return MB_FAILURE;
    }
    Element::SpectralHex * spcHex = (Element::SpectralHex *) _spectralSource;
    field = spcHex->evaluate_scalar_field(nat_coord, vx);
  }
  else
  {
    double vfields[27]; // will work for linear hex, quadratic hex or Tets
    moab::Element::Map *elemMap;
    int num_verts = 0;
    // get the EntityType
    // get the tag values at the vertices
    const EntityHandle *connect;
    int num_connect;
    ErrorCode result = mbImpl->get_connectivity(elem, connect, num_connect);
    if (MB_SUCCESS != result)
    {
      return result;
    }
    EntityType etype = mbImpl->type_from_handle(elem);
    if (etype == MBHEX) {
      if (num_connect==8) {
      elemMap = new moab::Element::LinearHex();
      num_verts = 8;
    }
      else { /* (etype == MBHEX && num_connect==27) */
      elemMap = new moab::Element::QuadraticHex();
      num_verts = 27;
    }
    }
    else if (etype == MBTET) {
      elemMap = new moab::Element::LinearTet();
      num_verts = 4;
    }
    else if (etype == MBQUAD) {
      elemMap = new moab::Element::LinearQuad();
      num_verts = 4;
    }
    else {
      return MB_FAILURE;
    }

    result = mbImpl->tag_get_data(tag, connect, std::min(num_verts, num_connect), vfields);
    if (MB_SUCCESS != result) {
      delete elemMap;
      return result;
    }

    //function for the interpolation
    field = 0;

    // check the number of vertices
    assert(num_connect >= num_verts);

    //calculate the field
    try {
      field = elemMap->evaluate_scalar_field(nat_coord, vfields);
    }
    catch (moab::Element::Map::EvaluationError)
    {
      delete elemMap;
      return MB_FAILURE;
    }

    delete elemMap;
  }
  return MB_SUCCESS;
}


//Simplest "interpolation" for element-based source fields. Set the value of the field
//at the target point to that of the field in the source element it lies in.
ErrorCode Coupler::constant_interp(EntityHandle elem,
                                   Tag tag,
                                   double &field)
{
  double tempField;

  // get the tag values at the vertices
  ErrorCode result = mbImpl->tag_get_data(tag, &elem, 1, &tempField);
  if (MB_SUCCESS != result) return result;

  field = tempField;

  return MB_SUCCESS;
}

// Normalize a field over the entire mesh represented by the root_set.
int Coupler::normalize_mesh(iBase_EntitySetHandle &root_set,
                            const char            *norm_tag,
                            Coupler::IntegType    integ_type,
                            int                   num_integ_pts)
{
  int err = iBase_SUCCESS;

  // Get an iMesh_Instance from MBCoupler::mbImpl.
  iMesh_Instance iMeshInst = reinterpret_cast<iMesh_Instance>(mbImpl);

  // SLAVE START ****************************************************************
  // Search for entities based on tag_handles and tag_values

  std::vector< std::vector<iBase_EntitySetHandle> > entity_sets;
  std::vector< std::vector<iBase_EntityHandle> >    entity_groups;

  // put the root_set into entity_sets
  std::vector<iBase_EntitySetHandle> ent_set;
  ent_set.push_back(root_set);
  entity_sets.push_back(ent_set);

  // get all entities from root_set and put into entity_groups
  std::vector<iBase_EntityHandle> entities;
  iBase_EntityHandle *ents = NULL;
  int ents_alloc = 0;
  int ents_size = 0;

  iMesh_getEntities(iMeshInst, root_set, iBase_ALL_TYPES, iMesh_ALL_TOPOLOGIES,
                    &ents, &ents_alloc, &ents_size, &err);
  ERRORR("iMesh_getEntities failed on root_set.", err);

  // put all of the entities from the entity set into ent_set and free the memory for ents.
  for (int k = 0; k < ents_size; k++) {
    entities.push_back(ents[k]);
  }
  free(ents);

  entity_groups.push_back(entities);

  // Call do_normalization() to continue common normalization processing
  err = do_normalization(norm_tag, entity_sets, entity_groups, integ_type, num_integ_pts);
  ERRORR("Failure in do_normalization().", err);
  // SLAVE END   ****************************************************************

  return err;
}

// Normalize a field over the subset of entities identified by the tags and values passed
int Coupler::normalize_subset(iBase_EntitySetHandle &root_set,
                              const char            *norm_tag,
                              const char            **tag_names,
                              int                   num_tags,
                              const char            **tag_values,
                              Coupler::IntegType    integ_type,
                              int                   num_integ_pts)
{
  iMesh_Instance iMeshInst = reinterpret_cast<iMesh_Instance>(mbImpl);
  int err;
  std::vector<iBase_TagHandle> tag_handles;
  
  // Lookup tag handles from tag names
  for (int t = 0; t < num_tags; t++) {
    // get tag handle & size
    iBase_TagHandle th;
    iMesh_getTagHandle(iMeshInst, tag_names[t], &th, &err, strlen(tag_names[t]));
    ERRORR("Failed to get tag handle.", err);
    tag_handles.push_back(th);
  }

  return normalize_subset(root_set, 
                          norm_tag, 
                          &tag_handles[0], 
                          num_tags, 
                          tag_values, 
                          integ_type, 
                          num_integ_pts);
}

int Coupler::normalize_subset(iBase_EntitySetHandle &root_set,
                              const char            *norm_tag,
                              iBase_TagHandle       *tag_handles,
                              int                   num_tags,
                              const char            **tag_values,
                              Coupler::IntegType    integ_type,
                              int                   num_integ_pts)
{
  int err = iBase_SUCCESS;

  // SLAVE START ****************************************************************
  // Search for entities based on tag_handles and tag_values

  std::vector< std::vector<iBase_EntitySetHandle> > entity_sets;
  std::vector< std::vector<iBase_EntityHandle> >    entity_groups;
  err = get_matching_entities(root_set, tag_handles, tag_values, num_tags, 
                              &entity_sets, &entity_groups);
  ERRORR("Failed to get matching entities.", err);

  // Call do_normalization() to continue common normalization processing
  err = do_normalization(norm_tag, entity_sets, entity_groups, integ_type, num_integ_pts);
  ERRORR("Failure in do_normalization().", err);
  // SLAVE END   ****************************************************************

  return err;
}

int Coupler::do_normalization(const char                                        *norm_tag,
                              std::vector< std::vector<iBase_EntitySetHandle> > &entity_sets,
                              std::vector< std::vector<iBase_EntityHandle> >    &entity_groups,
                              Coupler::IntegType                                integ_type,
                              int                                               num_integ_pts)
{
  // SLAVE START ****************************************************************
  int err = iBase_SUCCESS;

  // Setup data for parallel computing
  int nprocs, rank;
  err = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  ERRORMPI("Getting number of procs failed.", err);
  err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ERRORMPI("Getting rank failed.", err);

  // Get the integrated field value for each group(vector) of entities.
  // If no entities are in a group then a zero will be put in the list
  // of return values.
  unsigned int num_ent_grps = entity_groups.size();
  std::vector<double> integ_vals(num_ent_grps);

  err = get_group_integ_vals(entity_groups, integ_vals, norm_tag, num_integ_pts, integ_type);
  ERRORR("Failed to get integrated field values for groups in mesh.", err);
  // SLAVE END   ****************************************************************

  // SLAVE/MASTER START #########################################################
  // Send list of integrated values back to master proc.  The ordering of the 
  // values will match the ordering of the entity groups (i.e. vector of vectors)
  // sent from master to slaves earlier.  The values for each entity group will
  // be summed during the transfer.
  std::vector<double> sum_integ_vals(num_ent_grps);

  if (nprocs > 1) {
    // If parallel then send the values back to the master.
    err = MPI_Reduce(&integ_vals[0], &sum_integ_vals[0], num_ent_grps, MPI_DOUBLE, 
                     MPI_SUM, MASTER_PROC, myPc->proc_config().proc_comm());
    ERRORMPI("Transfer and reduction of integrated values failed.", err);
  }
  else {
    // Otherwise just copy the vector
    sum_integ_vals = integ_vals;
  }
  // SLAVE/MASTER END   #########################################################

  // MASTER START ***************************************************************
  // Calculate the normalization factor for each group by taking the
  // inverse of each integrated field value.  Put the normalization factor 
  // for each group back into the list in the same order.

  for (unsigned int i = 0; i < num_ent_grps; i++) {
    double val = sum_integ_vals[i];
    if (val != 0) sum_integ_vals[i] = 1.0/val;
  }
  // MASTER END   ***************************************************************

  // MASTER/SLAVE START #########################################################
  if (nprocs > 1) {
    // If parallel then broadcast the normalization factors to the procs.
    err = MPI_Bcast(&sum_integ_vals[0], num_ent_grps, MPI_DOUBLE, MASTER_PROC, 
                    myPc->proc_config().proc_comm());
    ERRORMPI("Broadcast of normalization factors failed.", err);
  }
  // MASTER/SLAVE END   #########################################################

  // SLAVE START ****************************************************************
  // Save the normalization factors to a new tag with name of norm_tag's value
  // and the string "_normF" appended.  This new tag will be created on the entity
  // set that contains all of the entities from a group.

  err = apply_group_norm_factor(entity_sets, sum_integ_vals, norm_tag, integ_type);
  ERRORR("Failed to set the normalization factor for groups in mesh.", err);
  // SLAVE END   ****************************************************************

  return err;
}

// Functions supporting the subset normalization function

// Retrieve groups of entities matching tags and values if present
int Coupler::get_matching_entities(iBase_EntitySetHandle                             root_set,
                                   const char                                        **tag_names,
                                   const char                                        **tag_values,
                                   int                                               num_tags,
                                   std::vector< std::vector<iBase_EntitySetHandle> > *entity_sets,
                                   std::vector< std::vector<iBase_EntityHandle> >    *entity_groups)
{
  iMesh_Instance iMeshInst = reinterpret_cast<iMesh_Instance>(mbImpl);
  int err;
  std::vector<iBase_TagHandle> tag_handles;
  
  for (int t = 0; t < num_tags; t++) {
    // get tag handle & size
    iBase_TagHandle th;
    iMesh_getTagHandle(iMeshInst, tag_names[t], &th, &err, strlen(tag_names[t]));
    ERRORR("Failed to get tag handle.", err);
    tag_handles.push_back(th);
  }

  return get_matching_entities(root_set, &tag_handles[0], tag_values, num_tags,
                               entity_sets, entity_groups);
}

// Retrieve groups of entities matching tags and values if present
int Coupler::get_matching_entities(iBase_EntitySetHandle                          root_set,
                                   iBase_TagHandle                                *tag_handles,
                                   const char                                     **tag_values,
                                   int                                            num_tags,
                                   std::vector< std::vector<iBase_EntitySetHandle> > *entity_sets,
                                   std::vector< std::vector<iBase_EntityHandle> > *entity_groups)
{                                        

  // SLAVE START ****************************************************************
  // Get an iMesh_Instance from MBCoupler::mbImpl.
  iMesh_Instance iMeshInst = reinterpret_cast<iMesh_Instance>(mbImpl);

  // Setup data for parallel computing
  int err = iBase_SUCCESS;
  int nprocs, rank;
  err = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  ERRORMPI("Getting number of procs failed.", err);
  err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ERRORMPI("Getting rank failed.", err);

  int ent_sets_size;
  int ent_sets_alloc=0;
  iBase_EntitySetHandle *ent_sets = NULL;  // free at end

  // Get Entity Sets that match the tags and values.
  iMesh_getEntSetsByTagsRec(iMeshInst, root_set, tag_handles, 
                            tag_values, num_tags, 0,
                            &ent_sets, &ent_sets_alloc, &ent_sets_size, &err);
  ERRORR("iMesh_getEntSetsByTagsRec failed.", err);

  TupleList *tag_list = NULL;
  err = create_tuples(ent_sets, ent_sets_size, tag_handles, num_tags, &tag_list);
  ERRORR("Failed to create tuples from entity sets.", err);

  // Free up array memory from iMesh call
  free(ent_sets);
  ent_sets = NULL;
  ent_sets_alloc = 0;
  ent_sets_size = 0;
  // SLAVE END   ****************************************************************

  // If we are running in a mult-proc session then send tuple list back to master 
  // proc for consolidation.  Otherwise just copy the pointer to the tuple_list.
  TupleList *cons_tuples;
  if (nprocs > 1) {
    // SLAVE/MASTER START #########################################################

    // pack the tuple_list in a buffer.
    uint *tuple_buf;
    int tuple_buf_len;
    tuple_buf_len = pack_tuples(tag_list, (void**)&tuple_buf);

    // Free tag_list here as its not used again if nprocs > 1
    tag_list->reset();

    // Send back the buffer sizes to the master proc
    int *recv_cnts = (int*) malloc(nprocs * sizeof(int));
    int *offsets   = (int*) malloc(nprocs * sizeof(int));
    uint *all_tuples_buf = 0;

    MPI_Gather(&tuple_buf_len, 1, MPI_INT, recv_cnts, 1, MPI_INT, MASTER_PROC, 
               myPc->proc_config().proc_comm());
    ERRORMPI("Gathering buffer sizes failed.", err);

    // Allocate a buffer large enough for all the data
    if (rank == MASTER_PROC) {
      int all_tuples_len = recv_cnts[0];
      offsets[0] = 0;
      for (int i = 1; i < nprocs; i++) {
        offsets[i] = offsets[i-1] + recv_cnts[i-1];
        all_tuples_len += recv_cnts[i];
      }

      all_tuples_buf = (uint*) malloc(all_tuples_len * sizeof(uint));
    }

    // Send all buffers to the master proc for consolidation
    MPI_Gatherv(tuple_buf, tuple_buf_len, MPI_INT, 
                all_tuples_buf, recv_cnts, offsets, MPI_INT, MASTER_PROC,
                myPc->proc_config().proc_comm());
    ERRORMPI("Gathering tuple_lists failed.", err);
    free(tuple_buf);  // malloc'd in pack_tuples

    if (rank == MASTER_PROC) {
      // unpack the tuple_list from the buffer.
      TupleList **tl_array = (TupleList **) malloc(nprocs * sizeof(TupleList*));
      for (int i = 0; i < nprocs; i++)
        unpack_tuples((void*) &all_tuples_buf[offsets[i]], &tl_array[i]);

      // Free all_tuples_buf here as it is only allocated on the MASTER_PROC
      free(all_tuples_buf);
      // SLAVE/MASTER END   #########################################################
      
      // MASTER START ***************************************************************
      // Consolidate all tuple_lists into one tuple_list with no duplicates.
      err = consolidate_tuples(tl_array, nprocs, &cons_tuples);
      ERRORR("Failed to consolidate tuples.", err);

      for (int i = 0; i < nprocs; i++)
        tl_array[i]->reset();
      free(tl_array);
      // MASTER END   ***************************************************************
    }

    // Free offsets and recv_cnts as they are allocated on all procs
    free(offsets);
    free(recv_cnts);

    // MASTER/SLAVE START #########################################################
    // Broadcast condensed tuple list back to all procs.
    uint *ctl_buf;
    int ctl_buf_sz;
    if (rank == MASTER_PROC) {
      ctl_buf_sz = pack_tuples(cons_tuples, (void**) &ctl_buf);
    }

    // Send buffer size
    err = MPI_Bcast(&ctl_buf_sz, 1, MPI_INT, MASTER_PROC, myPc->proc_config().proc_comm());
    ERRORMPI("Broadcasting tuple_list size failed.", err);

    // Allocate a buffer in the other procs
    if (rank != MASTER_PROC) {
      ctl_buf = (uint*) malloc(ctl_buf_sz*sizeof(uint));
    }

    err = MPI_Bcast(ctl_buf, ctl_buf_sz, MPI_INT, MASTER_PROC, myPc->proc_config().proc_comm());
    ERRORMPI("Broadcasting tuple_list failed.", err);

    if (rank != MASTER_PROC) {
      unpack_tuples(ctl_buf, &cons_tuples);
    }
    free(ctl_buf);
    // MASTER/SLAVE END   #########################################################
  }
  else {
    cons_tuples = tag_list;
  }

  // SLAVE START ****************************************************************
  // Loop over the tuple list getting the entities with the tags in the tuple_list entry
  uint mi, ml, mul, mr;
  cons_tuples->getTupleSize(mi, ml, mul, mr);

  for (unsigned int i = 0; i < cons_tuples->get_n(); i++) {
    // Get Entity Sets that match the tags and values.

    // Convert the data in the tuple_list to an array of pointers to the data
    // in the tuple_list as that is what the iMesh API call is expecting.
    int **vals = (int**) malloc(mi * sizeof(int*));
    for (unsigned int j = 0; j < mi; j++)
      vals[j] = (int *)&(cons_tuples->vi_rd[(i*mi) + j]);

    iMesh_getEntSetsByTagsRec(iMeshInst, root_set, tag_handles, 
                              (const char * const *) vals,
                              mi, 0,
                              &ent_sets, &ent_sets_alloc, &ent_sets_size, &err);

    ERRORR("iMesh_getEntSetsByTagsRec failed.", err);
    if (debug) std::cout << "ent_sets_size=" << ent_sets_size << std::endl;

    // Free up the array of pointers
    free(vals);

    // Loop over the entity sets and then free the memory for ent_sets.
    std::vector<iBase_EntitySetHandle> ent_set_hdls;
    std::vector<iBase_EntityHandle> ent_hdls;
    for (int j = 0; j < ent_sets_size; j++) {
      // Save the entity set
      ent_set_hdls.push_back(ent_sets[j]);

      // Get all entities for the entity set
      iBase_EntityHandle *ents = NULL;
      int ents_alloc = 0;
      int ents_size = 0;

      iMesh_getEntities(iMeshInst, ent_sets[j], iBase_ALL_TYPES, iMesh_ALL_TOPOLOGIES,
                        &ents, &ents_alloc, &ents_size, &err);
      ERRORR("iMesh_getEntities failed.", err);
      if (debug) std::cout << "ents_size=" << ents_size << std::endl;

      // Save all of the entities from the entity set and free the memory for ents.
      for (int k = 0; k < ents_size; k++) {
        ent_hdls.push_back(ents[k]);
      }
      free(ents);
      if (debug) std::cout << "ent_hdls.size=" << ent_hdls.size() << std::endl;
    }

    // Free the entity set list for next tuple iteration.
    free(ent_sets);
    ent_sets = NULL;
    ent_sets_alloc = 0;
    ent_sets_size = 0;

    // Push ent_set_hdls onto entity_sets, ent_hdls onto entity_groups
    // and clear both ent_set_hdls and ent_hdls.
    entity_sets->push_back(ent_set_hdls);
    ent_set_hdls.clear();
    entity_groups->push_back(ent_hdls);
    ent_hdls.clear();
    if (debug) std::cout << "entity_sets->size=" << entity_sets->size() 
                         << ", entity_groups->size=" << entity_groups->size() << std::endl;
  }

  cons_tuples->reset();
  // SLAVE END   ****************************************************************

  return err;
}


// Return a tuple_list containing  tag values for each Entity Set
// The tuple_list will have a column for each tag and a row for each
// Entity Set.  It is assumed all of the tags are integer tags.
int Coupler::create_tuples(iBase_EntitySetHandle *ent_sets,
                           int                   num_sets, 
                           const char            **tag_names,
                           int                   num_tags,
                           TupleList            **tuple_list)
{
  iMesh_Instance iMeshInst = reinterpret_cast<iMesh_Instance>(mbImpl);
  int err;
  std::vector<iBase_TagHandle> tag_handles;
  
  for (int t = 0; t < num_tags; t++) {
    // get tag handle & size
    iBase_TagHandle th;
    iMesh_getTagHandle(iMeshInst, tag_names[t], &th, &err, strlen(tag_names[t]));
    ERRORR("Failed to get tag handle.", err);
    tag_handles.push_back(th);
  }

  return create_tuples(ent_sets, num_sets, &tag_handles[0], num_tags, tuple_list);
}

// Return a tuple_list containing  tag values for each Entity Set
// The tuple_list will have a column for each tag and a row for each
// Entity Set.  It is assumed all of the tags are integer tags.
int Coupler::create_tuples(iBase_EntitySetHandle *ent_sets, 
                           int                   num_sets, 
                           iBase_TagHandle       *tag_handles,
                           int                   num_tags,
                           TupleList            **tuples)
{
  // Get an iMesh_Instance from MBCoupler::mbImpl.
  iMesh_Instance iMeshInst = reinterpret_cast<iMesh_Instance>(mbImpl);

  int err = iBase_SUCCESS;

  // ASSUMPTION: All tags are of type integer.  This may need to be expanded in future.

  // Allocate a tuple_list for the number of entity sets passed in
  TupleList *tag_tuples = new TupleList(num_tags, 0, 0, 0, num_sets);
  //tag_tuples->initialize(num_tags, 0, 0, 0, num_sets);
  uint mi, ml, mul, mr;
  tag_tuples->getTupleSize(mi, ml, mul, mr);
  tag_tuples->enableWriteAccess();

  if (mi == 0)
    ERRORR("Failed to initialize tuple_list.", iBase_FAILURE);

  // Loop over the filtered entity sets retrieving each matching tag value one by one.
  int val;
  for (int i = 0; i < num_sets; i++) {
    for (int j = 0; j < num_tags; j++) {
      iMesh_getEntSetIntData(iMeshInst, ent_sets[i], tag_handles[j], &val, &err);
      ERRORR("Failed to get integer tag data.", err);
      tag_tuples->vi_wr[i*mi + j] = val;
    }

    // If we get here there was no error so increment n in the tuple_list
    tag_tuples->inc_n();
  }
  tag_tuples->disableWriteAccess();
  *tuples = tag_tuples;

  return err;
}

// Consolidate tuple_lists into one list with no duplicates
int Coupler::consolidate_tuples(TupleList **all_tuples, 
                                int        num_tuples,
                                TupleList **unique_tuples)
{
  int err = iBase_SUCCESS;

  int total_rcv_tuples = 0;
  int offset = 0, copysz = 0;
  unsigned num_tags = 0;

  uint ml, mul, mr;
  uint *mi = (uint *)malloc(sizeof(uint) * num_tuples);

  for(int i=0; i<num_tuples; i++){
    all_tuples[i]->getTupleSize(mi[i], ml, mul, mr);
  }

  for (int i = 0; i < num_tuples; i++) {
    if (all_tuples[i] != NULL) {
      total_rcv_tuples += all_tuples[i]->get_n();
      num_tags = mi[i];
    }
  }
  const unsigned int_size = sizeof(sint);
  const unsigned int_width = num_tags * int_size;

  // Get the total size of all of the tuple_lists in all_tuples.
  for (int i = 0; i < num_tuples; i++) {
    if (all_tuples[i] != NULL)
      total_rcv_tuples += all_tuples[i]->get_n();
  }

  // Copy the tuple_lists into a single tuple_list.
  TupleList *all_tuples_list = new TupleList(num_tags, 0, 0, 0, total_rcv_tuples);
  all_tuples_list->enableWriteAccess();
  //all_tuples_list->initialize(num_tags, 0, 0, 0, total_rcv_tuples);
  for (int i = 0; i < num_tuples; i++) {
    if (all_tuples[i] != NULL) {
      copysz = all_tuples[i]->get_n() * int_width;
      memcpy(all_tuples_list->vi_wr+offset, all_tuples[i]->vi_rd, copysz);
      offset = offset + (all_tuples[i]->get_n() * mi[i]);
      all_tuples_list->set_n( all_tuples_list->get_n() + all_tuples[i]->get_n() );
    }
  }

  // Sort the new tuple_list.  Use a radix type sort, starting with the last (or least significant)
  // tag column in the vi array and working towards the first (or most significant) tag column.
  TupleList::buffer sort_buffer;
  sort_buffer.buffer_init(2 * total_rcv_tuples * int_width);
  for (int i = num_tags - 1; i >= 0; i--) {
    all_tuples_list->sort(i, &sort_buffer);
  }

  // Cycle through the sorted list eliminating duplicates.
  // Keep counters to the current end of the tuple_list (w/out dups) and the last tuple examined.
  unsigned int end_idx = 0, last_idx = 1;
  while (last_idx < all_tuples_list->get_n()) {
    if (memcmp(all_tuples_list->vi_rd+(end_idx*num_tags), all_tuples_list->vi_rd+(last_idx*num_tags), int_width) == 0) {
      // Values equal - skip
      last_idx += 1;
    }
    else {
      // Values different - copy
      // Move up the end index
      end_idx += 1;
      memcpy(all_tuples_list->vi_wr+(end_idx*num_tags), all_tuples_list->vi_rd+(last_idx*num_tags), int_width);
      last_idx += 1;
    }
  }
  // Update the count in all_tuples_list
  all_tuples_list->set_n( end_idx + 1 );

  // Resize the tuple_list
  all_tuples_list->resize(all_tuples_list->get_n());

  // Set the output parameter
  *unique_tuples = all_tuples_list;

  return err;
}

// Calculate integrated field values for groups of entities
int Coupler::get_group_integ_vals(std::vector< std::vector<iBase_EntityHandle> > &groups,
                                  std::vector<double> &integ_vals,
                                  const char *norm_tag,
                                  int /*num_integ_vals*/,
                                  Coupler::IntegType integ_type)
{
  // Get an iMesh_Instance from MBCoupler::mbImpl.
  iMesh_Instance iMeshInst = reinterpret_cast<iMesh_Instance>(mbImpl);

  int err = iBase_SUCCESS;

  std::vector< std::vector<iBase_EntityHandle> >::iterator iter_i;
  std::vector<iBase_EntityHandle>::iterator iter_j;
  double grp_intrgr_val, intgr_val;

  // Get the tag handle for norm_tag
  iBase_TagHandle norm_hdl;
  iMesh_getTagHandle(iMeshInst, norm_tag, &norm_hdl, &err, strlen(norm_tag));
  ERRORR("Failed to get norm_tag handle.", err);

  // Check size of integ_vals vector
  if (integ_vals.size() != groups.size())
    integ_vals.resize(groups.size());

  // Loop over the groups(vectors) of entities
  unsigned int i;
  for (i = 0, iter_i = groups.begin(); iter_i != groups.end(); i++, iter_i++) {
    grp_intrgr_val = 0;

    // Loop over the all the entities in the group, integrating 
    // the field_fn over the entity in iter_j
    for (iter_j = (*iter_i).begin(); iter_j != (*iter_i).end(); iter_j++) {
      // Check that the entity in iter_j is of the same dimension as the 
      // integ_type we are performing
      int j_type;
      iMesh_getEntType(iMeshInst, (*iter_j), &j_type, &err);
      ERRORR("Failed to get entity type.", err);
      // Skip any entities in the group that are not of the type being considered
      if ((integ_type == VOLUME) && (j_type != iBase_REGION))
        continue;
      
      intgr_val = 0;

      // Check to see if this is a HEXAHEDRON or a TETRAHEDRON
      // and setup parameters accordingly.
      int topo_type;
      iMesh_getEntTopo(iMeshInst, (*iter_j), &topo_type, &err);
      ERRORR("Failed to get topology for entity.", err);

      moab::Element::Map *elemMap;
      int num_verts = 0;
      if (topo_type == iMesh_HEXAHEDRON) {
        elemMap = new moab::Element::LinearHex();
        num_verts = 8;
      }
      else if (topo_type == iMesh_TETRAHEDRON) {
        elemMap = new moab::Element::LinearTet();
        num_verts = 4;
      }
      else
        ERRORR("Unknown topology type.", iBase_NOT_SUPPORTED);

      // Get coordinates of all corner vertices (in normal order) and
      // put in array of CartVec.
      std::vector<CartVect> vertices(num_verts);

      // Retrieve the vertices from the element
      iBase_EntityHandle *verts = NULL;
      int verts_alloc = 0;
      int verts_size = 0;

      iMesh_getEntAdj(iMeshInst, (*iter_j), iBase_VERTEX, &verts, &verts_alloc, &verts_size, &err);
      if (iBase_SUCCESS != err) {
        std::cerr << "Failed to get vertices from entity." << std::endl;
        delete(elemMap);
        return err;
      }

      if (verts_size != num_verts) {
        std::cerr << "Failed to get correct number of vertices." << std::endl;
        delete(elemMap);
        free(verts);
        return iBase_FAILURE;
      }

      // Get the vertex coordinates and the field values at the vertices.
      double *coords = NULL;
      int coords_alloc = 0;
      int coords_size = 0;
      iMesh_getVtxArrCoords(iMeshInst, verts, verts_size, iBase_INTERLEAVED, &coords, &coords_alloc, &coords_size, &err);
      if (iBase_SUCCESS != err) {
        std::cerr << "Failed to get vertex coordinates." << std::endl;
        delete(elemMap);
        free(verts);
        return err;
      }

      double *vfield = NULL;
      int vfield_alloc = 0;
      int vfield_size = 0;
      iMesh_getDblArrData(iMeshInst, verts, verts_size, norm_hdl, &vfield, &vfield_alloc, &vfield_size, &err);
      if (iBase_SUCCESS != err) {
        std::cerr << "Failed to get vertex double data for norm_tag." << std::endl;
        delete(elemMap);
        free(verts);
        free(coords);
        return err;
      }

      // Put the vertices into a CartVect vector
      double *x = coords;
      for (int j = 0; j < verts_size; j++, x+=3) {
        vertices[i] = CartVect(x);
      }
      free(verts);
      free(coords);

      // Set the vertices in the Map and perform the integration
      try {
        elemMap->set_vertices(vertices);
        intgr_val = elemMap->integrate_scalar_field(vfield);

        // Combine the result with those of the group
        grp_intrgr_val += intgr_val;
      }
      catch (moab::Element::Map::ArgError) {
        std::cerr << "Failed to set vertices on Element::Map." << std::endl;
        delete(elemMap);
        free(vfield);
      }
      catch (moab::Element::Map::EvaluationError) {
        std::cerr << "Failed to get inverse evaluation of coordinate on Element::Map." << std::endl;
        delete(elemMap);
        free(vfield);
      }

      delete(elemMap);
      free(vfield);
    }

    // Set the group integrated value in the vector
    integ_vals[i] = grp_intrgr_val;
  }

  return err;
}

// Apply a normalization factor to group of entities
int Coupler::apply_group_norm_factor(std::vector< std::vector<iBase_EntitySetHandle> > &entity_sets,
                                     std::vector<double> &norm_factors, 
                                     const char *norm_tag,
                                     Coupler::IntegType /*integ_type*/)
{
  // Get an iMesh_Instance from MBCoupler::mbImpl.
  iMesh_Instance iMeshInst = reinterpret_cast<iMesh_Instance>(mbImpl);

  int err = iBase_SUCCESS;

  // Construct the new tag for the normalization factor from the norm_tag name
  // and "_normf".
  int norm_tag_len = strlen(norm_tag);
  const char* normf_appd = "_normf";
  int normf_appd_len = strlen(normf_appd);

  char* normf_tag = (char*) malloc(norm_tag_len + normf_appd_len + 1);
  char* tmp_ptr = normf_tag;

  memcpy(tmp_ptr, norm_tag, norm_tag_len);
  tmp_ptr += norm_tag_len;
  memcpy(tmp_ptr, normf_appd, normf_appd_len);
  tmp_ptr += normf_appd_len;
  *tmp_ptr = '\0';

  iBase_TagHandle normf_hdl = NULL;
  // Check to see if the tag exists.  If not then create it.
  iMesh_getTagHandle(iMeshInst, normf_tag, &normf_hdl, &err, strlen(normf_tag));
  if (normf_hdl == NULL) {
    iMesh_createTag(iMeshInst, normf_tag, 1, iBase_DOUBLE, &normf_hdl, &err, strlen(normf_tag));
    ERRORR("Failed to create normalization factor tag.", err);
  }
  free(normf_tag);

  std::vector< std::vector<iBase_EntitySetHandle> >::iterator iter_i;
  std::vector<iBase_EntitySetHandle>::iterator iter_j;
  std::vector<double>::iterator iter_f;
  double grp_norm_factor = 0.0;

  // Loop over the entity sets
  for (iter_i = entity_sets.begin(), iter_f = norm_factors.begin(); 
       (iter_i != entity_sets.end()) && (iter_f != norm_factors.end()); 
       iter_i++, iter_f++) {
    grp_norm_factor = *iter_f;

    // Loop over the all the entity sets in iter_i and set the 
    // new normf_tag with the norm factor value on each
    for (iter_j = (*iter_i).begin(); iter_j != (*iter_i).end(); iter_j++) {

      iMesh_setEntSetDblData(iMeshInst, (*iter_j), normf_hdl, grp_norm_factor, &err);
      ERRORR("Failed to set normalization factor on entity set.", err);
    }
  }

  return err;
}

#define UINT_PER_X(X) ((sizeof(X)+sizeof(uint)-1)/sizeof(uint))
#define UINT_PER_REAL UINT_PER_X(real)
#define UINT_PER_LONG UINT_PER_X(slong)
#define UINT_PER_UNSIGNED UINT_PER_X(unsigned)

// Function for packing tuple_list.  Returns number of uints copied into buffer.
int pack_tuples(TupleList* tl, void **ptr)
{
  uint mi, ml, mul, mr;
  tl->getTupleSize(mi, ml, mul, mr);

  uint n = tl->get_n();

  int sz_buf = 1 + 4*UINT_PER_UNSIGNED +
               tl->get_n() * (mi + 
			      ml*UINT_PER_LONG + 
			      mul*UINT_PER_LONG + 
			      mr*UINT_PER_REAL);
  
  uint *buf = (uint*) malloc(sz_buf*sizeof(uint));
  *ptr = (void*) buf;

  // copy n
  memcpy(buf, &n,   sizeof(uint)),                buf+=1;
  // copy mi
  memcpy(buf, &mi,  sizeof(unsigned)),            buf+=UINT_PER_UNSIGNED;
  // copy ml
  memcpy(buf, &ml,  sizeof(unsigned)),            buf+=UINT_PER_UNSIGNED;
  // copy mul
  memcpy(buf, &mul, sizeof(unsigned)),            buf+=UINT_PER_UNSIGNED;
  // copy mr
  memcpy(buf, &mr,  sizeof(unsigned)),            buf+=UINT_PER_UNSIGNED;
  // copy vi_wr
  memcpy(buf, tl->vi_rd,     tl->get_n()*mi*sizeof(sint)),   buf+=tl->get_n()*mi;
  // copy vl_wr
  memcpy(buf, tl->vl_rd,     tl->get_n()*ml*sizeof(slong)),  buf+=tl->get_n()*ml*UINT_PER_LONG;
  // copy vul_wr
  memcpy(buf, tl->vul_rd,    tl->get_n()*mul*sizeof(ulong)), buf+=tl->get_n()*mul*UINT_PER_LONG;
  // copy vr_wr
  memcpy(buf, tl->vr_rd,     tl->get_n()*mr*sizeof(real)),   buf+=tl->get_n()*mr*UINT_PER_REAL;

  return sz_buf;
}

// Function for packing tuple_list
void unpack_tuples(void *ptr, TupleList** tlp)
{
  TupleList *tl = new TupleList();
  *tlp = tl;

  uint nt;
  unsigned mit, mlt, mult, mrt;
  uint *buf = (uint*)ptr;

  // get n
  memcpy(&nt,   buf, sizeof(uint)),          buf+=1;
  // get mi
  memcpy(&mit,  buf, sizeof(unsigned)),      buf+=UINT_PER_UNSIGNED;
  // get ml
  memcpy(&mlt,  buf, sizeof(unsigned)),      buf+=UINT_PER_UNSIGNED;
  // get mul
  memcpy(&mult, buf, sizeof(unsigned)),      buf+=UINT_PER_UNSIGNED;
  // get mr
  memcpy(&mrt,  buf, sizeof(unsigned)),      buf+=UINT_PER_UNSIGNED;

  // initalize tl
  tl->initialize(mit, mlt, mult, mrt, nt);
  tl->enableWriteAccess();
  tl->set_n( nt );

  uint mi, ml, mul, mr;
  tl->getTupleSize(mi, ml, mul, mr);

  // get vi_rd
  memcpy(tl->vi_wr,     buf, tl->get_n()*mi*sizeof(sint)),   buf+=tl->get_n()*mi;
  // get vl_rd
  memcpy(tl->vl_wr,     buf, tl->get_n()*ml*sizeof(slong)),  buf+=tl->get_n()*ml*UINT_PER_LONG;
  // get vul_rd
  memcpy(tl->vul_wr,    buf, tl->get_n()*mul*sizeof(ulong)), buf+=tl->get_n()*mul*UINT_PER_LONG;
  // get vr_rd
  memcpy(tl->vr_wr,     buf, tl->get_n()*mr*sizeof(real)),   buf+=tl->get_n()*mr*UINT_PER_REAL;

  tl->disableWriteAccess();
  return;
}

ErrorCode Coupler::get_gl_points_on_elements(Range & targ_elems, std::vector<double> & vpos, int & numPointsOfInterest)
{
  numPointsOfInterest = targ_elems.size() * _ntot;//
  vpos.resize(3*numPointsOfInterest);
  int ielem=0;
  for (Range::iterator eit = targ_elems.begin(); eit!=targ_elems.end(); eit++, ielem +=_ntot*3)
  {
    EntityHandle eh = *eit;
    const double * xval;
    const double * yval;
    const double * zval;
    ErrorCode rval = mbImpl-> tag_get_by_ptr(_xm1Tag, &eh, 1,(const void **) &xval );
    if (moab::MB_SUCCESS != rval)
    {
      std::cout << "can't get xm1 values \n";
      return MB_FAILURE;
    }
    rval = mbImpl-> tag_get_by_ptr(_ym1Tag, &eh, 1, (const void **)&yval );
    if (moab::MB_SUCCESS != rval)
    {
      std::cout << "can't get ym1 values \n";
      return MB_FAILURE;
    }
    rval = mbImpl-> tag_get_by_ptr(_zm1Tag, &eh, 1, (const void **)&zval );
    if (moab::MB_SUCCESS != rval)
    {
      std::cout << "can't get zm1 values \n";
      return MB_FAILURE;
    }
    // now, in a stride, populate vpos
    for (int i=0; i<_ntot; i++)
    {
      vpos[ielem + 3*i    ] = xval[i];
      vpos[ielem + 3*i + 1] = yval[i];
      vpos[ielem + 3*i + 2] = zval[i];
    }
  }
  return MB_SUCCESS;
}
} // namespace_moab
