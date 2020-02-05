//
// Usage:
// tools/mbslavepart -d 2 -m mpas/x1.2562.grid.h5m -s mpas/x1.10242.grid.h5m -o mpas_slave.h5m -e 1e-8 -b 1e-6 -O
//
#include <iostream>
#include <exception>
#include <cmath>
#include <vector>
#include <string>

#include "moab/ProgOptions.hpp"
#include "moab/Core.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "moab/BVHTree.hpp"
#include "moab/nanoflann.hpp"

#include "moab/IntxMesh/IntxUtils.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#endif

using namespace moab;

// #define MOAB_USE_NANOFLANN


// And this is the "dataset to kd-tree" adaptor class:
template <typename Derived, typename coord_t>
struct PointCloudAdaptor
{
  //typedef typename Derived::coord_t coord_t;

  const Derived &obj; //!< A const ref to the data set origin

  /// The constructor that sets the data set source
  PointCloudAdaptor(const Derived &obj_) : obj(obj_) { }

  /// CRTP helper method
  inline const Derived& derived() const { return obj; }

  // Must return the number of data points
  inline size_t kdtree_get_point_count() const { return obj.size()/3; }

  // Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
  inline coord_t kdtree_distance_euclidean(const coord_t *p1, const size_t idx_p2,size_t /*size*/) const
  {
    assert(idx_p2 < obj.size());
    size_t offset = idx_p2*3;
    const coord_t d0=p1[0]-obj[offset];
    const coord_t d1=p1[1]-obj[offset+1];
    const coord_t d2=p1[2]-obj[offset+2];
    return d0*d0+d1*d1+d2*d2;
  }

  // Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
  inline coord_t kdtree_distance_angle(const coord_t *p1, const size_t idx_p2,size_t /*size*/) const
  {
    assert(idx_p2 < obj.size());
    size_t offset = idx_p2*3;
    const coord_t d0=p1[0]-obj[offset];
    const coord_t d1=p1[1]-obj[offset+1];
    const coord_t d2=p1[2]-obj[offset+2];
    return atan2((d0*d0+d1*d1+d2*d2),(std::sqrt(p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2])*std::sqrt(obj[offset]*obj[offset]+obj[offset+1]*obj[offset+1]+obj[offset+2]*obj[offset+2])))*180/M_PI;
  }

  inline coord_t kdtree_distance(const coord_t *p1, const size_t idx_p2, size_t size) const
  {
    // return kdtree_distance_euclidean(p1,idx_p2,size);
    return kdtree_distance_angle(p1,idx_p2,size);
  }

  // Returns the dim'th component of the idx'th point in the class:
  // Since this is inlined and the "dim" argument is typically an immediate value, the
  //  "if/else's" are actually solved at compile time.
  inline coord_t kdtree_get_pt(const size_t idx, int dim) const
  {
    assert(idx < obj.size());
    if (dim==0) return obj[idx*3];
    else if (dim==1) return obj[idx*3+1];
    else return obj[idx*3+2];
  }

  // Optional bounding-box computation: return false to default to a standard bbox computation loop.
  //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
  //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
  template <class BBOX>
  bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }

}; // end of PointCloudAdaptor


// A function to get the non-default value from a std::map
template <typename K, typename V>
static V get_map_value(const  std::map <K,V> & m, const K & key, const V & defval ) {
   typename std::map<K,V>::const_iterator it = m.find( key );
   if ( it == m.end() ) {
      return defval;
   }
   else {
      return it->second;
   }
}

int main(int argc, char* argv[])
{
  int proc_id = 0, size = 1, dimension=3;
#ifdef MOAB_HAVE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &proc_id );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
#endif

  int defaultpart=0;
  double tolerance=1e-13, btolerance=1e-4;
  std::string masterfile, slavefile, outfile("slavemesh.h5m");
  bool keepsparts=false;
  bool use_spherical=false;
  ProgOptions opts;

  opts.addOpt<std::string>("master,m", "Master mesh filename", &masterfile);
  opts.addOpt<std::string>("slave,s", "Slave mesh filename", &slavefile);
  opts.addOpt<std::string>("output,o", "Output partitioned mesh filename", &outfile);
  opts.addOpt<int>("dim,d", "Dimension of entities to use for partitioning", &dimension);
  opts.addOpt<int>("defaultpart,p", "Default partition number if target element is not found on source grid", &defaultpart);
  opts.addOpt<double>("eps,e", "Tolerance for the point search", &tolerance);
  opts.addOpt<double>("beps,b", "Tolerance for the bounding box search", &btolerance);
  opts.addOpt<void>("keep,K", "Keep the existing partitions in the slave mesh (use PARALLEL_PARTITION_SLAVE instead)", &keepsparts);
  opts.addOpt<void>("spherical", "Hint that the meshes are defined on a spherical surface (Climate problems)", &use_spherical);
  opts.parseCommandLine(argc, argv);

  if (masterfile.empty() || slavefile.empty())
  {
    opts.printHelp();
#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif
    exit(1);
  }

  ErrorCode error;
  Core *mbCore = new Core();

  //Set the read options for parallel file loading
  const std::string partition_set_name = "PARALLEL_PARTITION";
  const std::string global_id_name = "GLOBAL_ID";
  std::vector<std::string> read_opts, write_opts;
  std::string read_options(""), write_options("");

  if (size > 1) {
    read_options = "PARALLEL=READ_PART;PARTITION="+partition_set_name+";PARALLEL_RESOLVE_SHARED_ENTS";
    write_options = "PARALLEL=WRITE_PART";
  }

  EntityHandle masterfileset, slavefileset;
  error = mbCore->create_meshset(moab::MESHSET_TRACK_OWNER | moab::MESHSET_SET, masterfileset);MB_CHK_ERR(error);
  error = mbCore->create_meshset(moab::MESHSET_TRACK_OWNER | moab::MESHSET_SET, slavefileset);MB_CHK_ERR(error);

  //Load file
  error = mbCore->load_file(masterfile.c_str(), &masterfileset, read_options.c_str());MB_CHK_ERR(error);
  error = mbCore->load_file(slavefile.c_str(), &slavefileset, read_options.c_str());MB_CHK_ERR(error);
  // if (error != MB_SUCCESS && size > 1)
  // {
  //   std::string newread_options = "PARALLEL=BCAST_DELETE;PARALLEL_RESOLVE_SHARED_ENTS";
  //   error = mbCore->load_file(slavefile.c_str(), &slavefileset, newread_options.c_str());
  // }
  // else MB_CHK_ERR(error);

  Tag gidtag=0, parttag=0, sparttag=0;
  int dum_id = -1;
  error = mbCore->tag_get_handle(partition_set_name.c_str(), 1, MB_TYPE_INTEGER, parttag, MB_TAG_SPARSE | MB_TAG_CREAT, &dum_id);MB_CHK_ERR(error);
  gidtag = mbCore->globalId_tag();
  if (keepsparts) {
    error = mbCore->tag_get_handle(std::string(partition_set_name+"_SLAVE").c_str(), 1, MB_TYPE_INTEGER, sparttag, MB_TAG_CREAT | MB_TAG_SPARSE, &dum_id);MB_CHK_ERR(error);
  }

  Range melems, msets, selems, ssets;

  // Get the partition sets on the master mesh
  std::map<int, int> mpartvals;
  error = mbCore->get_entities_by_type_and_tag(masterfileset, MBENTITYSET, &parttag, NULL, 1, msets, moab::Interface::UNION, true);MB_CHK_ERR(error);
  if (msets.size() == 0) {
    std::cout << "No partition sets found in the master mesh. Quitting..." << std::endl;
    exit (1);
  }

  for (unsigned i=0; i < msets.size(); ++i) {
    EntityHandle mset=msets[i];

    moab::Range msetelems;
    error = mbCore->get_entities_by_dimension(mset, dimension, msetelems);MB_CHK_ERR(error);
    melems.merge(msetelems);

    int partID;
    error = mbCore->tag_get_data(parttag, &mset, 1, &partID);MB_CHK_ERR(error);

    // Get the global ID and use that as the indicator
    std::vector<int> gidMelems(msetelems.size());
    error = mbCore->tag_get_data(gidtag, msetelems, gidMelems.data());MB_CHK_ERR(error);

    for (unsigned j=0; j < msetelems.size(); ++j)
      mpartvals[gidMelems[j]]=partID;
      // mpartvals[msetelems[j]]=partID;
    
    std::cout << "Part " << partID << " has " << msetelems.size() << " elements." << std::endl;
  }

  // Get information about the slave file set
  error = mbCore->get_entities_by_type_and_tag(slavefileset, MBENTITYSET, &parttag, NULL, 1, ssets, moab::Interface::UNION);MB_CHK_ERR(error);
  // TODO: expand and add other dimensional elements
  error = mbCore->get_entities_by_dimension(slavefileset, dimension, selems);MB_CHK_ERR(error);

  std::cout << "Master (elements, parts) : (" << melems.size() << ", " << msets.size() << "), Slave elements : (" << selems.size() << ", " << ssets.size() << ")" << std::endl;

  double master_radius = 1.0, slave_radius = 1.0;
  std::vector<double> mastercoords;
  Range masterverts, slaveverts;
  {
    error = mbCore->get_entities_by_dimension(masterfileset, 0, masterverts);MB_CHK_ERR(error);
    error = mbCore->get_entities_by_dimension(slavefileset, 0, slaveverts);MB_CHK_ERR(error);
#ifdef MOAB_USE_NANOFLANN
    // mastercoords.resize(3*masterverts.size());
    // error = mbCore->get_coords(masterverts, &mastercoords[0]);MB_CHK_ERR(error);
    mastercoords.resize(3*melems.size());
    error = mbCore->get_coords(melems, &mastercoords[0]);MB_CHK_ERR(error);
#endif
  }
  if (use_spherical)
  {
    double points[6];
    EntityHandle mfrontback[2] = {masterverts[0], masterverts[masterverts.size()-1]};
    error = mbCore->get_coords(&mfrontback[0], 2, points);MB_CHK_ERR(error);
    master_radius = 0.5*(std::sqrt(points[0]*points[0]+points[1]*points[1]+points[2]*points[2]) + std::sqrt(points[3]*points[3]+points[4]*points[4]+points[5]*points[5]));
    EntityHandle sfrontback[2] = {slaveverts[0], slaveverts[slaveverts.size()-1]};
    error = mbCore->get_coords(&sfrontback[0], 2, points);MB_CHK_ERR(error);
    slave_radius = 0.5*(std::sqrt(points[0]*points[0]+points[1]*points[1]+points[2]*points[2]) + std::sqrt(points[3]*points[3]+points[4]*points[4]+points[5]*points[5]));
    // Let us rescale both master and slave meshes to a unit sphere
    error = ScaleToRadius(mbCore, masterfileset, 1.0);MB_CHK_ERR(error);
    error = ScaleToRadius(mbCore, slavefileset, 1.0);MB_CHK_ERR(error);
  }

  try {
    std::map<int, moab::Range > spartvals;
    int npoints_notfound=0;
    {
      FileOptions fopts((use_spherical ? "SPHERICAL" : ""));
      moab::EntityHandle tree_root = masterfileset;

#ifdef MOAB_USE_NANOFLANN
      typedef PointCloudAdaptor<std::vector<double>, double > PC2KD;

      const PC2KD  pc2kd(mastercoords); // The adaptor

      // construct a kd-tree index:
      typedef nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<double, PC2KD>,
        PC2KD,
        3 /* dim */
        > NanoFlannKDTree;

      NanoFlannKDTree   tree(dimension /*dim*/, pc2kd, nanoflann::KDTreeSingleIndexAdaptorParams(50 /* max leaf */) );
      tree.buildIndex();

      // const double avgarea = 1.0 / melems.size() * (use_spherical ? 2.0 * M_PI : 1.0); // crude estimate of average element surface area in master mesh
      // const double search_radius = static_cast<double>(std::pow(avgarea, 1.0/dimension) * tolerance);
      // const double search_radius = static_cast<double>(1.0);

      /*
      const double search_radius = static_cast<double>(0.25);
      std::vector<std::pair<size_t,double> >   ret_matches;

      nanoflann::SearchParams params;
      params.sorted = true;

      const size_t nMatches = tree.radiusSearch(&query_pt[0], search_radius, ret_matches, params);
      */

#else

      moab::AdaptiveKDTree tree(mbCore);
      error = tree.build_tree (melems, &tree_root, &fopts);

      // moab::BVHTree tree(mbCore);
      // error = tree.build_tree(melems, &tree_root);MB_CHK_ERR(error);

#endif

      for (size_t ie=0; ie < selems.size(); ie++) {
        moab::EntityHandle selem, leaf;
        double point[3];
        selem = selems[ie];

        // Get the element centroid to be queried
        error = mbCore->get_coords(&selem, 1, point);MB_CHK_ERR(error);

        std::vector<moab::EntityHandle> leaf_elems;
#ifdef MOAB_USE_NANOFLANN

        leaf = masterfileset;
        std::vector<std::pair<size_t,double> >   ret_matches;

        // nanoflann::SearchParams params;
        // params.sorted = true;
        // const size_t nMatches = tree.radiusSearch(point, search_radius, ret_matches, params);
        // std::cout << ie << ": " << nMatches << " matches and ret_matches = " << ret_matches.size() << std::endl;


        const size_t num_results = 1;
        size_t ret_index;
        double out_dist_sqr;
        nanoflann::KNNResultSet<double> resultSet(num_results);
        resultSet.init(&ret_index, &out_dist_sqr );
        tree.findNeighbors(resultSet, point, nanoflann::SearchParams(10)); 
        ret_matches.push_back(std::make_pair(ret_index, out_dist_sqr));

        // We only care about the dimension that the user specified.
        // MOAB partitions are ordered by elements anyway.
        for (unsigned ir = 0; ir < ret_matches.size(); ++ir)
        {
          leaf_elems.push_back(melems[ret_matches[ir].first]);
        }

#else

        // Search for the closest source element in the master mesh corresponding
        // to the target element centroid in the slave mesh
        error = tree.point_search( point, leaf, tolerance, btolerance );MB_CHK_ERR(error);

        // We only care about the dimension that the user specified.
        // MOAB partitions are ordered by elements anyway.
        error = mbCore->get_entities_by_dimension( leaf, dimension, leaf_elems, true);MB_CHK_ERR(error);

#endif

        if (leaf != 0 && leaf_elems.size()) {

          // Now get the master element centroids so that we can compute
          // the minimum distance to the target point
          std::vector<double> centroids(leaf_elems.size()*3);
          error = mbCore->get_coords(&leaf_elems[0], leaf_elems.size(), &centroids[0]);MB_CHK_ERR(error);

          if (!leaf_elems.size())
            std::cout << ie << ": " << " No leaf elements found." << std::endl;
          // else std::cout << ie << " found " << leaf_elems.size() << " leaves for current point " << std::endl;

          double dist=1e10;
          int pinelem=-1;
          for (size_t il=0; il < leaf_elems.size(); ++il) {
#ifdef MOAB_USE_NANOFLANN
            const double locdist = ret_matches[il].second;
#else
            const double *centroid = &centroids[il*3];
            const double locdist = std::pow(point[0]-centroid[0],2)+std::pow(point[1]-centroid[1],2)+std::pow(point[2]-centroid[2],2);
#endif
            if (locdist < dist && locdist < 1.0E-1) {
              dist = locdist;
              pinelem = il;

#ifndef VERBOSE
              int gidMelem;
              error = mbCore->tag_get_data(gidtag, &leaf_elems[il], 1, &gidMelem);MB_CHK_ERR(error);
              std::cout << "\t Trial leaf " << il << " set " << gidMelem << " and part = " << get_map_value(mpartvals, gidMelem, -1) << " with distance = " << locdist << std::endl;
#endif
            }
          }

          if (pinelem < 0) {
            std::cout << ie << ": [Error] - Could not find a minimum distance within the leaf nodes." << std::endl;
            npoints_notfound++;
          }
          else {
            int gidMelem;
            error = mbCore->tag_get_data(gidtag, &leaf_elems[pinelem], 1, &gidMelem);MB_CHK_ERR(error);

            int mpartid = get_map_value(mpartvals, gidMelem, -1);
            if (mpartid < 0)
              std::cout << "[WARNING]: Part number for element " << leaf_elems[pinelem] << " with global ID = " << gidMelem << " not found.\n";

#ifdef VERBOSE
            std::cout << "Adding element " << ie << " set " << mpartid << " with distance = " << dist << std::endl;
#endif
            spartvals[ mpartid ].insert(selems[ie]);
          }
        }
        else
        {
          std::cout << "[WARNING]: Adding element " << ie << " to default (part) set " << defaultpart << std::endl;

          spartvals[ defaultpart ].insert(selems[ie]);
        }

      }
#ifndef MOAB_USE_NANOFLANN
      error = tree.reset_tree();MB_CHK_ERR(error);
#endif
    }
    if (npoints_notfound) std::cout << "Could not find " << npoints_notfound << " points in the master mesh" << std::endl;

    if (use_spherical)
    {
      error = ScaleToRadius(mbCore, slavefileset, slave_radius);MB_CHK_ERR(error);
    }

    error = mbCore->delete_entities(&masterfileset, 1);MB_CHK_ERR(error);
    // Find parallel partition sets in the slave mesh - and delete it since we are going to overwrite the sets
    if (!keepsparts) {
      std::cout << "Deleting " << ssets.size() << " sets in the slave mesh" << std::endl;
      error = mbCore->remove_entities(slavefileset, ssets);MB_CHK_ERR(error);
      ssets.clear();
    }

    for (std::map<int, moab::Range >::iterator it = spartvals.begin(); it != spartvals.end(); ++it) {
      int partID = it->first;
      moab::EntityHandle pset;
      error = mbCore->create_meshset(moab::MESHSET_SET, pset);MB_CHK_ERR(error);
      error = mbCore->add_entities(pset, it->second);MB_CHK_ERR(error);
      error = mbCore->add_parent_child(slavefileset, pset);MB_CHK_ERR(error);

      std::cout << "Slave Part " << partID << " has " << it->second.size() << " elements." << std::endl;

      if (keepsparts) {
        error = mbCore->tag_set_data(sparttag, &pset, 1, &partID);MB_CHK_ERR(error);
      }
      else {
        error = mbCore->tag_set_data(parttag, &pset, 1, &partID);MB_CHK_ERR(error);
      }
    }

    // mbCore->print_database();

    // Write the re-partitioned slave mesh to disk
    if (size == 1) {
      error = mbCore->write_file("slavemesh.vtk", "VTK", NULL, &slavefileset, 1);MB_CHK_ERR(error);
    }
    error = mbCore->write_file(outfile.c_str(), NULL, write_options.c_str(), &slavefileset, 1);MB_CHK_ERR(error);
    // error = mbCore->write_file(outfile.c_str(), NULL, write_options.c_str());MB_CHK_ERR(error);
  }
  catch (std::exception & e)
    {
      std::cout << " exception caught during tree initialization " << e.what() << std::endl;
    }
  delete mbCore;

#ifdef MOAB_HAVE_MPI
  MPI_Finalize();
#endif
  exit(0);
}
