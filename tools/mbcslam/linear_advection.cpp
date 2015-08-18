
#include <iostream>
#include <sstream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "Intx2MeshOnSphere.hpp"
#include <math.h>
#include "moab/ParallelComm.hpp"
#include "moab/ProgOptions.hpp"
#include "MBParallelConventions.h"
#include "moab/ReadUtilIface.hpp"
#include "MBTagConventions.hpp"
#include "TestUtil.hpp"
#include "CslamUtils.hpp"

#ifdef MESHDIR
std::string TestDir( STRINGIFY(MESHDIR) );
#else
std::string TestDir(".");
#endif

//std::string file_name("./uniform_30.g");
//std::string file_name("./uniform_120.g");
//std::string file_name("./eulerHomme.vtk");

using namespace moab;

moab::ErrorCode update_density(moab::Interface * mb, moab::EntityHandle euler_set,
    moab::EntityHandle lagr_set, moab::EntityHandle out_set, moab::Tag & rhoTag, 
    moab::Tag & areaTag, moab::Tag & rhoCoefsTag, moab::Tag & weightsTag, moab::Tag & planeTag);

moab::ErrorCode get_departure_grid(moab::Interface * mb, moab::EntityHandle euler_set,
                                   moab::EntityHandle lagr_set, moab::EntityHandle covering_set, int tStep, moab::Range & connecVerts);

double gtol = 1.e-9; // this is for geometry tolerance

double radius = 1.;

bool writeFiles = true;
bool parallelWrite = false;
bool velocity = false;

int numSteps = 200; // number of times with velocity displayed at points
double T = 5;

Intx2MeshOnSphere * pworker = NULL;

int main(int argc, char *argv[]) {

   MPI_Init(&argc, &argv);

  // set up MOAB interface and parallel communication
   moab::Core moab;
   moab::Interface& mb = moab;
   moab::ParallelComm mb_pcomm(&mb, MPI_COMM_WORLD);

   //int rank = mb_pcomm->proc_config().proc_rank();
   int rank = mb_pcomm.proc_config().proc_rank();

  // create meshset
   moab::EntityHandle euler_set;
   moab::ErrorCode rval = mb.create_meshset(MESHSET_SET, euler_set);
   CHECK_ERR(rval);

   //std::stringstream opts;
   //opts << "PARALLEL=READ_PART;PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;GATHER_SET=0;PARTITION_METHOD=TRIVIAL_PARTITION;VARIABLE=";
   //opts << "PARALLEL=READ_PART;PARTITION;PARALLEL_RESOLVE_SHARED_ENTS";
   //rval = mb.load_file(file_name.c_str(), &euler_set, opts.str().c_str());
   std::string fileN = TestDir + "/mbcslam/fine4.h5m";

   rval = mb.load_file(fileN.c_str(), &euler_set);
   CHECK_ERR(rval);

   // Create tag for cell density
    moab::Tag rhoTag = 0;
    rval=mb.tag_get_handle("Density",1,moab::MB_TYPE_DOUBLE,rhoTag,moab::MB_TAG_CREAT|moab::MB_TAG_DENSE);
    CHECK_ERR(rval);

   // Create tag for cell area
    moab::Tag areaTag = 0;
    rval=mb.tag_get_handle("Area",1,moab::MB_TYPE_DOUBLE,areaTag,moab::MB_TAG_CREAT|moab::MB_TAG_DENSE);
    CHECK_ERR(rval);
   // Create tag for cell barycenters in 3D Cartesian space
    moab::Tag barycenterTag = 0;
    rval=mb.tag_get_handle("CellBarycenter",3,moab::MB_TYPE_DOUBLE,barycenterTag,moab::MB_TAG_CREAT|moab::MB_TAG_DENSE);
    CHECK_ERR(rval);
   // Create tag for cell density reconstruction coefficients
    moab::Tag rhoCoefTag = 0;
    rval=mb.tag_get_handle("LinearCoefRho",3,moab::MB_TYPE_DOUBLE,rhoCoefTag,moab::MB_TAG_CREAT|moab::MB_TAG_DENSE);
    CHECK_ERR(rval);
   // Create tag for index of gnomonic plane for each cell
    moab::Tag planeTag = 0;
    rval=mb.tag_get_handle("gnomonicPlane",1,moab::MB_TYPE_INTEGER,planeTag,moab::MB_TAG_CREAT|moab::MB_TAG_DENSE);
    CHECK_ERR(rval);
   // Create tag for intersection weights
    moab::Tag weightsTag = 0;
    rval=mb.tag_get_handle("Weights",3,moab::MB_TYPE_DOUBLE,weightsTag,moab::MB_TAG_CREAT|moab::MB_TAG_DENSE);
    CHECK_ERR(rval);
   // get cell plane
    get_gnomonic_plane(&mb, euler_set, planeTag);

   // get cell barycenters (and cell area)
    get_barycenters(&mb, euler_set, planeTag, areaTag, barycenterTag);

   // Set density distributions
    set_density(&mb,euler_set, barycenterTag, rhoTag, 1);

   // Get initial values for use in error computation
    moab::Range redEls;
    rval = mb.get_entities_by_dimension(euler_set, 2, redEls);
    CHECK_ERR(rval);
    std::vector<double> iniValsRho(redEls.size());
    rval = mb.tag_get_data(rhoTag, redEls, &iniValsRho[0]);
    CHECK_ERR(rval);
   // Get Lagrangian set
    moab::EntityHandle out_set, lagrange_set, covering_set;  
    rval = mb.create_meshset(MESHSET_SET, out_set);
    CHECK_ERR(rval);
    rval = mb.create_meshset(MESHSET_SET, lagrange_set);
    CHECK_ERR(rval);
    rval = mb.create_meshset(MESHSET_SET, covering_set);
    CHECK_ERR(rval);
    rval = deep_copy_set(&mb, euler_set, lagrange_set);
    CHECK_ERR(rval);
    moab::EntityHandle dum = 0;
    moab::Tag corrTag;
    rval = mb.tag_get_handle(CORRTAGNAME,1, MB_TYPE_HANDLE, corrTag,
                                           MB_TAG_DENSE|MB_TAG_CREAT, &dum);
    CHECK_ERR(rval);
   //Set up intersection of two meshes

/*
    moab::Intx2MeshOnSphere worker(&mb);
    worker.SetRadius(radius);
    worker.SetErrorTolerance(gtol);
*/


    // pworker is global here; maybe we should pass it around
    pworker = new Intx2MeshOnSphere(&mb);
    pworker->SetErrorTolerance(gtol);
    pworker->SetRadius(radius);
    pworker->set_box_error(100*gtol);


    // these stay fixed for one run
    moab::Range local_verts;
    rval = pworker->build_processor_euler_boxes(euler_set, local_verts);// output also the local_verts
    CHECK_ERR(rval);
    //rval = worker.build_processor_euler_boxes(euler_set, local_verts);// output also the local_verts

    // loop over time to update density
    for (int ts=1; ts < numSteps + 1; ts++){
 
       if (ts  == 1)  // output initial condition
       {
           std::stringstream newDensity;
           newDensity << "Density" << rank << "_" << ts-1 << ".vtk";
           rval = mb.write_file(newDensity.str().c_str(), 0, 0, &euler_set, 1);
       }

      // get linear reconstruction coefficients
       get_linear_reconstruction(&mb, euler_set, rhoTag, planeTag, barycenterTag, rhoCoefTag);

      // get depature grid
       rval =  get_departure_grid(&mb, euler_set, lagrange_set, covering_set, ts, local_verts);
       CHECK_ERR(rval);

      // intersect the meshes
       rval = pworker->intersect_meshes(lagrange_set, euler_set, out_set);
       CHECK_ERR(rval);

      // intersection weights (i.e. area, x integral, and y integral over cell intersections)
       get_intersection_weights3(&mb, euler_set, lagrange_set, out_set, planeTag, weightsTag);

      // update the density
       rval = update_density(&mb, euler_set, lagrange_set, out_set, rhoTag, areaTag,
                              rhoCoefTag, weightsTag, planeTag);
       CHECK_ERR(rval);

       if (writeFiles && (ts % 5 == 0)) // so if write
       {
           std::stringstream newDensity;
           newDensity << "Density" << rank << "_" << ts << ".vtk";
           rval = mb.write_file(newDensity.str().c_str(), 0, 0, &euler_set, 1);
           CHECK_ERR(rval);
       }

       // delete the polygons and elements of out_set
        moab::Range allVerts;
        rval = mb.get_entities_by_dimension(0, 0, allVerts);
        CHECK_ERR(rval);

        moab::Range allElems;
        rval = mb.get_entities_by_dimension(0, 2, allElems);
        CHECK_ERR(rval);

        // get Eulerian and lagrangian cells
        moab::Range polys;
        rval = mb.get_entities_by_dimension(euler_set, 2, polys); 
        CHECK_ERR(rval);
        rval = mb.get_entities_by_dimension(lagrange_set, 2, polys); // do not delete lagr set either, with its vertices
        CHECK_ERR(rval);

        // add to the connecVerts range all verts, from all initial polys
        moab::Range vertsToStay;
        rval = mb.get_connectivity(polys, vertsToStay);
        CHECK_ERR(rval);

        moab::Range todeleteVerts = subtract(allVerts, vertsToStay);

        moab::Range todeleteElem = subtract(allElems, polys);
       // empty the out mesh set
        rval = mb.clear_meshset(&out_set, 1);
        CHECK_ERR(rval);

        rval = mb.delete_entities(todeleteElem);
        CHECK_ERR(rval);
        rval = mb.delete_entities(todeleteVerts);
        CHECK_ERR(rval);
        if (rank==0)
            std::cout << " step: " << ts << "\n";

    }

    //final vals and errors
    moab::Range::iterator iter = redEls.begin();
    double norm1 = 0.;
    double norm2 = 0.;
    double exact2 = 0.;
    double exact1 = 0.;
    int count =0;
    void * data;
    int j=0;// index in iniVals
    while (iter != redEls.end())
    {
       rval = mb.tag_iterate(rhoTag, iter, redEls.end(), count, data);
       CHECK_ERR(rval);
       double * ptrTracer=(double*)data;

       rval = mb.tag_iterate(areaTag, iter, redEls.end(), count, data);
       CHECK_ERR(rval);
       double * ptrArea=(double*)data;
       for (int i=0; i<count; i++, ++iter, ptrTracer++, ptrArea++, j++)
       {
          //double area = *ptrArea;
          norm1+=fabs(*ptrTracer - iniValsRho[j])* (*ptrArea);
          norm2+=(*ptrTracer - iniValsRho[j])*(*ptrTracer - iniValsRho[j])* (*ptrArea);
          exact1+=(iniValsRho[j])* (*ptrArea);
          exact2+=(iniValsRho[j])*(iniValsRho[j])* (*ptrArea);
       }
    }

   double total_norm1=0;
   double total_norm2=0;
   double total_exact1=0;
   double total_exact2=0;
   int mpi_err = MPI_Reduce(&norm1, &total_norm1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   if (mpi_err)
     std::cout <<" error in MPI_reduce:" << mpi_err << "\n";
   mpi_err = MPI_Reduce(&norm2, &total_norm2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   if (mpi_err)
     std::cout <<" error in MPI_reduce:" << mpi_err << "\n";
   mpi_err = MPI_Reduce(&exact1, &total_exact1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   if (mpi_err)
     std::cout <<" error in MPI_reduce:" << mpi_err << "\n";
   mpi_err = MPI_Reduce(&exact2, &total_exact2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   if (mpi_err)
     std::cout <<" error in MPI_reduce:" << mpi_err << "\n";
   if (0==rank)
     std::cout << " numSteps:" << numSteps << " 1-norm:" << total_norm1/total_exact1 << " 2-norm:" << total_norm2/total_exact2 << "\n";

   MPI_Finalize();
   return 0;
}

moab::ErrorCode get_departure_grid(moab::Interface * mb, moab::EntityHandle euler_set,
    moab::EntityHandle lagr_set, moab::EntityHandle covering_set, int tStep, Range & connecVerts)
{
  EntityHandle dum=0;
  Tag corrTag;
  mb->tag_get_handle(CORRTAGNAME, 1, MB_TYPE_HANDLE, corrTag,
                                             MB_TAG_DENSE, &dum);

  double t = tStep * T / numSteps; // numSteps is global; so is T
  double delta_t = T / numSteps; // this is global too, actually
  // double delta_t = 0.0001;
  // double t = delta_t;

  Range polys;
  ErrorCode rval = mb->get_entities_by_dimension(euler_set, 2, polys);

  // change coordinates of lagr mesh vertices
  for (Range::iterator vit = connecVerts.begin(); vit != connecVerts.end();
      ++vit)
  {
    moab::EntityHandle oldV = *vit;
    CartVect posi;
    rval = mb->get_coords(&oldV, 1, &(posi[0]));
    if (MB_SUCCESS != rval)
      return rval;
    // cslam utils, case 1
    CartVect newPos;
    departure_point_swirl_rot(posi, t, T, delta_t, newPos);
    newPos = radius * newPos; // do we need this? the radius should be 1
    moab::EntityHandle new_vert;
    rval = mb->tag_get_data(corrTag, &oldV, 1, &new_vert);
    if (MB_SUCCESS != rval)
      return rval;
    // set the new position for the new vertex
    rval = mb->set_coords(&new_vert, 1, &(newPos[0]));
    if (MB_SUCCESS != rval)
      return rval;
  }

  // if in parallel, we have to move some elements to another proc, and receive other cells
  // from other procs
  rval = pworker->create_departure_mesh_3rd_alg(lagr_set, covering_set);


  return rval;
}

// !!! For now serial !!!
moab::ErrorCode update_density(moab::Interface * mb, moab::EntityHandle euler_set,
    moab::EntityHandle lagr_set, moab::EntityHandle out_set, moab::Tag & rhoTag,
    moab::Tag & areaTag, moab::Tag & rhoCoefsTag, moab::Tag & weightsTag,
    moab::Tag & planeTag)
{

//  moab::ParallelComm * parcomm = ParallelComm::get_pcomm(mb, 0);

  double R = 1.0;

  moab::EntityHandle dum=0;
  Tag corrTag;
  mb->tag_get_handle(CORRTAGNAME, 1, MB_TYPE_HANDLE, corrTag,
                                             MB_TAG_DENSE, &dum);

  moab::Tag gid;
  ErrorCode rval = mb->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, gid, MB_TAG_DENSE);
    if (MB_SUCCESS != rval)
      return rval;

  // get all polygons out of out_set; then see where are they coming from
  moab::Range polys;
  rval = mb->get_entities_by_dimension(out_set, 2, polys);
    if (MB_SUCCESS != rval)
      return rval;

  // get all Lagrangian cells 
  moab::Range rs1;
  rval = mb->get_entities_by_dimension(lagr_set, 2, rs1);
    if (MB_SUCCESS != rval)
      return rval;

  // get all Eulerian cells 
  moab::Range rs2;
  rval = mb->get_entities_by_dimension(euler_set, 2, rs2);
    if (MB_SUCCESS != rval)
      return rval;

  // get gnomonic plane for Eulerian cells
  std::vector<int>  plane(rs2.size());
  rval = mb->tag_get_data(planeTag, rs2, &plane[0]);
    if (MB_SUCCESS != rval)
      return rval;

  // get Eulerian cell reconstruction coefficients
  std::vector<double>  rhoCoefs(3*rs2.size());
  rval = mb->tag_get_data(rhoCoefsTag, rs2, &rhoCoefs[0]);
    if (MB_SUCCESS != rval)
      return rval;

  // get intersection weights 
  std::vector<double>  weights(3*polys.size());
  rval = mb->tag_get_data(weightsTag, polys, &weights[0]);
    if (MB_SUCCESS != rval)
      return rval;

  // Initialize the new values
  std::vector<double> newValues(rs2.size(), 0.);// initialize with 0 all of them

  // For each polygon get red/blue parent
  moab::Tag redParentTag;
  moab::Tag blueParentTag;
  rval = mb->tag_get_handle("RedParent", 1, MB_TYPE_INTEGER, redParentTag, MB_TAG_DENSE);
    if (MB_SUCCESS != rval)
      return rval;
  rval = mb->tag_get_handle("BlueParent", 1, MB_TYPE_INTEGER, blueParentTag, MB_TAG_DENSE);
    if (MB_SUCCESS != rval)
      return rval;

  // mass_lagr = (\sum_intx \int rho^h(x,y) dV)
  // rho_eul^n+1 = mass_lagr/area_eul
  double check_intx_area = 0.;
  int polyIndex = 0;
  for (Range::iterator it= polys.begin(); it!=polys.end(); ++it)
  {

    moab::EntityHandle poly=*it;
    int blueIndex, redIndex;
    rval =  mb->tag_get_data(blueParentTag, &poly, 1, &blueIndex);
    if (MB_SUCCESS != rval)
      return rval;
    moab::EntityHandle blue = rs1[blueIndex];
       
    rval = mb->tag_get_data(redParentTag, &poly, 1, &redIndex);
    if (MB_SUCCESS != rval)
      return rval;

    moab::EntityHandle redArr;
    rval = mb->tag_get_data(corrTag, &blue, 1, &redArr);
    if (MB_SUCCESS != rval)
      return rval;
    int arrRedIndex = rs2.index(redArr);

    // sum into new density values
    newValues[arrRedIndex] += rhoCoefs[redIndex*3]*weights[polyIndex*3] + rhoCoefs[redIndex*3+1]*weights[polyIndex*3+1] 
                                + rhoCoefs[redIndex*3+2]*weights[polyIndex*3+2];

    check_intx_area += weights[polyIndex*3+2];

    polyIndex++;

  }


 // now divide by red area (current)
  int j=0;
  Range::iterator iter = rs2.begin();
  void * data=NULL; //used for stored area
  int count =0;
  double total_mass_local=0.;
  while (iter != rs2.end())
  {
    rval = mb->tag_iterate(areaTag, iter, rs2.end(), count, data);
    if (MB_SUCCESS != rval)
      return rval;
    double * ptrArea=(double*)data;
    for (int i=0; i<count; i++, ++iter, j++, ptrArea++)
    {
      total_mass_local+=newValues[j];
      newValues[j]/= (*ptrArea);
    }
  }

  rval = mb->tag_set_data(rhoTag, rs2, &newValues[0]);
  if (MB_SUCCESS != rval)
    return rval;

  std::cout <<"total mass now:" << total_mass_local << "\n";
  std::cout <<"check: total intersection area: (4 * M_PI * R^2): "  << 4 * M_PI * R*R << " " << check_intx_area << "\n";

  return MB_SUCCESS;
}


<<<<<<< HEAD

/*
 *  Deformational flow 
 */
void departure_point_swirl(moab::CartVect & arrival_point, double t, double delta_t, moab::CartVect & departure_point)
{

  // always assume radius is 1 here?
  moab::SphereCoords sph = cart_to_spherical(arrival_point);
  double k = 2.4; //flow parameter
  /*     radius needs to be within some range   */
  double  sl2 = sin(sph.lon/2);
  double pit = M_PI * t / T;
  double omega = M_PI/T;
  double costheta = cos(sph.lat);
  //double u = k * sl2*sl2 * sin(2*sph.lat) * cos(pit);
  double v = k * sin(sph.lon) * costheta * cos(pit);
  //double psi = k * sl2 * sl2 *costheta * costheta * cos(pit);
  double u_tilda = 2*k*sl2*sl2*sin(sph.lat)*cos(pit);

  // formula 35, page 8
  // this will approximate dep point using a Taylor series with up to second derivative
  // this will be O(delta_t^3) exact.
  double lon_dep = sph.lon - delta_t * u_tilda -delta_t*delta_t * k * sl2 *
      ( sl2 * sin (sph.lat) * sin(pit) * omega
          - u_tilda * sin(sph.lat) * cos(pit) * cos (sph.lon/2)
          - v * sl2 * costheta * cos(pit)   );
  // formula 36, page 8 again
  double lat_dep = sph.lat - delta_t*v - delta_t * delta_t/4* k *
      ( sin(sph.lon)* cos(sph.lat) * sin(pit) * omega
          - u_tilda * cos(sph.lon) * cos(sph.lat) * cos(pit)
          + v * sin(sph.lon) * sin(sph.lat) * cos(pit)  );
  moab::SphereCoords sph_dep;
  sph_dep.R = 1.; // radius
  sph_dep.lat = lat_dep;
  sph_dep.lon = lon_dep;

  departure_point = spherical_to_cart(sph_dep);
  return;
}

/*
 *  Deformational flow with rotation
 */
void departure_point_swirl_rot(moab::CartVect & arrival_point, double t, double delta_t, moab::CartVect & departure_point)
{

  moab::SphereCoords sph = cart_to_spherical(arrival_point);
  double omega = M_PI/T;
  double gt = cos(M_PI*t/T);

  double lambda = sph.lon - 2.0*omega*t;
  double u_tilda = 4.0*sin(lambda)*sin(lambda)*sin(sph.lat)*gt + 2.0*omega;
  double v = 2.0*sin(2.0*lambda)*cos(sph.lat)*gt;

  double lon_dep = sph.lon - delta_t*u_tilda - delta_t*delta_t*2.0*sin(lambda) *
                   (   sin(lambda)*sin(sph.lat)*sin(omega*t)*omega
                     - sin(lambda)*cos(sph.lat)*cos(omega*t)*v
                     - 2.0*cos(lambda)*sin(sph.lat)*cos(omega*t)*u_tilda);

  double lat_dep = sph.lat - delta_t*v - delta_t*delta_t*2.0*
                   (  cos(sph.lat)*sin(omega*t)*omega*sin(lambda)*cos(lambda)
                    - 2.0*u_tilda*cos(sph.lat)*cos(omega*t)*cos(lambda)*cos(lambda)
                    + u_tilda*cos(sph.lat)*cos(omega*t) 
                    + v*sin(sph.lat)*cos(omega*t)*sin(lambda)*cos(lambda));

  moab::SphereCoords sph_dep;
  sph_dep.R = 1.; // radius
  sph_dep.lat = lat_dep;
  sph_dep.lon = lon_dep;

  departure_point = spherical_to_cart(sph_dep);
  return;
}

/*
 *  Zonal flow
 */

void get_intersection_weights(moab::Interface * mb, moab::EntityHandle euler_set, moab::EntityHandle lagr_set, 
                              moab::EntityHandle intx_set, moab::Tag &planeTag, moab::Tag &weightsTag)
{
  // get all intersection polygons
  moab::Range polys;
  ErrorCode rval = mb->get_entities_by_dimension(intx_set, 2, polys);
  if (MB_SUCCESS != rval)
    return;

  // get all Eulerian cells
  moab::Range eul_cells;
  rval = mb->get_entities_by_dimension(euler_set, 2, eul_cells);
  if (MB_SUCCESS != rval)
    return;

  // get all Lagrangian cells
  moab::Range lagr_cells;
  rval = mb->get_entities_by_dimension(lagr_set, 2, lagr_cells);
  if (MB_SUCCESS != rval)
    return;

  // get tag for Eulerian parent cell of intersection polygon
  moab::Tag redParentTag;
  rval = mb->tag_get_handle("RedParent", 1, MB_TYPE_INTEGER, redParentTag, MB_TAG_DENSE);
  if (MB_SUCCESS != rval)
    return;

  // get tag for Lagrangian parent cell of intersection polygon
  moab::Tag blueParentTag;
  rval = mb->tag_get_handle("BlueParent", 1, MB_TYPE_INTEGER, blueParentTag, MB_TAG_DENSE);
  if (MB_SUCCESS != rval)
    return;

  // get gnomonic plane for Eulerian cells
  std::vector<int>  plane(eul_cells.size());
  rval = mb->tag_get_data(planeTag, eul_cells, &plane[0]);
    if (MB_SUCCESS != rval)
      return;

  double total_area = 0.;
  for (moab::Range::iterator it = polys.begin(); it != polys.end(); ++it)
  {
    moab::EntityHandle poly = *it;

    // get the nodes 
    const moab::EntityHandle * verts;
    int num_nodes;
    rval = mb->get_connectivity(poly, verts, num_nodes);
    if (MB_SUCCESS != rval)
      return;

    // get coordinates
    std::vector<double> coords(3 * num_nodes);
    rval = mb->get_coords(verts, num_nodes, &coords[0]);
    if (MB_SUCCESS != rval)
      return;

    // get index of Eulerian parent cell for polygon
    int redIndex;
    rval = mb->tag_get_data(redParentTag, &poly, 1, &redIndex);
    if (MB_SUCCESS != rval)
      return;

    // get index of Lagrangian parent cell for polygon
    int blueIndex;
    rval = mb->tag_get_data(blueParentTag, &poly, 1, &blueIndex);
    if (MB_SUCCESS != rval)
      return;

    std::vector<double> x(num_nodes);
    std::vector<double> y(num_nodes);
    double poly_area = 0;
    double poly_intx = 0;
    double poly_inty = 0;
    double R = 1.0;
    for (int inode = 0; inode < num_nodes; inode++){
         double rad = sqrt(coords[inode*3]*coords[inode*3] + coords[inode*3+1]*coords[inode*3+1] + coords[inode*3+2]*coords[inode*3+2]);
         moab::CartVect xyzcoord(coords[inode*3]/rad,coords[inode*3+1]/rad,coords[inode*3+2]/rad);
         gnomonic_projection(xyzcoord, R, plane[redIndex], x[inode],y[inode]);
    }

    std::vector<double> weights(3);
    for (int inode = 0; inode < num_nodes; inode++){
        int inode2 = inode+1;
        if (inode2 >= num_nodes) inode2 = 0;
        double xmid = 0.5*(x[inode] + x[inode2]);
        double ymid = 0.5*(y[inode] + y[inode2]);
        double r1 = sqrt(1 + x[inode]*x[inode] + y[inode]*y[inode]);
        double rm = sqrt(1 + xmid*xmid + ymid*ymid);
        double r2 = sqrt(1 + x[inode2]*x[inode2] + y[inode2]*y[inode2]);
        double hx=x[inode2]-x[inode];

        poly_area += -hx*(y[inode]/(r1*(1 + x[inode]*x[inode]))
                            + 4.0*ymid/(rm*(1+xmid*xmid))
                            + y[inode2]/(r2*(1 + x[inode2]*x[inode2])))/6.0;

        poly_intx += -hx*(x[inode]*y[inode]/(r1*(1 + x[inode]*x[inode]))
                            + 4.0*xmid*ymid/(rm*(1+xmid*xmid))
                            + x[inode2]*y[inode2]/(r2*(1 + x[inode2]*x[inode2])))/6.0;

        poly_inty += hx*(1.0/r1 + 4.0/rm + 1.0/r2)/6.0;
    }
     weights[0] = poly_intx;
     weights[1] = poly_inty;
     weights[2] = poly_area;

     total_area += poly_area;

     rval = mb->tag_set_data(weightsTag, &poly, 1, &weights[0]);
     if (MB_SUCCESS != rval)
       return;

  }

      std::cout << "polygon area = " << total_area << "\n";
  return;
}

ErrorCode  create_lagr_mesh(moab::Interface * mb, moab::EntityHandle euler_set, moab::EntityHandle lagr_set)
{
  // create the handle tag for the corresponding element / vertex

  moab::EntityHandle dum = 0;

  moab::Tag corrTag;
  ErrorCode rval = mb->tag_get_handle(CORRTAGNAME,
                                           1, MB_TYPE_HANDLE, corrTag,
                                           MB_TAG_DENSE|MB_TAG_CREAT, &dum);

  CHECK_ERR(rval);
  // give the same global id to new verts and cells created in the lagr(departure) mesh
  moab::Tag gid;
  rval = mb->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, gid, MB_TAG_DENSE);
  CHECK_ERR(rval);

  moab::Range polys;
  rval = mb->get_entities_by_dimension(euler_set, 2, polys);
  CHECK_ERR(rval);

  moab::Range connecVerts;
  rval = mb->get_connectivity(polys, connecVerts);
  CHECK_ERR(rval);

  std::map<moab::EntityHandle, moab::EntityHandle> newNodes;
  for (Range::iterator vit = connecVerts.begin(); vit != connecVerts.end(); ++vit)
  {
    moab::EntityHandle oldV = *vit;
    moab::CartVect posi;
    rval = mb->get_coords(&oldV, 1, &(posi[0]));
    CHECK_ERR(rval);
    int global_id;
    rval = mb->tag_get_data(gid, &oldV, 1, &global_id);
    CHECK_ERR(rval);
    moab::EntityHandle new_vert;
    rval = mb->create_vertex(&(posi[0]), new_vert); // duplicate the position
    CHECK_ERR(rval);
    newNodes[oldV] = new_vert;
    // set also the correspondent tag :)
    rval = mb->tag_set_data(corrTag, &oldV, 1, &new_vert);
    CHECK_ERR(rval);
    // also the other side
    // need to check if we really need this; the new vertex will never need the old vertex
    // we have the global id which is the same
    /*rval = mb->tag_set_data(corrTag, &new_vert, 1, &oldV);
    CHECK_ERR(rval);*/
    // set the global id on the corresponding vertex the same as the initial vertex
    rval = mb->tag_set_data(gid, &new_vert, 1, &global_id);
    CHECK_ERR(rval);

  }
 for (Range::iterator it = polys.begin(); it != polys.end(); ++it)
  {
    moab::EntityHandle q = *it;
    int nnodes;
    const moab::EntityHandle * conn;
    rval = mb->get_connectivity(q, conn, nnodes);
    CHECK_ERR(rval);
    int global_id;
    rval = mb->tag_get_data(gid, &q, 1, &global_id);
    CHECK_ERR(rval);
    EntityType typeElem = mb->type_from_handle(q);
    std::vector<moab::EntityHandle> new_conn(nnodes);
    for (int i = 0; i < nnodes; i++)
    {
      moab::EntityHandle v1 = conn[i];
      new_conn[i] = newNodes[v1];
    }
    moab::EntityHandle newElement;
    rval = mb->create_element(typeElem, &new_conn[0], nnodes, newElement);
    CHECK_ERR(rval);
    //set the corresponding tag; not sure we need this one, from old to new
    /*rval = mb->tag_set_data(corrTag, &q, 1, &newElement);
    CHECK_ERR(rval);*/
    rval = mb->tag_set_data(corrTag, &newElement, 1, &q);
    CHECK_ERR(rval);
    // set the global id
    rval = mb->tag_set_data(gid, &newElement, 1, &global_id);
    CHECK_ERR(rval);

    rval = mb->add_entities(lagr_set, &newElement, 1);
    CHECK_ERR(rval);
  }

  return MB_SUCCESS;
}

=======
>>>>>>> translate the linear_tracer branch
