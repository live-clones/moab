/*
 * intx_on_sphere_test.cpp
 *
 *  Created on: Oct 3, 2012
 *      Author: iulian
 */
#include <iostream>
#include <sstream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "moab/Core.hpp"
#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#endif
#include "moab/IntxMesh/Intx2MeshOnSphere.hpp"
#include "moab/IntxMesh/IntxUtils.hpp"
#include "TestUtil.hpp"
#include "moab/ProgOptions.hpp"
#include <math.h>

using namespace moab;

int main(int argc, char* argv[])
{

  std::string firstModel, secondModel, outputFile;

  firstModel = TestDir + "/mbcslam/lagrangeHomme.vtk";
  secondModel = TestDir + "/mbcslam/eulerHomme.vtk";

  ProgOptions opts;
  opts.addOpt<std::string>("first,t", "first mesh filename (source)", &firstModel);
  opts.addOpt<std::string>("second,m", "second mesh filename (target)", &secondModel);
  opts.addOpt<std::string>("outputFile,o", "output intersection file", &outputFile);

  double R = 1.; // input
  double epsrel=1.e-12;
  double boxeps=1.e-4;
  outputFile = "intx.h5m";
  opts.addOpt<double>("radius,R", "radius for model intx", &R);
  opts.addOpt<double>("epsilon,e", "relative error in intx", &epsrel);
  opts.addOpt<double>("boxerror,b", "relative error for box boundaries", &boxeps);


  int output_fraction = 0;
  int write_files_rank = 0;
  int brute_force = 0;

  opts.addOpt<int>("outputFraction,f", "output fraction of areas", &output_fraction);
  opts.addOpt<int>("writeFiles,w", "write files of interest", &write_files_rank);
  opts.addOpt<int>("kdtreeOption,k", "use kd tree for intersection", &brute_force);

  opts.parseCommandLine(argc, argv);
  int rank=0, size=1;
#ifdef MOAB_HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif



  // check command line arg second grid is red, arrival, first mesh is blue, departure
  // will will keep the
  std::string optsRead = (size == 1 ? "" : std::string("PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION")+
        std::string(";PARALLEL_RESOLVE_SHARED_ENTS"));

  // read meshes in 2 file sets
  ErrorCode rval;
  Core moab;
  Interface * mb = &moab; // global
  EntityHandle sf1, sf2, outputSet;

  // create meshsets and load files

  rval = mb->create_meshset(MESHSET_SET, sf1);MB_CHK_ERR(rval);
  rval = mb->create_meshset(MESHSET_SET, sf2);MB_CHK_ERR(rval);
  if (0==rank)
    std::cout << "Loading mesh file " <<firstModel << "\n";
  rval = mb->load_file(firstModel.c_str(), &sf1, optsRead.c_str());MB_CHK_ERR(rval);
  if (0==rank)
    std::cout << "Loading mesh file " << secondModel << "\n";
  rval = mb->load_file(secondModel.c_str(), &sf2, optsRead.c_str());MB_CHK_ERR(rval);

  if (0==rank)
  {
    std::cout << "Radius:  " << R << "\n";
    std::cout << "relative eps:  " << epsrel << "\n";
    std::cout << "box eps:  " << boxeps << "\n";
    std::cout << " use kd tree for intersection: " << brute_force << "\n";
  }
  rval = mb->create_meshset(MESHSET_SET, outputSet);MB_CHK_ERR(rval);

  // fix radius of both meshes, to be consistent with input R
  rval = ScaleToRadius(mb, sf1, R); MB_CHK_ERR(rval);
  rval = ScaleToRadius(mb, sf2, R); MB_CHK_ERR(rval);

#if 0
  // std::cout << "Fix orientation etc ..\n";
  //IntxUtils; those calls do nothing for a good mesh
  rval = fix_degenerate_quads(mb, sf1);MB_CHK_ERR(rval);
  rval = fix_degenerate_quads(mb, sf2);MB_CHK_ERR(rval);

  rval = positive_orientation(mb, sf1, R);MB_CHK_ERR(rval);
  rval = positive_orientation(mb, sf2, R);MB_CHK_ERR(rval);
#endif

#ifdef MOAB_HAVE_MPI
  ParallelComm* pcomm = ParallelComm::get_pcomm(mb, 0);
#endif
  Intx2MeshOnSphere  worker(mb);

  worker.set_error_tolerance(R*epsrel);
  worker.set_box_error(boxeps);
#ifdef MOAB_HAVE_MPI
  worker.set_parallel_comm(pcomm);
#endif
  //worker.SetEntityType(moab::MBQUAD);
  worker.set_radius_source_mesh(R);
  worker.set_radius_destination_mesh(R);
  //worker.enable_debug();

  rval = worker.FindMaxEdges(sf1, sf2);MB_CHK_ERR(rval);

  EntityHandle covering_set;
#ifdef MOAB_HAVE_MPI
  if (size>1)
  {
    Range local_verts;
    rval = worker.build_processor_euler_boxes(sf2, local_verts); MB_CHK_ERR(rval);// output also the local_verts
    if (write_files_rank)
    {
      std::stringstream outf;
      outf<<"second_mesh" << rank<<".h5m";
      rval = mb->write_file(outf.str().c_str(), 0, 0, &sf2, 1); MB_CHK_ERR(rval);
    }
  }
  if (size>1)
  {
    double elapsed = MPI_Wtime();
    rval = mb->create_meshset(moab::MESHSET_SET, covering_set);MB_CHK_SET_ERR(rval, "Can't create new set");
    rval = worker.construct_covering_set(sf1, covering_set); MB_CHK_ERR(rval);// lots of communication if mesh is distributed very differently
    elapsed = MPI_Wtime() - elapsed;
    if (0==rank)
      std::cout << "\nTime to communicate the mesh = " << elapsed << std::endl;
    // area fraction of the covering set that needed to be communicated from other processors
    // number of elements in the covering set communicated, compared to total number of elements in the covering set
    if (output_fraction)
    {
      EntityHandle comm_set; // set with elements communicated from other tasks
      rval = mb->create_meshset(MESHSET_SET, comm_set); MB_CHK_ERR(rval);
      // see how much more different is compared to sf1
      rval = mb->unite_meshset(comm_set, covering_set) ;  MB_CHK_ERR(rval); // will have to subtract from covering set, initial set
      // subtract
      rval = mb->subtract_meshset(comm_set, sf1); MB_CHK_ERR(rval);
      // compute fractions
      double area_cov_set = area_on_sphere(mb, covering_set, R);
      assert(area_cov_set>0);
      double comm_area = area_on_sphere(mb, comm_set, R);
      // more important is actually the number of elements communicated
      int num_cov_cells, num_comm_cells;
      rval = mb->get_number_entities_by_dimension(covering_set, 2, num_cov_cells); MB_CHK_ERR(rval);
      rval = mb->get_number_entities_by_dimension(comm_set, 2, num_comm_cells); MB_CHK_ERR(rval);
      double fraction_area = comm_area/area_cov_set;
      double fraction_num_cells = (double) num_comm_cells/num_cov_cells; // determine min, max, average of these fractions

      double max_fraction_area, max_fraction_num_cells, min_fraction_area, min_fraction_num_cells;
      double average_fraction_area, average_fraction_num_cells;
      MPI_Reduce(&fraction_area, &max_fraction_area, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      MPI_Reduce(&fraction_num_cells, &max_fraction_num_cells, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      MPI_Reduce(&fraction_area, &min_fraction_area, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce(&fraction_num_cells, &min_fraction_num_cells, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce(&fraction_area, &average_fraction_area, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&fraction_num_cells, &average_fraction_num_cells, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      average_fraction_area /= size;
      average_fraction_num_cells /= size;
      if (rank==0)
      {
        std::cout << " fraction area:      min: " << min_fraction_area << " max: " << max_fraction_area << " average :" << average_fraction_area << " \n";
        std::cout << " fraction num cells: min: " << min_fraction_num_cells << " max: " << max_fraction_num_cells << " average: " <<average_fraction_num_cells <<  " \n";
      }

    }
    if (write_files_rank)
    {
      std::stringstream cof;
      cof<<"covering_mesh" << rank<<".h5m";
      rval = mb->write_file(cof.str().c_str(), 0, 0, &covering_set, 1); MB_CHK_ERR(rval);
    }
  }
  else
#endif
    covering_set = sf1;

  if (0==rank)
    std::cout << "Computing intersections ..\n";
#ifdef MOAB_HAVE_MPI
  double elapsed = MPI_Wtime();
#endif
  if (brute_force)
  {
    rval = worker.intersect_meshes_kdtree(covering_set, sf2, outputSet);MB_CHK_SET_ERR(rval,"failed to intersect meshes with slow method");
  }
  else
  {
    rval = worker.intersect_meshes(covering_set, sf2, outputSet);MB_CHK_SET_ERR(rval,"failed to intersect meshes");
  }
#ifdef MOAB_HAVE_MPI
  elapsed = MPI_Wtime() - elapsed;
  if (0==rank)
    std::cout << "\nTime to compute the intersection between meshes = " << elapsed << std::endl;
#endif
  // the output set does not have the intx vertices on the boundary shared, so they will be duplicated right now
  // we write this file just for checking it looks OK

  //compute total area with 2 methods
  //double initial_area = area_on_sphere_lHuiller(mb, sf1, R);
  //double area_method1 = area_on_sphere_lHuiller(mb, outputSet, R);
  //double area_method2 = area_on_sphere(mb, outputSet, R);

  //std::cout << "initial area: " << initial_area << "\n";
  //std::cout<< " area with l'Huiller: " << area_method1 << " with Girard: " << area_method2<< "\n";
  //std::cout << " relative difference areas " << fabs(area_method1-area_method2)/area_method1 << "\n";
  //std::cout << " relative error " << fabs(area_method1-initial_area)/area_method1 << "\n";

  if (write_files_rank)
  {
    std::stringstream outf;
    outf<<"intersect" << rank<<".h5m";
    rval = mb->write_file(outf.str().c_str(), 0, 0, &outputSet, 1);
  }
  double intx_area = area_on_sphere(mb, outputSet, R);
  double arrival_area = area_on_sphere(mb, sf2, R) ;
  std::cout<< "On rank : " << rank << " arrival area: " << arrival_area<<
      "  intersection area:" << intx_area << " rel error: " << fabs((intx_area-arrival_area)/arrival_area) << "\n";

#ifdef MOAB_HAVE_MPI
#ifdef MOAB_HAVE_HDF5_PARALLEL
  rval = mb->write_file(outputFile.c_str(), 0, "PARALLEL=WRITE_PART", &outputSet, 1);MB_CHK_SET_ERR(rval,"failed to write intx file");
#else
  // write intx set on rank 0, in serial; we cannot write in parallel
  if (0==rank)
  {
    rval = mb->write_file(outputFile.c_str(), 0, 0, &outputSet, 1);MB_CHK_SET_ERR(rval,"failed to write intx file");
  }
#endif
  MPI_Finalize();
#else
  rval = mb->write_file(outputFile.c_str(), 0, 0, &outputSet, 1);MB_CHK_SET_ERR(rval,"failed to write intx file");
#endif
  return 0;
}

