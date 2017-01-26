/*
 * Usage: MOAB-Tempest tool
 *
 * Generate a Cubed-Sphere mesh: ./mbtempest -t 0 -res 25 -f cubed_sphere_mesh.exo
 * Generate a RLL mesh: ./mbtempest -t 1 -res 25 -f rll_mesh.exo
 * Generate a Icosahedral-Sphere mesh: ./mbtempest -t 2 -res 25 <-dual> -f icosa_mesh.exo
 *
 * Now you can compute the intersections between the meshes too!
 *
 * Generate the overlap mesh: ./mbtempest -t 3 -l cubed_sphere_mesh.exo -l rll_mesh.exo -f overlap_mesh.exo
 *
 */

#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <cassert>

#include "moab/Core.hpp"
#include "moab/TempestRemapper.hpp"
#include "moab/TempestOfflineMap.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/CpuTimer.hpp"
#include "DebugOutput.hpp"

#ifndef MOAB_HAVE_MPI
#error mbtempest tool requires MPI configuration
#endif

// MPI includes
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"


struct ToolContext {
  int blockSize;
  std::vector<std::string> inFilenames;
  std::vector<Mesh*> meshes;
  std::vector<moab::EntityHandle> meshsets;
  std::vector<int> disc_orders;
  std::string outFilename;
  moab::TempestRemapper::TempestMeshType meshType;
  const int proc_id, n_procs;
  bool computeDual;
  bool computeWeights;
  moab::DebugOutput outStream;

  ToolContext(int procid, int nprocs) : 
    blockSize(5), outFilename("output.exo"), meshType(moab::TempestRemapper::DEFAULT), 
    proc_id(procid), n_procs(nprocs),
    computeDual(false), computeWeights(false),
    outStream(std::cout, procid)
  {
    inFilenames.resize(2);
    timer = new moab::CpuTimer();
  }

  ~ToolContext()
  {
    inFilenames.clear();
    outFilename.clear();
    delete timer;
  }

  void clear()
  {
    meshes.clear();
  }

  void timer_push(std::string operation)
  {
    timer_ops = timer->time_since_birth();
    opName = operation;
  }

  void timer_pop()
  {
    outStream.printf(0, "[LOG] Time taken to %s = %f\n", opName.c_str(), timer->time_since_birth() - timer_ops);
    // std::cout << "\n[LOG" << proc_id << "] Time taken to " << opName << " = " << timer->time_since_birth() - timer_ops << std::endl;
    opName.clear();
  }

  void ParseCLOptions(int argc, char* argv[]) {
    ProgOptions opts;
    int imeshType = 0;
    std::string expectedFName = "output.exo";
    int expectedOrder = 1;

    opts.addOpt<int>("res,r", "Resolution of the mesh (default=5)", &blockSize);
    opts.addOpt<int>("type,t", "Type of mesh (default=CS; Choose from [CS=0, RLL=1, ICO=2, OVERLAP_FILES=3, OVERLAP_MEMORY=4, OVERLAP_MOAB=5])", &imeshType);
    opts.addOpt<std::string>("file,f", "Output mesh filename (default=output.exo)", &outFilename);
    opts.addOpt<void>("dual,d", "Output the dual of the mesh (generally relevant only for ICO mesh)", &computeDual);
    opts.addOpt<void>("weights,w", "Compute and output the weights using the overlap mesh (generally relevant only for OVERLAP mesh)", &computeWeights);
    opts.addOpt<std::string>("load,l", "Input mesh filenames (a source and target mesh)", &expectedFName);
    opts.addOpt<int>("order,o", "Discretization orders for the source and target solution fields", &expectedOrder);

    opts.parseCommandLine(argc, argv);

    switch (imeshType) {
    case 0:
      meshType = moab::TempestRemapper::CS;
      break;
    case 1:
      meshType = moab::TempestRemapper::RLL;
      break;
    case 2:
      meshType = moab::TempestRemapper::ICO;
      break;
    case 3:
      meshType = moab::TempestRemapper::OVERLAP_FILES;
      break;
    case 4:
      meshType = moab::TempestRemapper::OVERLAP_MEMORY;
      break;
    case 5:
      meshType = moab::TempestRemapper::OVERLAP_MOAB;
      break;
    default:
      meshType = moab::TempestRemapper::DEFAULT;
      break;
    }

    if (meshType > moab::TempestRemapper::ICO) {
      opts.getOptAllArgs("load,l", inFilenames);
      opts.getOptAllArgs("order,o", disc_orders);
      if (disc_orders.size() == 0)
        disc_orders.resize(2, 1);
      if (disc_orders.size() == 1)
        disc_orders.push_back(1);
      assert(inFilenames.size() == 2);
      assert(disc_orders.size() == 2);
    }
    expectedFName.clear();
  }
private:
  moab::CpuTimer *timer;
  double timer_ops;
  std::string opName;
};

// Forward declare some methods
moab::ErrorCode CreateTempestMesh(ToolContext&, moab::TempestRemapper& remapper, Mesh** );

int main(int argc, char* argv[])
{
  moab::ErrorCode rval;
  NcError error(NcError::verbose_nonfatal);
  std::stringstream sstr;

  int proc_id = 0, nprocs = 1;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &proc_id );
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs );

  ToolContext ctx(proc_id, nprocs);
  ctx.ParseCLOptions(argc, argv);

  moab::Interface* mbCore = new (std::nothrow) moab::Core;
  if (NULL == mbCore) return 1;

  moab::ParallelComm* pcomm = new moab::ParallelComm(mbCore, MPI_COMM_WORLD, 0);

  moab::TempestRemapper remapper (mbCore, pcomm);
  remapper.meshValidate = true;
  remapper.constructEdgeMap = true;
  remapper.initialize();

  Mesh* tempest_mesh = NULL;
  ctx.timer_push("create Tempest mesh");
  rval = CreateTempestMesh(ctx, remapper, &tempest_mesh); MB_CHK_ERR(rval);
  ctx.timer_pop();

  if (ctx.meshType == moab::TempestRemapper::OVERLAP_MEMORY)
  {
    // Compute intersections with MOAB
    // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
    assert(ctx.meshes.size() == 3);
 
    rval = pcomm->check_all_shared_handles();MB_CHK_ERR(rval);

    // Load the meshes and validate
    rval = remapper.ConvertTempestMesh(moab::Remapper::SourceMesh); MB_CHK_ERR(rval);
    rval = remapper.ConvertTempestMesh(moab::Remapper::TargetMesh); MB_CHK_ERR(rval);
    rval = remapper.ConvertTempestMesh(moab::Remapper::IntersectedMesh); MB_CHK_ERR(rval);
    rval = mbCore->write_mesh("tempest_intersection.h5m", &ctx.meshsets[2], 1); MB_CHK_ERR(rval);

    // print verbosely about the problem setting
    {
      moab::Range rintxverts, rintxelems;
      rval = mbCore->get_entities_by_dimension(ctx.meshsets[0], 0, rintxverts); MB_CHK_ERR(rval);
      rval = mbCore->get_entities_by_dimension(ctx.meshsets[0], 2, rintxelems); MB_CHK_ERR(rval);
      ctx.outStream.printf(0, "The red set contains %lu vertices and %lu elements \n", rintxverts.size(), rintxelems.size());

      moab::Range bintxverts, bintxelems;
      rval = mbCore->get_entities_by_dimension(ctx.meshsets[1], 0, bintxverts); MB_CHK_ERR(rval);
      rval = mbCore->get_entities_by_dimension(ctx.meshsets[1], 2, bintxelems); MB_CHK_ERR(rval);
      ctx.outStream.printf(0, "The blue set contains %lu vertices and %lu elements \n", bintxverts.size(), bintxelems.size());
    }

    const double epsrel = 1.e-8;
    const double radius = 1.0 /*2.0*acos(-1.0)*/;
    const double boxeps = 0.1;
    moab::EntityHandle intxset; // == remapper.GetMeshSet(moab::Remapper::IntersectedMesh);

    // Compute intersections with MOAB
    {
      // Create the intersection on the sphere object
      ctx.timer_push("setup the intersector");

      moab::Intx2MeshOnSphere *mbintx = new moab::Intx2MeshOnSphere(mbCore);
      mbintx->SetErrorTolerance(epsrel);
      mbintx->set_box_error(boxeps);
      mbintx->SetRadius(radius);
      mbintx->set_parallel_comm(pcomm);

      rval = mbintx->FindMaxEdges(ctx.meshsets[0], ctx.meshsets[1]);MB_CHK_ERR(rval);

      moab::Range local_verts;
      rval = mbintx->build_processor_euler_boxes(ctx.meshsets[1], local_verts); MB_CHK_ERR(rval);
      
      ctx.timer_pop();
      
      moab::EntityHandle covering_set;
      ctx.timer_push("communicate the mesh");
      rval = mbintx->construct_covering_set(ctx.meshsets[0], covering_set); MB_CHK_ERR(rval);// lots of communication if mesh is distributed very differently
      ctx.timer_pop();

      // Now let's invoke the MOAB intersection algorithm in parallel with a 
      // source and target mesh set representing two different decompositions
      ctx.timer_push("compute intersections with MOAB");
      rval = mbCore->create_meshset(moab::MESHSET_SET, intxset);MB_CHK_SET_ERR(rval, "Can't create new set");
      rval = mbintx->intersect_meshes(covering_set, ctx.meshsets[1], intxset); MB_CHK_SET_ERR(rval, "Can't compute the intersection of meshes on the sphere");
      ctx.timer_pop();

      rval = fix_degenerate_quads(mbCore, ctx.meshsets[2]);MB_CHK_ERR(rval);
      rval = positive_orientation(mbCore, ctx.meshsets[2], radius);MB_CHK_ERR(rval);

      // free the memory
      delete mbintx;
    }

    {
      moab::Range intxelems, intxverts;
      rval = mbCore->get_entities_by_dimension(intxset, 2, intxelems); MB_CHK_ERR(rval);
      rval = mbCore->get_entities_by_dimension(intxset, 0, intxverts,true); MB_CHK_ERR(rval);
      ctx.outStream.printf(0, "The intersection set contains %lu elements and %lu vertices \n", intxelems.size(), intxverts.size());

      double initial_area = area_on_sphere_lHuiller(mbCore, ctx.meshsets[0], radius);
      double area_method1 = area_on_sphere_lHuiller(mbCore, intxset, radius);
      double area_method2 = area_on_sphere(mbCore, intxset, radius);

      ctx.outStream.printf(0, "initial area: %12.10f\n", initial_area);
      ctx.outStream.printf(0, " area with l'Huiller: %12.10f with Girard: %12.10f\n", area_method1, area_method2);
      ctx.outStream.printf(0, " relative difference areas = %12.10e\n", fabs(area_method1 - area_method2) / area_method1);
      ctx.outStream.printf(0, " relative error = %12.10e\n", fabs(area_method1 - initial_area) / area_method1);
    }

    // Write out our computed intersection file
    rval = mbCore->write_mesh("moab_intersection.h5m", &intxset, 1); MB_CHK_ERR(rval);

    if (ctx.computeWeights) {
      ctx.timer_push("compute weights with the Tempest meshes");
      // Call to generate an offline map with the tempest meshes
      OfflineMap* weightMap = GenerateOfflineMapWithMeshes(  NULL, *ctx.meshes[0], *ctx.meshes[1], *ctx.meshes[2],
                                                            "", "",     // std::string strInputMeta, std::string strOutputMeta,
                                                            "fv", "fv", // std::string strInputType, std::string strOutputType,
                                                            ctx.disc_orders[0], ctx.disc_orders[1]  // int nPin=4, int nPout=4,
                                                          );
      ctx.timer_pop();
      weightMap->Write("outWeights.nc");
      delete weightMap;
    }
  }
  else if (ctx.meshType == moab::TempestRemapper::OVERLAP_MOAB)
  {
    const double epsrel = 1.e-8;
    const double radius = 1.0 /*2.0*acos(-1.0)*/;
    const double boxeps = 0.1;
    // Usage: mpiexec -n 2 tools/mbtempest -t 5 -l mycs_2.h5m -l myico_2.h5m -f myoverlap_2.h5m

    rval = pcomm->check_all_shared_handles();MB_CHK_ERR(rval);

    // print verbosely about the problem setting
    {
      moab::Range rintxverts, rintxelems;
      rval = mbCore->get_entities_by_dimension(ctx.meshsets[0], 0, rintxverts); MB_CHK_ERR(rval);
      rval = mbCore->get_entities_by_dimension(ctx.meshsets[0], 2, rintxelems); MB_CHK_ERR(rval);
      rval = fix_degenerate_quads(mbCore, ctx.meshsets[0]);MB_CHK_ERR(rval);
      rval = positive_orientation(mbCore, ctx.meshsets[0], radius);MB_CHK_ERR(rval);
      ctx.outStream.printf(0, "The red set contains %lu vertices and %lu elements \n", rintxverts.size(), rintxelems.size());
      
      moab::Range bintxverts, bintxelems;
      rval = mbCore->get_entities_by_dimension(ctx.meshsets[1], 0, bintxverts); MB_CHK_ERR(rval);
      rval = mbCore->get_entities_by_dimension(ctx.meshsets[1], 2, bintxelems); MB_CHK_ERR(rval);
      rval = fix_degenerate_quads(mbCore, ctx.meshsets[1]);MB_CHK_ERR(rval);
      rval = positive_orientation(mbCore, ctx.meshsets[1], radius);MB_CHK_ERR(rval);
      ctx.outStream.printf(0, "The blue set contains %lu vertices and %lu elements \n", bintxverts.size(), bintxelems.size());

      // Assign fully new and compatible Global identifiers for vertices and elements
      // rval = pcomm->assign_global_ids(ctx.meshsets[0], 2, 1, false, true, false);MB_CHK_ERR(rval);
      // rval = pcomm->assign_global_ids(ctx.meshsets[1], 2, rintxverts.size(), false, true, false);MB_CHK_ERR(rval);
      // do not assign new global ids; will make debugging easier , as we want the same results serial/parallel
      // rval = pcomm->assign_global_ids(0, 2, 1, false, true, false);MB_CHK_ERR(rval);
    }

    // Compute intersections with MOAB
    ctx.timer_push("setup and compute mesh intersections");
    rval = remapper.ComputeOverlapMesh(epsrel);MB_CHK_ERR(rval);
    ctx.timer_pop();

    if (false)
    {
      // Create the intersection on the sphere object
      ctx.timer_push("setup the intersector");

      moab::Intx2MeshOnSphere *mbintx = new moab::Intx2MeshOnSphere(mbCore);
      mbintx->SetErrorTolerance(epsrel);
      mbintx->set_box_error(boxeps);
      mbintx->SetRadius(radius);
      mbintx->set_parallel_comm(pcomm);

      rval = mbintx->FindMaxEdges(ctx.meshsets[0], ctx.meshsets[1]);MB_CHK_ERR(rval);

      moab::Range local_verts;
      rval = mbintx->build_processor_euler_boxes(ctx.meshsets[1], local_verts); MB_CHK_ERR(rval);
      // rval = mbintx->build_processor_euler_boxes(ctx.meshsets[0], local_verts); MB_CHK_ERR(rval);
      // ctx.outStream.printf(0, "--Local verts = %lu\n",local_verts.size());

      ctx.timer_pop();

      moab::EntityHandle covering_set = remapper.GetCoveringSet();
      ctx.timer_push("communicate the mesh");
      rval = mbintx->construct_covering_set(ctx.meshsets[0], covering_set); MB_CHK_ERR(rval);// lots of communication if mesh is distributed very differently
      // rval = mbintx->construct_covering_set(ctx.meshsets[1], covering_set); MB_CHK_ERR(rval);// lots of communication if mesh is distributed very differently
      ctx.timer_pop();

      // Now let's invoke the MOAB intersection algorithm in parallel with a 
      // source and target mesh set representing two different decompositions
      ctx.timer_push("compute intersections with MOAB");
      rval = mbintx->intersect_meshes(covering_set, ctx.meshsets[1], ctx.meshsets[2]); MB_CHK_SET_ERR(rval, "Can't compute the intersection of meshes on the sphere");
      // rval = mbintx->intersect_meshes(ctx.meshsets[0], covering_set, ctx.meshsets[2]); MB_CHK_SET_ERR(rval, "Can't compute the intersection of meshes on the sphere");
      ctx.timer_pop();

      // rval = mbCore->add_entities(ctx.meshsets[2], &ctx.meshsets[0], 2);MB_CHK_ERR(rval);
      // rval = mbCore->add_entities(ctx.meshsets[2], &covering_set, 1);MB_CHK_ERR(rval);

      rval = fix_degenerate_quads(mbCore, ctx.meshsets[2]);MB_CHK_ERR(rval);
      rval = positive_orientation(mbCore, ctx.meshsets[2], radius);MB_CHK_ERR(rval);

      // free the memory
      delete mbintx;

      // Assign fully new and compatible Global identifiers for vertices and elements
      // pcomm->assign_global_ids(ctx.meshsets[2], 2, 1, false);
    }

    {
      // Now let us re-convert the MOAB mesh back to Tempest representation
      rval = remapper.ConvertMeshToTempest(moab::Remapper::IntersectedMesh);MB_CHK_ERR(rval);

      // moab::Tag redPtag,bluePtag;
      // rval = mbCore->tag_get_handle("RedParent", redPtag);MB_CHK_ERR(rval);
      // rval = mbCore->tag_get_handle("BlueParent", bluePtag);MB_CHK_ERR(rval);

      // moab::Range overlapEls, overlapVerts;
      // int locsize[2], globsize[2];
      // rval = mbCore->get_entities_by_dimension(ctx.meshsets[2], 2, overlapEls); MB_CHK_ERR(rval);
      // rval = mbCore->get_entities_by_dimension(ctx.meshsets[2], 0, overlapVerts, true); MB_CHK_ERR(rval);
      // rval = pcomm->filter_pstatus(overlapVerts, PSTATUS_NOT_OWNED, PSTATUS_NOT);MB_CHK_ERR(rval);
      // locsize[0] = overlapEls.size(); locsize[1] = overlapVerts.size();
      // ctx.outStream.printf(0, "The intersection set contains %d elements and %d vertices \n", locsize[0], locsize[1]);
      // MPI_Reduce(&locsize, &globsize, 2, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);
      // if (!proc_id) ctx.outStream.printf(0, "-- Global: Intersection set contains %d elements and %d vertices\n", globsize[0], globsize[1]);

      // // Overlap mesh: mesh[2]
      // ctx.meshes[2]->vecSourceFaceIx.resize(overlapEls.size());
      // ctx.meshes[2]->vecTargetFaceIx.resize(overlapEls.size());
      // ctx.meshes[2]->ConstructEdgeMap();

      // rval = mbCore->

      // rval = mbCore->tag_get_data(redPtag,  overlapEls, &ctx.meshes[2]->vecSourceFaceIx[0]); MB_CHK_ERR(rval);
      // rval = mbCore->tag_get_data(bluePtag, overlapEls, &ctx.meshes[2]->vecTargetFaceIx[0]); MB_CHK_ERR(rval);

      rval = remapper.AssociateSrcTargetInOverlap();MB_CHK_ERR(rval);
      ctx.meshes[2] = remapper.GetMesh(moab::Remapper::IntersectedMesh);
    }

    // Write out our computed intersection file
    // we should not add ini sets to the intersection set
    // rval = mbCore->add_entities(ctx.meshsets[2], &ctx.meshsets[0], 2);MB_CHK_ERR(rval);
    std::stringstream ste;
    ste<<"moab_intersection_" << proc_id<<".h5m";

    // write in serial, 2 meshes
    rval = mbCore->write_file(ste.str().c_str(), 0, 0, &ctx.meshsets[2], 1); MB_CHK_ERR(rval);

    {
      double local_areas[3], global_areas[3]; // Array for Initial area, and through Method 1 and Method 2
      local_areas[0] = area_on_sphere_lHuiller(mbCore, ctx.meshsets[1], radius);
      local_areas[1] = area_on_sphere_lHuiller(mbCore, ctx.meshsets[2], radius);
      local_areas[2] = area_on_sphere(mbCore, ctx.meshsets[2], radius);

      MPI_Allreduce(&local_areas, &global_areas, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      if (!proc_id) {
        ctx.outStream.printf(0, "initial area: %12.10f\n", global_areas[0]);
        ctx.outStream.printf(0, " area with l'Huiller: %12.10f with Girard: %12.10f\n", global_areas[1], global_areas[2]);
        ctx.outStream.printf(0, " relative difference areas = %12.10e\n", fabs(global_areas[1] - global_areas[2]) / global_areas[1]);
        ctx.outStream.printf(0, " relative error = %12.10e\n", fabs(global_areas[1] - global_areas[0]) / global_areas[1]);
      }
    }

    if (ctx.computeWeights) {

      ctx.timer_push("compute weights with the Tempest meshes");

      // Call to generate an offline map with the tempest meshes
      TempestOfflineMap* weightMap = new TempestOfflineMap(&remapper);

      weightMap->GenerateOfflineMap("fv", "fv",          // std::string strInputType, std::string strOutputType,
                                    ctx.disc_orders[0], ctx.disc_orders[1],  // int nPin=4, int nPout=4,
                                    false, 0,            // bool fBubble=false, int fMonotoneTypeID=0,
                                    false, false, false/*false*/ // bool fVolumetric=false, bool fNoConservation=false, bool fNoCheck=false,
                                    // std::string strVariables="", std::string strOutputMap="",
                                    // std::string strInputData="", std::string strOutputData="",
                                    // std::string strNColName="", bool fOutputDouble=false,
                                    // std::string strPreserveVariables="", bool fPreserveAll=false, double dFillValueOverride=0.0
                                                          );
      weightMap->GatherAllToRoot();

      // OfflineMap* weightMap;
      // weightMap = GenerateOfflineMapWithMeshes( NULL, *ctx.meshes[0], *ctx.meshes[1], *ctx.meshes[2],
      //                         "", "",              // std::string strInputMeta, std::string strOutputMeta,
      //                         "fv", "fv",          // std::string strInputType, std::string strOutputType,
      //                         ctx.disc_orders[0], ctx.disc_orders[1],  // int nPin=4, int nPout=4,
      //                         false, 0,            // bool fBubble=false, int fMonotoneTypeID=0,
      //                         false, false, (nprocs>1) // bool fVolumetric=false, bool fNoConservation=false, bool fNoCheck=false,
      //                         // std::string strVariables="", std::string strOutputMap="",
      //                         // std::string strInputData="", std::string strOutputData="",
      //                         // std::string strNColName="", bool fOutputDouble=false,
      //                         // std::string strPreserveVariables="", bool fPreserveAll=false, double dFillValueOverride=0.0
      //                                                     );

      // weightMap->m_vecSourceDimSizes.resize(ctx.meshes[0]->faces.size());
      // weightMap->m_vecTargetDimSizes.resize(ctx.meshes[1]->faces.size());

      rval = remapper.ExchangeGhostWeights(weightMap);MB_CHK_ERR(rval);
      ctx.timer_pop();
      sstr.str("");
      sstr << "outWeights_" << proc_id << ".nc";
      weightMap->Write(sstr.str().c_str());
      delete weightMap;
    }
  }

  // Clean up
  ctx.clear();
  delete mbCore;

  MPI_Finalize();
  exit(0);
}


moab::ErrorCode CreateTempestMesh(ToolContext& ctx, moab::TempestRemapper& remapper, Mesh** tempest_mesh)
{
  moab::ErrorCode rval = moab::MB_SUCCESS;

  if (!ctx.proc_id) ctx.outStream.printf(0, "Creating TempestRemap Mesh object ...\n");
  if (ctx.meshType == moab::TempestRemapper::OVERLAP_FILES) {
    // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
    *tempest_mesh = GenerateOverlapMesh(ctx.inFilenames[0], ctx.inFilenames[1], ctx.outFilename, "exact", true);
    ctx.meshes.push_back(*tempest_mesh);
  }
  else if (ctx.meshType == moab::TempestRemapper::OVERLAP_MEMORY) {
    // Load the meshes and validate
    ctx.meshsets.resize(3);
    ctx.meshes.resize(3);
    ctx.meshsets[0] = remapper.GetMeshSet(moab::Remapper::SourceMesh);
    ctx.meshsets[1] = remapper.GetMeshSet(moab::Remapper::TargetMesh);
    ctx.meshsets[2] = remapper.GetMeshSet(moab::Remapper::IntersectedMesh);

    // First the source
    rval = remapper.LoadMesh(moab::Remapper::SourceMesh, ctx.inFilenames[0], moab::TempestRemapper::DEFAULT); MB_CHK_ERR(rval);
    ctx.meshes[0] = remapper.GetMesh(moab::Remapper::SourceMesh);
    
    // Next the target
    rval = remapper.LoadMesh(moab::Remapper::TargetMesh, ctx.inFilenames[1], moab::TempestRemapper::DEFAULT); MB_CHK_ERR(rval);
    ctx.meshes[1] = remapper.GetMesh(moab::Remapper::TargetMesh);

    // Now let us construct the overlap mesh, by calling TempestRemap interface directly
    // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
    *tempest_mesh = GenerateOverlapWithMeshes(*ctx.meshes[0], *ctx.meshes[1], "" /*ctx.outFilename*/, "exact", false);
    remapper.SetMesh(moab::Remapper::IntersectedMesh, *tempest_mesh);
    // ctx.meshes.push_back(*tempest_mesh);
  }
  else if (ctx.meshType == moab::TempestRemapper::OVERLAP_MOAB) {
    ctx.meshsets.resize(3);
    ctx.meshes.resize(3);
    ctx.meshsets[0] = remapper.GetMeshSet(moab::Remapper::SourceMesh);
    ctx.meshsets[1] = remapper.GetMeshSet(moab::Remapper::TargetMesh);
    ctx.meshsets[2] = remapper.GetMeshSet(moab::Remapper::IntersectedMesh);

    // Load the source mesh and validate
    rval = remapper.LoadNativeMesh(ctx.inFilenames[0], ctx.meshsets[0], 0); MB_CHK_ERR(rval);
    rval = remapper.ConvertMeshToTempest(moab::Remapper::SourceMesh); MB_CHK_ERR(rval);
    ctx.meshes[0] = remapper.GetMesh(moab::Remapper::SourceMesh);

    // Load the target mesh and validate 
    rval = remapper.LoadNativeMesh(ctx.inFilenames[1], ctx.meshsets[1], 0); MB_CHK_ERR(rval);
    rval = remapper.ConvertMeshToTempest(moab::Remapper::TargetMesh); MB_CHK_ERR(rval);
    ctx.meshes[1] = remapper.GetMesh(moab::Remapper::TargetMesh);

    // Set the references for the overlap mesh
    // ctx.meshes[2] = remapper.GetMesh(moab::Remapper::IntersectedMesh);
  }
  else if (ctx.meshType == moab::TempestRemapper::ICO) {
    *tempest_mesh = GenerateICOMesh(ctx.blockSize, ctx.computeDual, ctx.outFilename);
    ctx.meshes.push_back(*tempest_mesh);
  }
  else if (ctx.meshType == moab::TempestRemapper::RLL) {
    *tempest_mesh = GenerateRLLMesh(ctx.blockSize * 2, ctx.blockSize, 0.0, 360.0, -90.0, 90.0, false, ctx.outFilename);
    ctx.meshes.push_back(*tempest_mesh);
  }
  else { // default
    *tempest_mesh = GenerateCSMesh(ctx.blockSize, false, ctx.outFilename);
    ctx.meshes.push_back(*tempest_mesh);
  }

  if (ctx.meshType != moab::TempestRemapper::OVERLAP_MOAB && !*tempest_mesh) {
    std::cout << "Tempest Mesh is not a complete object; Quitting...";
    exit(-1);
  }

  return rval;
}

