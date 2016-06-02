
#include "Mesquite.hpp"
#include "MsqIBase.hpp"
#include "MsqIGeom.hpp"
#include "MsqIMesh.hpp"
#include "MBiMesh.hpp"
#include "MeshImpl.hpp"
#include "moab/Core.hpp"
#include "moab/Skinner.hpp"
#include "FacetModifyEngine.hpp"
#include "MsqError.hpp"
#include "InstructionQueue.hpp"
#include "PatchData.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"
#include "SphericalDomain.hpp"
#include "PlanarDomain.hpp"
#include "MeshWriter.hpp"

// algorithms
#include "IdealWeightInverseMeanRatio.hpp"
#include "ConditionNumberQualityMetric.hpp"
#include "VertexConditionNumberQualityMetric.hpp"
#include "EdgeLengthQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "PMeanPTemplate.hpp"
#include "SteepestDescent.hpp"
#include "FeasibleNewton.hpp"
#include "ConjugateGradient.hpp"
#include "SmartLaplacianSmoother.hpp"


#include "ReferenceMesh.hpp"
#include "RefMeshTargetCalculator.hpp"
#include "TShapeB1.hpp"
#include "TQualityMetric.hpp"
#include "IdealShapeTarget.hpp"

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

#include "iBase.h"

using namespace MBMesquite;

static int print_iGeom_error( const char* desc,
                         int err,
                         iGeom_Instance geom,
                         const char* file,
                         int line )
{
  char buffer[1024];
  iGeom_getDescription( geom, buffer, sizeof(buffer) );
  buffer[sizeof(buffer)-1] = '\0';

  std::cerr << "ERROR: " << desc << std::endl
            << "  Error code: " << err << std::endl
            << "  Error desc: " << buffer << std::endl
            << "  At        : " << file << ':' << line << std::endl
            ;

  return -1; // must always return false or CHECK macro will break
}

static int print_iMesh_error( const char* desc,
                         int err,
                         iMesh_Instance mesh,
                         const char* file,
                         int line )
{
  char buffer[1024];
  iMesh_getDescription( mesh, buffer, sizeof(buffer) );
  buffer[sizeof(buffer)-1] = '\0';

  std::cerr << "ERROR: " << desc << std::endl
            << "  Error code: " << err << std::endl
            << "  Error desc: " << buffer << std::endl
            << "  At        : " << file << ':' << line << std::endl
            ;

  return -1; // must always return false or CHECK macro will break
}

#define CHECK_IGEOM( STR ) if (err != iBase_SUCCESS) return print_iGeom_error( STR, err, geom, __FILE__, __LINE__ )

#define CHECK_IMESH( STR ) if (err != iBase_SUCCESS) return print_iMesh_error( STR, err, instance, __FILE__, __LINE__ )

std::string default_file_name = "../../../MeshFiles/mesquite/3D/vtk/large_box_hex_1000.vtk";

void usage()
{
  cout << "main [-N] [filename]" << endl;
  cout << "  -N : Use native representation instead of TSTT implementation\n";
  cout << "  If no file name is specified, will use \""
       << default_file_name << '"' << endl;
  exit (1);
}

  // Construct a MeshTSTT from the file
int get_imesh_mesh( MBMesquite::Mesh** , const char* file_name, int dimension );

  // Construct a MeshImpl from the file
int get_native_mesh( MBMesquite::Mesh** , const char* file_name, int dimension );

int get_itaps_domain(MeshDomain**, const char*);
int get_mesquite_domain(MeshDomain**, const char*);

  // Run FeasibleNewton solver
int run_global_smoother( MeshDomainAssoc& mesh, MsqError& err );

  // Run SmoothLaplacian solver
int run_local_smoother( MeshDomainAssoc& mesh, MsqError& err );

int run_quality_optimizer( MeshDomainAssoc& mesh_and_domain, MsqError& err );

int main(int argc, char* argv[])
{
  MBMesquite::MsqPrintError err(cout);

    // command line arguments
  std::string file_name, geometry_file_name;
  bool use_native = false, opts_done = false;
  int dimension = 3;
  for (int arg = 1; arg < argc; ++arg)
  {
    if (!opts_done && argv[arg][0] == '-')
    {
      if (!strcmp( argv[arg], "-N"))
        use_native = true;
      else if(!strcmp( argv[arg], "-d2"))
        dimension = 2;
      else if(!strcmp( argv[arg], "--"))
        opts_done = true;
      else
        usage();
    }
    else if (!file_name.length()) {
        file_name = std::string(argv[arg]) + ".vtk";
        geometry_file_name = std::string(argv[arg]) + ".facet";
    }
    else
      usage();
  }
  if (!file_name.length())
  {
    file_name = default_file_name + ".vtk";
    geometry_file_name = default_file_name + ".facet";
    cout << "No file specified: using default: " << default_file_name << endl;
  }

  // Vector3D pnt(0,0,0);
  // Vector3D s_norm(0,0,1);
  // PlanarDomain plane(s_norm, pnt);

//  PlanarDomain plane( PlanarDomain::XY );
  MeshDomain* plane;
  int ierr;
  ierr = get_itaps_domain(&plane, geometry_file_name.c_str());
  if (ierr) return 1;

    // Try running a global smoother on the mesh
  Mesh* mesh=0;
  if(use_native)
    ierr = get_native_mesh(&mesh, file_name.c_str(), dimension);
  else
    ierr = get_imesh_mesh(&mesh, file_name.c_str(), dimension);

  if (!mesh || ierr) {
    std::cerr << "Failed to load input file.  Aborting." << std::endl;
    return 1;
  }

  MeshDomainAssoc mesh_and_domain(mesh, plane);

  MeshWriter::write_vtk(mesh, "original.vtk", err);
  if (err) return 1;
  cout << "Wrote \"original.vtk\"" << endl;

  run_global_smoother( mesh_and_domain, err );
  if (err) return 1;

    // Try running a local smoother on the mesh
//  Mesh* meshl=0;
//  if(use_native)
//    ierr = get_native_mesh(&meshl, file_name.c_str(), dimension);
//  else
//    ierr = get_imesh_mesh(&meshl, file_name.c_str(), dimension);
//  if (!mesh || ierr) {
//    std::cerr << "Failed to load input file.  Aborting." << std::endl;
//    return 1;
//  }

//  MeshDomainAssoc mesh_and_domain_local(meshl, plane);
  run_local_smoother( mesh_and_domain, err );
  if (err) return 1;

  run_quality_optimizer( mesh_and_domain, err );
  if (err) return 1;

  return 0;
}


int run_global_smoother( MeshDomainAssoc& mesh_and_domain, MsqError& err )
{
  double OF_value = 0.00001;

  Mesh* mesh = mesh_and_domain.get_mesh();

  // creates an intruction queue
  InstructionQueue queue1;

  // creates a mean ratio quality metric ...
  // IdealWeightInverseMeanRatio* mean_ratio = new IdealWeightInverseMeanRatio(err);
  // ConditionNumberQualityMetric* mean_ratio = new ConditionNumberQualityMetric();

  VertexConditionNumberQualityMetric* mean_ratio = new VertexConditionNumberQualityMetric();
  if (err) return 1;
  mean_ratio->set_averaging_method(QualityMetric::RMS, err);
  if (err) return 1;

  // ... and builds an objective function with it
  LPtoPTemplate* obj_func = new LPtoPTemplate(mean_ratio, 1, err);
  if (err) return 1;

  // creates the feas newt optimization procedures
  ConjugateGradient* pass1 = new ConjugateGradient( obj_func, err );
  // FeasibleNewton* pass1 = new FeasibleNewton( obj_func );
  pass1->use_global_patch();
  if (err) return 1;

  QualityAssessor stop_qa( mean_ratio );

  // **************Set stopping criterion****************
  TerminationCriterion tc_inner;
  tc_inner.add_absolute_vertex_movement( OF_value );
  if (err) return 1;
  TerminationCriterion tc_outer;
  tc_inner.add_iteration_limit( 15 );
  tc_outer.add_iteration_limit( 2 );
  pass1->set_inner_termination_criterion(&tc_inner);
  pass1->set_outer_termination_criterion(&tc_outer);

  queue1.add_quality_assessor(&stop_qa, err);
  if (err) return 1;

  // adds 1 pass of pass1 to mesh_set1
  queue1.set_master_quality_improver(pass1, err);
  if (err) return 1;

  queue1.add_quality_assessor(&stop_qa, err);
  if (err) return 1;

  // launches optimization on mesh_set
  queue1.run_instructions(&mesh_and_domain, err);
  if (err) return 1;

  MeshWriter::write_vtk(mesh, "feasible-newton-result.vtk", err);
  if (err) return 1;
  cout << "Wrote \"feasible-newton-result.vtk\"" << endl;

  //print_timing_diagnostics( cout );
  return 0;
}

int run_local_smoother( MeshDomainAssoc& mesh_and_domain, MsqError& err )
{
  double OF_value = 0.000001;

  Mesh* mesh = mesh_and_domain.get_mesh();

  // creates an intruction queue
  InstructionQueue queue1;

  // creates a mean ratio quality metric ...
  IdealWeightInverseMeanRatio* mean_ratio = new IdealWeightInverseMeanRatio(err);
  if (err) return 1;
  //  mean_ratio->set_gradient_type(QualityMetric::NUMERICAL_GRADIENT);
  //   mean_ratio->set_hessian_type(QualityMetric::NUMERICAL_HESSIAN);
  mean_ratio->set_averaging_method(QualityMetric::RMS, err);
  if (err) return 1;

  // ... and builds an objective function with it
  LPtoPTemplate* obj_func = new LPtoPTemplate(mean_ratio, 1, err);
  if (err) return 1;

  // creates the smart laplacian optimization procedures
  SmartLaplacianSmoother* pass1 = new SmartLaplacianSmoother( obj_func );

  QualityAssessor stop_qa( mean_ratio );

  // **************Set stopping criterion****************
  TerminationCriterion tc_inner;
  tc_inner.add_absolute_vertex_movement( OF_value );
  TerminationCriterion tc_outer;
  tc_outer.add_iteration_limit( 1 );
  pass1->set_inner_termination_criterion(&tc_inner);
  pass1->set_outer_termination_criterion(&tc_outer);

  queue1.add_quality_assessor(&stop_qa, err);
  if (err) return 1;

  // adds 1 pass of pass1 to mesh_set
  queue1.set_master_quality_improver(pass1, err);
  if (err) return 1;

  queue1.add_quality_assessor(&stop_qa, err);
  if (err) return 1;

  // launches optimization on mesh_set
  queue1.run_instructions(&mesh_and_domain, err);
  if (err) return 1;

  MeshWriter::write_vtk(mesh, "smart-laplacian-result.vtk", err);
  if (err) return 1;
  cout << "Wrote \"smart-laplacian-result.vtk\"" << endl;

  //print_timing_diagnostics( cout );
  return 0;
}


int run_quality_optimizer( MeshDomainAssoc& mesh_and_domain, MsqError& err )
{
  Mesh* mesh = mesh_and_domain.get_mesh();

  // creates an intruction queue
  InstructionQueue q;

  // do optimization
  const double eps = 0.01;
  IdealShapeTarget w;
  TShapeB1 mu;
  TQualityMetric metric( &w, &mu );
  PMeanPTemplate func( 1.0, &metric );

  SteepestDescent solver( &func );
  solver.use_element_on_vertex_patch();
  solver.do_jacobi_optimization();

  TerminationCriterion inner, outer;
  inner.add_absolute_vertex_movement( 0.5*eps );
  outer.add_absolute_vertex_movement( 0.5*eps );

  QualityAssessor qa( &metric );

  q.add_quality_assessor( &qa, err );
  if (err) return 1;
  q.set_master_quality_improver( &solver, err );
  if (err) return 1;
  q.add_quality_assessor( &qa, err );
  if (err) return 1;

  // launches optimization on mesh_set
  q.run_instructions( &mesh_and_domain, err );
  if (err) return 1;

  MeshWriter::write_vtk(mesh, "ideal-shape-result.vtk", err);
  if (err) return 1;
  cout << "Wrote \"ideal-shape-result.vtk\"" << endl;

  //print_timing_diagnostics( cout );
  return 0;
}


int get_imesh_mesh( MBMesquite::Mesh** mesh, const char* file_name, int dimension )
{
  int err;
  iMesh_Instance instance = 0;
  iMesh_newMesh( NULL, &instance, &err, 0 ); CHECK_IMESH("Creation of mesh instance failed");

  iBase_EntitySetHandle root_set;
  iMesh_getRootSet( instance, &root_set, &err );CHECK_IMESH("Could not get root set");

  iMesh_load( instance, root_set, file_name, 0, &err, strlen(file_name), 0 );CHECK_IMESH("Could not load mesh from file");

  iBase_TagHandle fixed_tag;
  iMesh_getTagHandle( instance, "fixed", &fixed_tag, &err, strlen("fixed") );
  if (iBase_SUCCESS != err) {
    // get the skin vertices of those cells and mark them as fixed; we don't want to fix the vertices on a
    // part boundary, but since we exchanged a layer of ghost cells, those vertices aren't on the skin locally
    // ok to mark non-owned skin vertices too, I won't move those anyway
    // use MOAB's skinner class to find the skin

    // get all vertices and cells
    // make tag to specify fixed vertices, since it's input to the algorithm; use a default value of non-fixed
    // so we only need to set the fixed tag for skin vertices
    moab::Interface* mbi = reinterpret_cast<MBiMesh*>(instance)->mbImpl;
    moab::EntityHandle currset=0;
    moab::Tag fixed;
    int def_val = 0;
    err=0;
    moab::ErrorCode rval = mbi->tag_get_handle("fixed", 1, moab::MB_TYPE_INTEGER, fixed, moab::MB_TAG_CREAT | moab::MB_TAG_DENSE, &def_val); CHECK_IMESH("Getting tag handle failed");
    moab::Range verts, cells, skin_verts;
    rval = mbi->get_entities_by_type(currset, moab::MBVERTEX, verts); CHECK_IMESH("Querying vertices failed");
    rval = mbi->get_entities_by_dimension(currset, dimension, cells); CHECK_IMESH("Querying elements failed");
    std::cout << "Found " << verts.size() << " vertices and " << cells.size() << " elements" << std::endl;

    moab::Skinner skinner(mbi);
    rval = skinner.find_skin(currset, cells, true, skin_verts); CHECK_IMESH("Finding the skin of the mesh failed"); // 'true' param indicates we want vertices back, not cells

    std::vector<int> fix_tag(skin_verts.size(), 1); // initialized to 1 to indicate fixed
    rval = mbi->tag_set_data(fixed, skin_verts, &fix_tag[0]); CHECK_IMESH("Setting tag data failed");;

    iMesh_getTagHandle( instance, "fixed", &fixed_tag, &err, strlen("fixed") );
  }


  MsqError ierr;
  MBMesquite::MsqIMesh* imesh = new MBMesquite::MsqIMesh( instance, root_set, (dimension == 3 ? iBase_REGION : iBase_FACE), ierr, &fixed_tag );
  if (MSQ_CHKERR(ierr)) {
    delete imesh;
    cerr << err << endl;
    err = iBase_FAILURE;
    CHECK_IMESH("Creation of MsqIMesh instance failed");
    return 0;
  }

  *mesh = imesh;
  return iBase_SUCCESS;
}



int get_native_mesh( MBMesquite::Mesh** mesh, const char* file_name, int  )
{
  MsqError err;
  MBMesquite::MeshImpl* imesh = new MBMesquite::MeshImpl();
  imesh->read_vtk( file_name, err );
  if (err)
  {
    cerr << err << endl;
    exit(3);
  }
  *mesh = imesh;

  return iBase_SUCCESS;
}


int get_itaps_domain(MeshDomain** odomain, const char* filename)
{
//  const double EPS = 1e-2;
  int err;
  iGeom_Instance geom;
  iGeom_newGeom( "", &geom, &err, 0 ); CHECK_IGEOM("ERROR: iGeom creation failed");

#ifdef MOAB_HAVE_CGM_FACET
  FacetModifyEngine::set_modify_enabled(CUBIT_TRUE);
#endif

  iGeom_load( geom, filename, 0, &err, strlen(filename), 0 );
  CHECK_IGEOM( "Cannot load the geometry" );

  iBase_EntitySetHandle root_set;
  iGeom_getRootSet( geom, &root_set, &err );
  CHECK_IGEOM( "getRootSet failed!" );

    // print out the number of entities
  std::cout << "Model contents: " << std::endl;
  const char *gtype[] = {"vertices: ", "edges: ", "faces: ", "regions: "};
  int nents[4];
  for (int i = 0; i <= 3; ++i) {
    iGeom_getNumOfType( geom, root_set, i, &nents[i], &err );
    CHECK_IGEOM( "Error: problem getting entities after gLoad." );
    std::cout << gtype[i] << nents[i] << std::endl;
  }

  iBase_EntityHandle* hd_geom_ents;
  int csize=0, sizealloc=0;
  if (nents[3] > 0) {
    hd_geom_ents = (iBase_EntityHandle*)malloc(sizeof(iBase_EntityHandle)*nents[2]);
    csize = nents[2];
    iGeom_getEntities(geom, root_set, 2, &hd_geom_ents, &csize, &sizealloc, &err);
  }
  else {
    hd_geom_ents = (iBase_EntityHandle*)malloc(sizeof(iBase_EntityHandle)*nents[1]);
    csize = nents[1];
    iGeom_getEntities(geom, root_set, 1, &hd_geom_ents, &csize, &sizealloc, &err);
  }
  CHECK_IGEOM( "ERROR: Could not get entities" );

  *odomain = new MsqIGeom( geom, hd_geom_ents[0] );
  return iBase_SUCCESS;
}

