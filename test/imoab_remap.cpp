#include "moab/MOABConfig.h"

#ifdef MOAB_HAVE_MPI
#  include "moab_mpi.h"
#endif
#include "moab/iMOAB.h"
#include "TestUtil.hpp"
#include "moab/CpuTimer.hpp"
#include "moab/ProgOptions.hpp"

// for malloc, free:
#include <iostream>

#define CHECKRC(rc, message)  if (0!=rc) { printf ("%s", message); return 1;}

// this test will be run in serial only
int main(int argc, char * argv[])
{
  std::string atmFilename, ocnFilename, lndFilename;

  atmFilename = TestDir + "/wholeATM_T.h5m";
  ocnFilename = TestDir + "/recMeshOcn.h5m";
  lndFilename = TestDir + "/wholeLnd.h5m";

  ProgOptions opts;
  opts.addOpt<std::string>("atmosphere,t", "atm mesh filename (source)", &atmFilename);
  opts.addOpt<std::string>("ocean,m", "ocean mesh filename (target)", &ocnFilename);
  opts.addOpt<std::string>("land,l", "land mesh filename (target)", &lndFilename);

  opts.parseCommandLine(argc, argv);

  {
    std::cout << " atm file: " << atmFilename << 
        "\n ocn file: " << ocnFilename << 
        "\n lnd file: " << lndFilename << "\n";
  }

#ifdef MOAB_HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
#endif
  /*
   * MOAB needs to be initialized; A MOAB instance will be created, and will be used by each application
   * in this framework. There is no parallel context yet.
   */
  ErrCode rc = iMOAB_Initialize(argc, argv);
  CHECKRC(rc, "failed to initialize MOAB");

  int atmAppID, ocnAppID, lndAppID, atmocnAppID, atmlndAppID;
  iMOAB_AppID atmPID=&atmAppID;
  iMOAB_AppID ocnPID=&ocnAppID;
  iMOAB_AppID lndPID=&lndAppID;
  iMOAB_AppID atmocnPID=&atmocnAppID;
  iMOAB_AppID atmlndPID=&atmlndAppID;
  /*
   * Each application has to be registered once. A mesh set and a parallel communicator will be associated
   * with each application. A unique application id will be returned, and will be used for all future
   * mesh operations/queries.
   */
  int atmCompID = 10, ocnCompID = 20, lndCompID = 30, atmocnCompID = 100, atmlndCompID = 101;
  rc = iMOAB_RegisterApplication( "ATM-APP",
#ifdef MOAB_HAVE_MPI
      &comm,
#endif
      &atmCompID,
      atmPID);
  CHECKRC(rc, "failed to register application1");

  rc = iMOAB_RegisterApplication( "OCN-APP",
#ifdef MOAB_HAVE_MPI
      &comm,
#endif
      &ocnCompID,
      ocnPID);
  CHECKRC(rc, "failed to register application2");

  rc = iMOAB_RegisterApplication( "LND-APP",
#ifdef MOAB_HAVE_MPI
      &comm,
#endif
      &lndCompID,
      lndPID);
  CHECKRC(rc, "failed to register application3");

  rc = iMOAB_RegisterApplication( "ATM-OCN-CPL",
#ifdef MOAB_HAVE_MPI
      &comm,
#endif
      &atmocnCompID,
      atmocnPID);
  CHECKRC(rc, "failed to register application4");

  rc = iMOAB_RegisterApplication( "ATM-LND-CPL",
#ifdef MOAB_HAVE_MPI
      &comm,
#endif
      &atmlndCompID,
      atmlndPID);
  CHECKRC(rc, "failed to register application5");

  const char *read_opts="";
  int num_ghost_layers=0;

  /*
   * Loading the mesh is a parallel IO operation. Ghost layers can be exchanged too, and default MOAB
   * sets are augmented with ghost elements. By convention, blocks correspond to MATERIAL_SET sets,
   * side sets with NEUMANN_SET sets, node sets with DIRICHLET_SET sets. Each element and vertex entity should have
   * a GLOBAL ID tag in the file, which will be available for visible entities
   */
  rc = iMOAB_LoadMesh(  atmPID, atmFilename.c_str(), read_opts, &num_ghost_layers, atmFilename.size(), strlen(read_opts) );
  CHECKRC(rc, "failed to load mesh");
  rc = iMOAB_LoadMesh(  ocnPID, ocnFilename.c_str(), read_opts, &num_ghost_layers, ocnFilename.size(), strlen(read_opts) );
  CHECKRC(rc, "failed to load mesh");
  rc = iMOAB_LoadMesh(  lndPID, lndFilename.c_str(), read_opts, &num_ghost_layers, lndFilename.size(), strlen(read_opts) );
  CHECKRC(rc, "failed to load mesh");

  int nverts[3], nelem[3];
  /*
   * Each process in the communicator will have access to a local mesh instance, which will contain the
   * original cells in the local partition and ghost entities. Number of vertices, primary cells, visible blocks,
   * number of sidesets and nodesets boundary conditions will be returned in size 3 arrays, for local, ghost and total
   * numbers.
   */
  rc = iMOAB_GetMeshInfo( atmPID, nverts, nelem, 0, 0, 0);
  CHECKRC(rc, "failed to get mesh info");
  printf("Atmosphere Component Mesh: %d vertices and %d elements\n", nverts[0], nelem[0]);

  rc = iMOAB_GetMeshInfo( ocnPID, nverts, nelem, 0, 0, 0);
  CHECKRC(rc, "failed to get mesh info");
  printf("Ocean Component Mesh: %d vertices and %d elements\n", nverts[0], nelem[0]);

  rc = iMOAB_GetMeshInfo( lndPID, nverts, nelem, 0, 0, 0);
  CHECKRC(rc, "failed to get mesh info");
  printf("Land Component Mesh: %d vertices and %d elements\n", nverts[0], nelem[0]);

  /*
   * The 2 tags used in this example exist in the file, already.
   * If a user needs a new tag, it can be defined with the same call as this one
   * this method, iMOAB_DefineTagStorage, will return a local index for the tag.
   * The name of the tag is case sensitive.
   * This method is collective.
   */
  // int disc_orders[2] = {1, 1};
  // const char* disc_methods[2] = {"fv", "fv"};
  // const char* dof_tag_names[2] = {"GLOBAL_ID", "GLOBAL_ID"};
  int disc_orders[3] = {4, 1, 1};
  const char* disc_methods[3] = {"cgll", "fv", "pcloud"};
  const char* dof_tag_names[3] = {"GLOBAL_DOFS", "GLOBAL_ID", "GLOBAL_ID"};
  
  const char* weights_identifiers[3] = {"scalar", "scalar_pointcloud", "scalar_conservative"};
  int fMonotoneTypeID=0, fVolumetric=0, fValidate=1, fNoConserve=0;

  const char* bottomTempField = "a2oTbot";
  const char* bottomTempProjectedField = "a2oTbot_proj";
  const char* bottomTempProjectedNCField = "a2oTbot_projnocons";
  int tagIndex[4];
  // int entTypes[2] = {1, 1}; /* both on elements; */
  int tagTypes[2] = { DENSE_DOUBLE, DENSE_DOUBLE } ;
  int atmCompNDoFs = disc_orders[0]*disc_orders[0], ocnCompNDoFs = disc_orders[1]*disc_orders[1];

  rc = iMOAB_DefineTagStorage(atmPID, bottomTempField, &tagTypes[0], &atmCompNDoFs, &tagIndex[0],  strlen(bottomTempField) );
  CHECKRC(rc, "failed to define the field tag");

  rc = iMOAB_DefineTagStorage(ocnPID, bottomTempProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndex[1],  strlen(bottomTempProjectedField) );
  CHECKRC(rc, "failed to define the field tag");

  rc = iMOAB_DefineTagStorage(ocnPID, bottomTempProjectedNCField, &tagTypes[1], &ocnCompNDoFs, &tagIndex[2],  strlen(bottomTempProjectedNCField) );
  CHECKRC(rc, "failed to define the field tag");

  rc = iMOAB_DefineTagStorage(lndPID, bottomTempProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndex[3],  strlen(bottomTempProjectedField) );
  CHECKRC(rc, "failed to define the field tag");

  /* Next compute the mesh intersection on the sphere between the source and target meshes */
  rc = iMOAB_ComputeMeshIntersectionOnSphere(atmPID, ocnPID, atmocnPID);
  CHECKRC(rc, "failed to compute mesh intersection");

  /* Next compute the mesh intersection on the sphere between the source and target meshes */
  // rc = iMOAB_ComputeMeshIntersectionOnSphere(atmPID, lndPID, atmlndPID);
  // CHECKRC(rc, "failed to compute mesh intersection");
  // No need to comoute intersection, but still need to compute coverage nad translate covering set to Tempest datastructure
  // m_covering_source = new Mesh();
  // rval = convert_mesh_to_tempest_private ( m_covering_source, m_covering_source_set, m_covering_source_entities, &m_covering_source_vertices ); MB_CHK_SET_ERR ( rval, "Can't convert source Tempest mesh" );
  

  rc = iMOAB_ComputePointDoFIntersection(atmPID, lndPID, atmlndPID, 
                                          disc_methods[0], &disc_orders[0], dof_tag_names[0], 
                                          disc_methods[2], &disc_orders[2], dof_tag_names[2],
                                          strlen(disc_methods[0]), strlen(dof_tag_names[0]), 
                                          strlen(disc_methods[1]), strlen(dof_tag_names[1]));
  CHECKRC(rc, "failed to compute mesh intersection");

  /* We have the mesh intersection now. Let us compute the remapping weights */
  fNoConserve=1;
  /* Compute the weights to preoject the solution from ATM component to OCN compoenent */
  rc = iMOAB_ComputeScalarProjectionWeights ( atmocnPID,
                                              weights_identifiers[0],
                                              disc_methods[0], &disc_orders[0],
                                              disc_methods[1], &disc_orders[1],
                                              &fMonotoneTypeID, &fVolumetric, &fNoConserve, &fValidate,
                                              dof_tag_names[0], dof_tag_names[1],
                                              strlen(weights_identifiers[0]),
                                              strlen(disc_methods[0]), strlen(disc_methods[1]),
                                              strlen(dof_tag_names[0]), strlen(dof_tag_names[1])
                                            );
  CHECKRC(rc, "failed to compute remapping projection weights for ATM-OCN scalar non-conservative field");

  /* Compute the weights to preoject the solution from ATM component to LND compoenent */
  rc = iMOAB_ComputeScalarProjectionWeights ( atmlndPID,
                                              weights_identifiers[1],
                                              disc_methods[0], &disc_orders[0],
                                              disc_methods[2], &disc_orders[2],
                                              &fMonotoneTypeID, &fVolumetric, &fNoConserve, &fValidate,
                                              dof_tag_names[0], dof_tag_names[2],
                                              strlen(weights_identifiers[0]),
                                              strlen(disc_methods[0]), strlen(disc_methods[2]),
                                              strlen(dof_tag_names[0]), strlen(dof_tag_names[2])
                                            );
  CHECKRC(rc, "failed to compute remapping projection weights for ATM-LND scalar non-conservative field");

  /* We have the mesh intersection now. Let us compute the remapping weights */
  fNoConserve=0;
  rc = iMOAB_ComputeScalarProjectionWeights ( atmocnPID,
                                              weights_identifiers[2],
                                              disc_methods[0], &disc_orders[0],
                                              disc_methods[1], &disc_orders[1],
                                              &fMonotoneTypeID, &fVolumetric, &fNoConserve, &fValidate,
                                              dof_tag_names[0], dof_tag_names[1],
                                              strlen(weights_identifiers[1]),
                                              strlen(disc_methods[0]), strlen(disc_methods[1]),
                                              strlen(dof_tag_names[0]), strlen(dof_tag_names[1])
                                            );
  CHECKRC(rc, "failed to compute remapping projection weights for scalar conservative field");

  /* We have the remapping weights now. Let us apply the weights onto the tag we defined
     on the srouce mesh and get the projection on the target mesh */
  rc = iMOAB_ApplyScalarProjectionWeights ( atmocnPID,
                                            weights_identifiers[0],
                                            bottomTempField,
                                            bottomTempProjectedNCField,
                                            strlen(weights_identifiers[0]),
                                            strlen(bottomTempField),
                                            strlen(bottomTempProjectedNCField)
                                            );
  CHECKRC(rc, "failed to compute projection weight application for scalar non-conservative field");

  /* We have the remapping weights now. Let us apply the weights onto the tag we defined
     on the srouce mesh and get the projection on the target mesh */
  rc = iMOAB_ApplyScalarProjectionWeights ( atmocnPID,
                                            weights_identifiers[2],
                                            bottomTempField,
                                            bottomTempProjectedField,
                                            strlen(weights_identifiers[1]),
                                            strlen(bottomTempField),
                                            strlen(bottomTempProjectedField)
                                            );
  CHECKRC(rc, "failed to compute projection weight application for scalar conservative field");

  /* We have the remapping weights now. Let us apply the weights onto the tag we defined
     on the srouce mesh and get the projection on the target mesh */
  rc = iMOAB_ApplyScalarProjectionWeights ( atmlndPID,
                                            weights_identifiers[1],
                                            bottomTempField,
                                            bottomTempProjectedField,
                                            strlen(weights_identifiers[1]),
                                            strlen(bottomTempField),
                                            strlen(bottomTempProjectedField)
                                            );
  CHECKRC(rc, "failed to compute projection weight application for scalar conservative field");

  /*
   * the file can be written in parallel, and it will contain additional tags defined by the user
   * we may extend the method to write only desired tags to the file
   */
  {
    // free allocated data
    char outputFileOv[]  = "fIntxOverlap.h5m";
    char outputFileTgt[] = "fIntxTarget.h5m";
    char writeOptions[] ="";

    rc = iMOAB_WriteMesh(ocnPID, outputFileTgt, writeOptions,
      strlen(outputFileTgt), strlen(writeOptions) );

    rc = iMOAB_WriteMesh(lndPID, outputFileOv, writeOptions,
      strlen(outputFileOv), strlen(writeOptions) );
  }

  /*
   * deregistering application will delete all mesh entities associated with the application and will
   *  free allocated tag storage.
   */
  rc = iMOAB_DeregisterApplication(lndPID);
  CHECKRC(rc, "failed to de-register application3");
  rc = iMOAB_DeregisterApplication(ocnPID);
  CHECKRC(rc, "failed to de-register application2");
  rc = iMOAB_DeregisterApplication(atmPID);
  CHECKRC(rc, "failed to de-register application1");

  /*
   * this method will delete MOAB instance
   */
  rc = iMOAB_Finalize();
  CHECKRC(rc, "failed to finalize MOAB");

#ifdef MOAB_HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}

