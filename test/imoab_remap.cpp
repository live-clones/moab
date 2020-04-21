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

#define CHECKIERR(ierr, message)  if (0!=ierr) { printf ("%s", message); return 1;}

// this test will be run in serial only
int main(int argc, char * argv[])
{
  ErrCode ierr;
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
  ierr = iMOAB_Initialize(argc, argv);
  CHECKIERR(ierr, "failed to initialize MOAB");

  int atmAppID, ocnAppID, lndAppID, atmocnAppID, atmlndAppID, lndatmAppID;
  iMOAB_AppID atmPID=&atmAppID;
  iMOAB_AppID ocnPID=&ocnAppID;
  iMOAB_AppID lndPID=&lndAppID;
  iMOAB_AppID atmocnPID=&atmocnAppID;
  iMOAB_AppID atmlndPID=&atmlndAppID;
  iMOAB_AppID lndatmPID=&lndatmAppID;
  /*
   * Each application has to be registered once. A mesh set and a parallel communicator will be associated
   * with each application. A unique application id will be returned, and will be used for all future
   * mesh operations/queries.
   */
  int atmCompID = 10, ocnCompID = 20, lndCompID = 30, atmocnCompID = 100, atmlndCompID = 102, lndatmCompID = 103;
  ierr = iMOAB_RegisterApplication( "ATM-APP",
#ifdef MOAB_HAVE_MPI
      &comm,
#endif
      &atmCompID,
      atmPID);
  CHECKIERR(ierr, "failed to register application1");

  ierr = iMOAB_RegisterApplication( "OCN-APP",
#ifdef MOAB_HAVE_MPI
      &comm,
#endif
      &ocnCompID,
      ocnPID);
  CHECKIERR(ierr, "failed to register application2");

  ierr = iMOAB_RegisterApplication( "LND-APP",
#ifdef MOAB_HAVE_MPI
      &comm,
#endif
      &lndCompID,
      lndPID);
  CHECKIERR(ierr, "failed to register application3");

  ierr = iMOAB_RegisterApplication( "ATM-OCN-CPL",
#ifdef MOAB_HAVE_MPI
      &comm,
#endif
      &atmocnCompID,
      atmocnPID);
  CHECKIERR(ierr, "failed to register application4");

  ierr = iMOAB_RegisterApplication( "ATM-LND-CPL",
#ifdef MOAB_HAVE_MPI
      &comm,
#endif
      &atmlndCompID,
      atmlndPID);
  CHECKIERR(ierr, "failed to register application5");

  ierr = iMOAB_RegisterApplication( "LND-ATM-CPL",
#ifdef MOAB_HAVE_MPI
      &comm,
#endif
      &lndatmCompID,
      lndatmPID);
  CHECKIERR(ierr, "failed to register application5");

  const char *read_opts="";
  int num_ghost_layers=0;
  /*
   * Loading the mesh is a parallel IO operation. Ghost layers can be exchanged too, and default MOAB
   * sets are augmented with ghost elements. By convention, blocks correspond to MATERIAL_SET sets,
   * side sets with NEUMANN_SET sets, node sets with DIRICHLET_SET sets. Each element and vertex entity should have
   * a GLOBAL ID tag in the file, which will be available for visible entities
   */
  ierr = iMOAB_LoadMesh(  atmPID, atmFilename.c_str(), read_opts, &num_ghost_layers, atmFilename.size(), strlen(read_opts) );
  CHECKIERR(ierr, "failed to load mesh");
  ierr = iMOAB_LoadMesh(  ocnPID, ocnFilename.c_str(), read_opts, &num_ghost_layers, ocnFilename.size(), strlen(read_opts) );
  CHECKIERR(ierr, "failed to load mesh");
  ierr = iMOAB_LoadMesh(  lndPID, lndFilename.c_str(), read_opts, &num_ghost_layers, lndFilename.size(), strlen(read_opts) );
  CHECKIERR(ierr, "failed to load mesh");

  int nverts[3], nelem[3];
  /*
   * Each process in the communicator will have access to a local mesh instance, which will contain the
   * original cells in the local partition and ghost entities. Number of vertices, primary cells, visible blocks,
   * number of sidesets and nodesets boundary conditions will be returned in size 3 arrays, for local, ghost and total
   * numbers.
   */
  ierr = iMOAB_GetMeshInfo( atmPID, nverts, nelem, 0, 0, 0);
  CHECKIERR(ierr, "failed to get mesh info");
  printf("Atmosphere Component Mesh: %d vertices and %d elements\n", nverts[0], nelem[0]);

  ierr = iMOAB_GetMeshInfo( ocnPID, nverts, nelem, 0, 0, 0);
  CHECKIERR(ierr, "failed to get mesh info");
  printf("Ocean Component Mesh: %d vertices and %d elements\n", nverts[0], nelem[0]);

  ierr = iMOAB_GetMeshInfo( lndPID, nverts, nelem, 0, 0, 0);
  CHECKIERR(ierr, "failed to get mesh info");
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
  const char* bottomTempFieldATM = "a2oTbotATM";
  const char* bottomTempProjectedField = "a2oTbot_proj";
  const char* bottomTempProjectedNCField = "a2oTbot_projnocons";
  int tagIndex[4];
  // int entTypes[2] = {1, 1}; /* both on elements; */
  int tagTypes[2] = { DENSE_DOUBLE, DENSE_DOUBLE } ;
  int atmCompNDoFs = disc_orders[0]*disc_orders[0], ocnCompNDoFs = 1 /*FV*/;

  ierr = iMOAB_DefineTagStorage(atmPID, bottomTempField, &tagTypes[0], &atmCompNDoFs, &tagIndex[0],  strlen(bottomTempField) );
  CHECKIERR(ierr, "failed to define the field tag on ATM");

  ierr = iMOAB_DefineTagStorage(atmPID, bottomTempFieldATM, &tagTypes[0], &atmCompNDoFs, &tagIndex[0],  strlen(bottomTempFieldATM) );
  CHECKIERR(ierr, "failed to define the field tag on ATM");

  ierr = iMOAB_DefineTagStorage(ocnPID, bottomTempProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndex[1],  strlen(bottomTempProjectedField) );
  CHECKIERR(ierr, "failed to define the field tag on OCN");

  ierr = iMOAB_DefineTagStorage(ocnPID, bottomTempProjectedNCField, &tagTypes[1], &ocnCompNDoFs, &tagIndex[2],  strlen(bottomTempProjectedNCField) );
  CHECKIERR(ierr, "failed to define the field tag on OCN");

  ierr = iMOAB_DefineTagStorage(lndPID, bottomTempProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndex[3],  strlen(bottomTempProjectedField) );
  CHECKIERR(ierr, "failed to define the field tag on LND");

  /* Next compute the mesh intersection on the sphere between the source (ATM) and target (OCN) meshes */
  ierr = iMOAB_ComputeMeshIntersectionOnSphere(atmPID, ocnPID, atmocnPID);
  CHECKIERR(ierr, "failed to compute mesh intersection between ATM and OCN");

  /* Next compute the mesh intersection on the sphere between the source (ATM) and target (LND) meshes */
  ierr = iMOAB_ComputePointDoFIntersection(atmPID, lndPID, atmlndPID);
  CHECKIERR(ierr, "failed to compute point-cloud mapping ATM-LND");

  /* Next compute the mesh intersection on the sphere between the source (LND) and target (ATM) meshes */
  ierr = iMOAB_ComputePointDoFIntersection(lndPID, atmPID, lndatmPID);
  CHECKIERR(ierr, "failed to compute point-cloud mapping LND-ATM");

  /* We have the mesh intersection now. Let us compute the remapping weights */
  fNoConserve=1;
  /* Compute the weights to preoject the solution from ATM component to OCN compoenent */
  ierr = iMOAB_ComputeScalarProjectionWeights ( atmocnPID,
                                              weights_identifiers[0],
                                              disc_methods[0], &disc_orders[0],
                                              disc_methods[1], &disc_orders[1],
                                              &fMonotoneTypeID, &fVolumetric, &fNoConserve, &fValidate,
                                              dof_tag_names[0], dof_tag_names[1],
                                              strlen(weights_identifiers[0]),
                                              strlen(disc_methods[0]), strlen(disc_methods[1]),
                                              strlen(dof_tag_names[0]), strlen(dof_tag_names[1])
                                            );
  CHECKIERR(ierr, "failed to compute remapping projection weights for ATM-OCN scalar non-conservative field");

  /* Compute the weights to preoject the solution from ATM component to LND compoenent */
  ierr = iMOAB_ComputeScalarProjectionWeights ( atmlndPID,
                                              weights_identifiers[1],
                                              disc_methods[0], &disc_orders[0],
                                              disc_methods[2], &disc_orders[2],
                                              &fMonotoneTypeID, &fVolumetric, &fNoConserve, &fValidate,
                                              dof_tag_names[0], dof_tag_names[2],
                                              strlen(weights_identifiers[1]),
                                              strlen(disc_methods[0]), strlen(disc_methods[2]),
                                              strlen(dof_tag_names[0]), strlen(dof_tag_names[2])
                                            );
  CHECKIERR(ierr, "failed to compute remapping projection weights for ATM-LND scalar non-conservative field");

  /* Compute the weights to preoject the solution from ATM component to LND compoenent */
  ierr = iMOAB_ComputeScalarProjectionWeights ( lndatmPID,
                                              weights_identifiers[1],
                                              disc_methods[2], &disc_orders[2],
                                              disc_methods[0], &disc_orders[0],
                                              &fMonotoneTypeID, &fVolumetric, &fNoConserve, &fValidate,
                                              dof_tag_names[2], dof_tag_names[0],
                                              strlen(weights_identifiers[1]),
                                              strlen(disc_methods[2]), strlen(disc_methods[0]),
                                              strlen(dof_tag_names[2]), strlen(dof_tag_names[0])
                                            );
  CHECKIERR(ierr, "failed to compute remapping projection weights for LND-ATM scalar non-conservative field");

  /* We have the mesh intersection now. Let us compute the remapping weights */
  fNoConserve=0;
  ierr = iMOAB_ComputeScalarProjectionWeights ( atmocnPID,
                                              weights_identifiers[2],
                                              disc_methods[0], &disc_orders[0],
                                              disc_methods[1], &disc_orders[1],
                                              &fMonotoneTypeID, &fVolumetric, &fNoConserve, &fValidate,
                                              dof_tag_names[0], dof_tag_names[1],
                                              strlen(weights_identifiers[1]),
                                              strlen(disc_methods[0]), strlen(disc_methods[1]),
                                              strlen(dof_tag_names[0]), strlen(dof_tag_names[1])
                                            );
  CHECKIERR(ierr, "failed to compute remapping projection weights for scalar conservative field");

  /* We have the remapping weights now. Let us apply the weights onto the tag we defined
     on the srouce mesh and get the projection on the target mesh */
  ierr = iMOAB_ApplyScalarProjectionWeights ( atmocnPID,
                                            weights_identifiers[0],
                                            bottomTempField,
                                            bottomTempProjectedNCField,
                                            strlen(weights_identifiers[0]),
                                            strlen(bottomTempField),
                                            strlen(bottomTempProjectedNCField)
                                            );
  CHECKIERR(ierr, "failed to apply projection weights for scalar non-conservative field");

  /* We have the remapping weights now. Let us apply the weights onto the tag we defined
     on the srouce mesh and get the projection on the target mesh */
  ierr = iMOAB_ApplyScalarProjectionWeights ( atmocnPID,
                                            weights_identifiers[2],
                                            bottomTempField,
                                            bottomTempProjectedField,
                                            strlen(weights_identifiers[1]),
                                            strlen(bottomTempField),
                                            strlen(bottomTempProjectedField)
                                            );
  CHECKIERR(ierr, "failed to apply projection weights for scalar conservative field");

  /* We have the remapping weights now. Let us apply the weights onto the tag we defined
     on the srouce mesh and get the projection on the target mesh */
  ierr = iMOAB_ApplyScalarProjectionWeights ( atmlndPID,
                                            weights_identifiers[1],
                                            bottomTempField,
                                            bottomTempProjectedField,
                                            strlen(weights_identifiers[1]),
                                            strlen(bottomTempField),
                                            strlen(bottomTempProjectedField)
                                            );
  CHECKIERR(ierr, "failed to apply projection weights for ATM-LND scalar field");

  /* We have the remapping weights now. Let us apply the weights onto the tag we defined
     on the srouce mesh and get the projection on the target mesh */
  ierr = iMOAB_ApplyScalarProjectionWeights ( lndatmPID,
                                            weights_identifiers[1],
                                            bottomTempProjectedField,
                                            bottomTempFieldATM,
                                            strlen(weights_identifiers[1]),
                                            strlen(bottomTempField),
                                            strlen(bottomTempProjectedField)
                                            );
  CHECKIERR(ierr, "failed to apply projection weights for LND-ATM scalar field");

  /*
   * the file can be written in parallel, and it will contain additional tags defined by the user
   * we may extend the method to write only desired tags to the file
   */
  {
    // free allocated data
    char outputFileAtmTgt[] = "fIntxAtmTarget.h5m";
    char outputFileOcnTgt[] = "fIntxOcnTarget.h5m";
    char outputFileLndTgt[]  = "fIntxLndTarget.h5m";
    char writeOptions[] ="";

    ierr = iMOAB_WriteMesh(ocnPID, outputFileOcnTgt, writeOptions,
      strlen(outputFileOcnTgt), strlen(writeOptions) );

    ierr = iMOAB_WriteMesh(lndPID, outputFileLndTgt, writeOptions,
      strlen(outputFileLndTgt), strlen(writeOptions) );

    ierr = iMOAB_WriteMesh(atmPID, outputFileAtmTgt, writeOptions,
      strlen(outputFileAtmTgt), strlen(writeOptions) );
  }

  /*
   * deregistering application will delete all mesh entities associated with the application and will
   *  free allocated tag storage.
   */
  ierr = iMOAB_DeregisterApplication(lndPID);
  CHECKIERR(ierr, "failed to de-register application3");
  ierr = iMOAB_DeregisterApplication(ocnPID);
  CHECKIERR(ierr, "failed to de-register application2");
  ierr = iMOAB_DeregisterApplication(atmPID);
  CHECKIERR(ierr, "failed to de-register application1");

  /*
   * this method will delete MOAB instance
   */
  ierr = iMOAB_Finalize();
  CHECKIERR(ierr, "failed to finalize MOAB");

#ifdef MOAB_HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}

