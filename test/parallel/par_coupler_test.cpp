/*
 * par_coupler_test.cpp
 * this will exercise the coupler in parallel.
 * It is mainly the same test as mbcoupler_test default one, but it is launched by default on 2 processors,
 *  as part of make check
 *  Created on: Nov 14, 2012
 */

#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#include "moab/Core.hpp"
#include "Coupler.hpp"
#include "iMesh_extensions.h"
#include "moab_mpi.h"
#include "ElemUtil.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <assert.h>
#include "MBiMesh.hpp"

using namespace moab;

#define STRINGIFY_(A) #A
#define STRINGIFY(A) STRINGIFY_(A)
#ifdef MESHDIR
std::string TestDir( STRINGIFY(MESHDIR) );
#else
std::string TestDir(".");
#endif

#define RRA(a) if (MB_SUCCESS != result) {\
      std::string tmp_str; mbImpl->get_last_error(tmp_str);\
      tmp_str.append("\n"); tmp_str.append(a);\
      dynamic_cast<Core*>(mbImpl)->get_error_handler()->set_last_error(tmp_str.c_str()); \
      return result;}

#define PRINT_LAST_ERROR \
    if (MB_SUCCESS != result) {\
      std::string tmp_str;\
      std::cout << "Failure; message:" << std::endl;\
      mbImpl->get_last_error(tmp_str);\
      std::cout << tmp_str << std::endl;\
      MPI_Abort(MPI_COMM_WORLD, result);        \
      return result;\
    }

#define PRINT_LAST_ERR \
    if (iBase_SUCCESS != err) {\
      std::string tmp_str;\
      std::cout << "Failure; message:" << std::endl;\
      mbImpl->get_last_error(tmp_str);\
      std::cout << tmp_str << std::endl;\
      MPI_Abort(MPI_COMM_WORLD, result);        \
      return result;\
    }

void print_usage();

ErrorCode get_file_options(int argc, char **argv,
                           std::vector<std::string> &meshFiles,
                           Coupler::Method &method,
                           std::string &interpTag,
                           std::string &gNormTag,
                           std::string &ssNormTag,
                           std::vector<const char *> &ssTagNames,
                           std::vector<const char *> &ssTagValues,
                           std::string &readOpts,
                           std::string &outFile,
                           std::string &writeOpts,
                           std::string &dbgFile,
                           bool &help,
                           double & epsilon);

// ErrorCode get_file_options(int argc, char **argv,
//                              std::vector<const char *> &filenames,
//                              std::string &tag_name,
//                              std::string &out_fname,
//                              std::string &opts);

ErrorCode report_iface_ents(Interface *mbImpl,
                              std::vector<ParallelComm *> &pcs,
                              bool print_results);

ErrorCode test_interpolation(Interface *mbImpl,
                             Coupler::Method method,
                             std::string &interpTag,
                             std::string &gNormTag,
                             std::string &ssNormTag,
                             std::vector<const char *> &ssTagNames,
                             std::vector<const char *> &ssTagValues,
                             iBase_EntitySetHandle *roots,
                             std::vector<ParallelComm *> &pcs,
                             double &instant_time,
                             double &pointloc_time,
                             double &interp_time,
                             double &gnorm_time,
                             double &ssnorm_time,
                             double & toler);

void reduceMax(double &v){
  double buf;

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce( &v, &buf, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  v=buf;
}


int main(int argc, char **argv)
{
    // need to init MPI first, to tell how many procs and rank
  int err = MPI_Init(&argc, &argv);

  std::vector<const char *> ssTagNames, ssTagValues;
  std::vector<std::string> meshFiles;
  std::string interpTag, gNormTag, ssNormTag, readOpts, outFile, writeOpts, dbgFile;
  Coupler::Method method = Coupler::CONSTANT;

  ErrorCode result = MB_SUCCESS;
  bool help = false;
  double toler = 5.e-10;
  result = get_file_options(argc, argv, meshFiles, method, interpTag,
                            gNormTag, ssNormTag, ssTagNames, ssTagValues,
                            readOpts, outFile, writeOpts, dbgFile, help, toler);

  if (result != MB_SUCCESS || help) {
    print_usage();
    err = MPI_Finalize();
    return 1;
  }

  int nprocs, rank;
  err = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // redirect stdout and stderr if dbgFile is not null
  if (!dbgFile.empty()) {
    std::stringstream dfname;
    dfname << dbgFile << rank << ".txt";
    std::freopen(dfname.str().c_str(), "a", stdout);
    std::freopen(dfname.str().c_str(), "a", stderr);
  }

  // create MOAB instance based on that
  Interface *mbImpl = new Core();
  if (NULL == mbImpl) return 1;

  iMesh_Instance iMeshInst = reinterpret_cast<iMesh_Instance>( new MBiMesh(mbImpl) );

    // read in mesh(es)
  std::vector<ParallelComm *> pcs(meshFiles.size());

    // Create root sets for each mesh using the iMesh API.  Then pass these
    // to the load_file functions to be populated.
  iBase_EntitySetHandle *roots = (iBase_EntitySetHandle *) malloc(sizeof(iBase_EntitySetHandle) * meshFiles.size());

  for (unsigned int i = 0; i < meshFiles.size(); i++) {
    pcs[i] = new ParallelComm(mbImpl, MPI_COMM_WORLD);
    int index = pcs[i]->get_id();
    std::string newReadopts;
    std::ostringstream extraOpt;
    extraOpt  << ";PARALLEL_COMM=" << index;
    newReadopts = readOpts+extraOpt.str();

    iMesh_createEntSet(iMeshInst, 0, &(roots[i]), &err);
    result = mbImpl->load_file( meshFiles[i].c_str(), (EntityHandle *)&roots[i], newReadopts.c_str() );
    //result = rps[i]->load_file(meshFiles[i].c_str(), (EntityHandle *)&roots[i], FileOptions(readOpts.c_str()));
    PRINT_LAST_ERROR;
  }

  result = report_iface_ents(mbImpl, pcs, true);
  PRINT_LAST_ERROR;

  double instant_time=0.0, pointloc_time=0.0, interp_time=0.0, gnorm_time=0.0, ssnorm_time=0.0;
    // test interpolation and global normalization and subset normalization

  result = test_interpolation(mbImpl, method, interpTag, gNormTag, ssNormTag,
                              ssTagNames, ssTagValues, roots, pcs,
                              instant_time, pointloc_time, interp_time,
                              gnorm_time, ssnorm_time, toler);
  PRINT_LAST_ERROR;


  reduceMax(instant_time);
  reduceMax(pointloc_time);
  reduceMax(interp_time);

  if(0==rank)
    printf("\nMax time : %g %g %g (inst loc interp -- %d procs )\n", instant_time,
     pointloc_time, interp_time, nprocs);

    // output mesh
  if (!outFile.empty()) {
    Range partSets;
      // only save the target mesh
    partSets.insert((EntityHandle)roots[1]);
    std::string newwriteOpts;
    std::ostringstream extraOpt;
    extraOpt  << ";PARALLEL_COMM=" << 1;
    newwriteOpts = writeOpts+extraOpt.str();
    result = mbImpl->write_file(outFile.c_str(), NULL, newwriteOpts.c_str(), partSets);
    PRINT_LAST_ERROR;
    std::cout << "Wrote " << outFile << std::endl;
    std::cout << "mbcoupler_test complete." << std::endl;
  }

  for (unsigned int i = 0; i < meshFiles.size(); i++) {
    delete pcs[i];
  }

  delete mbImpl;
  //may be leaking iMeshInst, don't care since it's end of program. Remove above deletes?

  err = MPI_Finalize();

  return 0;
}

ErrorCode report_iface_ents(Interface *mbImpl,
                              std::vector<ParallelComm *> &pcs,
                              const bool print_results)
{
  Range iface_ents[6];
  ErrorCode result = MB_SUCCESS, tmp_result;

    // now figure out which vertices are shared
  for (unsigned int p = 0; p < pcs.size(); p++) {
    for (int i = 0; i < 4; i++) {
      tmp_result = pcs[p]->get_iface_entities(-1, i, iface_ents[i]);

      if (MB_SUCCESS != tmp_result) {
        std::cerr << "get_iface_entities returned error on proc "
                  << pcs[p]->proc_config().proc_rank() << "; message: " << std::endl;
        std::string last_error;
        result = mbImpl->get_last_error(last_error);
        if (last_error.empty()) std::cerr << "(none)" << std::endl;
        else std::cerr << last_error << std::endl;
        result = tmp_result;
      }
      if (0 != i) iface_ents[4].merge(iface_ents[i]);
    }
  }

    // report # iface entities
  result = mbImpl->get_adjacencies(iface_ents[4], 0, false, iface_ents[5],
                                   Interface::UNION);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (print_results || iface_ents[0].size() != iface_ents[5].size()) {
    std::cerr << "Proc " << rank << " iface entities: " << std::endl;
    for (int i = 0; i < 4; i++)
      std::cerr << "    " << iface_ents[i].size() << " "
                << i << "d iface entities." << std::endl;
    std::cerr << "    (" << iface_ents[5].size()
              << " verts adj to other iface ents)" << std::endl;
  }

  return result;
}

// Print usage
void print_usage() {
  std::cerr << "Usage: ";
  std::cerr << "mbcoupler_test -meshes <source_mesh> <target_mesh> -itag <interp_tag> [-gnorm <gnorm_tag>] [-ssnorm <ssnorm_tag> <ssnorm_selection>] [-ropts <roptions>] [-outfile <out_file> [-wopts <woptions>]] [-dbgout [<dbg_file>]]" << std::endl;
  std::cerr << "    -meshes" << std::endl;
  std::cerr << "        Read in mesh files <source_mesh> and <target_mesh>." << std::endl;
  std::cerr << "    -itag" << std::endl;
  std::cerr << "        Interpolate tag <interp_tag> from source mesh to target mesh." << std::endl;
  std::cerr << "    -gnorm" << std::endl;
  std::cerr << "        Normalize the value of tag <gnorm_tag> over then entire mesh and save to" << std::endl;
  std::cerr << "        tag \"<gnorm_tag>_normf\" on the mesh set.  Do this for all meshes." << std::endl;
  std::cerr << "    -ssnorm" << std::endl;
  std::cerr << "        Normalize the value of tag <ssnorm_tag> over subsets of a mesh and save to" << std::endl;
  std::cerr << "        tag \"<ssnorm_tag>_normf\" on the Entity Set for each subset.  Subsets are selected" << std::endl;
  std::cerr << "        using criteria in <ssnorm_selection>.  Do this for all meshes." << std::endl;
  std::cerr << "    -ropts" << std::endl;
  std::cerr << "        Read in the mesh files using options in <roptions>." << std::endl;
  std::cerr << "    -outfile" << std::endl;
  std::cerr << "        Write out target mesh to <out_file>." << std::endl;
  std::cerr << "    -wopts" << std::endl;
  std::cerr << "        Write out mesh files using options in <woptions>." << std::endl;
  std::cerr << "    -dbgout" << std::endl;
  std::cerr << "        Write stdout and stderr streams to the file \'<dbg_file>.txt\'." << std::endl;
  std::cerr << "    -eps" << std::endl;
  std::cerr << "        epsilon" << std::endl;
  std::cerr << "    -meth <method> (0=CONSTANT, 1=LINEAR_FE, 2=QUADRATIC_FE, 3=SPECTRAL)" << std::endl;
}

// Check first character for a '-'.
// Return true if one is found.  False otherwise.
bool check_for_flag(const char *str) {
  if (str[0] == '-')
    return true;
  else
    return false;
}

// New get_file_options() function with added possibilities for mbcoupler_test.
ErrorCode get_file_options(int argc, char **argv,
                           std::vector<std::string> &meshFiles,
                           Coupler::Method &method,
                           std::string &interpTag,
                           std::string &gNormTag,
                           std::string &ssNormTag,
                           std::vector<const char *> &ssTagNames,
                           std::vector<const char *> &ssTagValues,
                           std::string &readOpts,
                           std::string &outFile,
                           std::string &writeOpts,
                           std::string &dbgFile,
                           bool &help,
                           double & epsilon)
{
  // Initialize some of the outputs to null values indicating not present
  // in the argument list.
  gNormTag = "";
  ssNormTag = "";
  readOpts = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARTITION_DISTRIBUTE;PARALLEL_RESOLVE_SHARED_ENTS;PARALLEL_GHOSTS=3.0.1;CPUTIME";
  outFile = "";
  writeOpts = "PARALLEL=WRITE_PART;CPUTIME";
  dbgFile = "";
  std::string defaultDbgFile = argv[0];  // The executable name will be the default debug output file.

  // These will indicate if we've gotten our required parameters at the end of parsing.
  bool haveMeshes = false;
  bool haveInterpTag = false;

  // Loop over the values in argv pulling out an parsing each one
  int npos = 1;

  if (argc > 1 && argv[1] == std::string("-h")) {
    help = true;
    return MB_SUCCESS;
  }

  while (npos < argc) {
    if (argv[npos] == std::string("-meshes")) {
      // Parse out the mesh filenames
      npos++;
      int numFiles = 2;
      meshFiles.resize(numFiles);
      for (int i = 0; i < numFiles; i++) {
        if ((npos < argc) && (!check_for_flag(argv[npos])))
          meshFiles[i] = argv[npos++];
        else {
          std::cerr << "    ERROR - missing correct number of mesh filenames" << std::endl;
          return MB_FAILURE;
        }
      }

      haveMeshes = true;
    }
    else if (argv[npos] == std::string("-itag")) {
      // Parse out the interpolation tag
      npos++;
      if ((npos < argc) && (!check_for_flag(argv[npos])))
        interpTag = argv[npos++];
      else {
        std::cerr << "    ERROR - missing <interp_tag>" << std::endl;
        return MB_FAILURE;
      }

      haveInterpTag = true;
    }
    else if (argv[npos] == std::string("-meth")) {
      // Parse out the interpolation tag
      npos++;
      if (argv[npos][0] == '0') method = Coupler::CONSTANT;
      else if (argv[npos][0] == '1') method = Coupler::LINEAR_FE;
      else if (argv[npos][0] == '2') method = Coupler::QUADRATIC_FE;
      else if (argv[npos][0] == '3') method = Coupler::SPECTRAL;
      else {
        std::cerr << "    ERROR - unrecognized method number " << method << std::endl;
        return MB_FAILURE;
      }
      npos++;
    }
    else if (argv[npos] == std::string("-eps")) {
      // Parse out the tolerance
      npos++;
      if ((npos < argc) && (!check_for_flag(argv[npos])))
        epsilon = atof(argv[npos++]);
      else {
        std::cerr << "    ERROR - missing <epsilon>" << std::endl;
        return MB_FAILURE;
      }
    }
    else if (argv[npos] == std::string("-gnorm")) {
      // Parse out the global normalization tag
      npos++;
      if ((npos < argc) && (!check_for_flag(argv[npos])))
        gNormTag = argv[npos++];
      else {
        std::cerr << "    ERROR - missing <gnorm_tag>" << std::endl;
        return MB_FAILURE;
      }
    }
    else if (argv[npos] == std::string("-ssnorm")) {
      // Parse out the subset normalization tag and selection criteria
      npos++;
      if ((npos < argc) && (!check_for_flag(argv[npos])))
        ssNormTag = argv[npos++];
      else {
        std::cerr << "    ERROR - missing <ssnorm_tag>" << std::endl;
        return MB_FAILURE;
      }

      if ((npos < argc) && (!check_for_flag(argv[npos]))) {
        char* opts = argv[npos++];
        char sep1[1] = {';'};
        char sep2[1] = {'='};
        bool end_vals_seen = false;
        std::vector<char *> tmpTagOpts;

        // first get the options
        for (char* i = strtok(opts, sep1); i; i = strtok(0, sep1)) {
          tmpTagOpts.push_back(i);
        }

        // parse out the name and val or just name.
        for (unsigned int j = 0; j < tmpTagOpts.size(); j++) {
          char* e = strtok(tmpTagOpts[j], sep2);
          ssTagNames.push_back(e);
          e = strtok(0, sep2);
          if (e != NULL) {
            // We have a value
            if (end_vals_seen) {
              // ERROR we should not have a value after none are seen
              std::cerr << "    ERROR - new value seen after end of values in <ssnorm_selection>" << std::endl;
              return MB_FAILURE;
            }
            // Otherwise get the value string from e and convert it to an int
            int *valp = new int;
            *valp = atoi(e);
            ssTagValues.push_back((const char *) valp);
          }
          else {
            // Otherwise there is no '=' so push a null on the list
            end_vals_seen = true;
            ssTagValues.push_back((const char *) 0);
          }
        }
      }
      else {
        std::cerr << "    ERROR - missing <ssnorm_selection>" << std::endl;
        return MB_FAILURE;
      }
    }
    else if (argv[npos] == std::string("-ropts")) {
      // Parse out the mesh file read options
      npos++;
      if ((npos < argc) && (!check_for_flag(argv[npos])))
        readOpts = argv[npos++];
      else {
        std::cerr << "    ERROR - missing <roptions>" << std::endl;
        return MB_FAILURE;
      }
    }
    else if (argv[npos] == std::string("-outfile")) {
      // Parse out the output file name
      npos++;
      if ((npos < argc) && (!check_for_flag(argv[npos])))
        outFile = argv[npos++];
      else {
        std::cerr << "    ERROR - missing <out_file>" << std::endl;
        return MB_FAILURE;
      }
    }
    else if (argv[npos] == std::string("-wopts")) {
      // Parse out the output file write options
      npos++;
      if ((npos < argc) && (!check_for_flag(argv[npos])))
        writeOpts = argv[npos++];
      else {
        std::cerr << "    ERROR - missing <woptions>" << std::endl;
        return MB_FAILURE;
      }
    }
    else if (argv[npos] == std::string("-dbgout")) {
      // Parse out the debug output file name.
      // If no name then use the default.
      npos++;
      if ((npos < argc) && (!check_for_flag(argv[npos])))
        dbgFile = argv[npos++];
      else
        dbgFile = defaultDbgFile;
    }
    else {
      // Unrecognized parameter.  Skip it and move along.
      std::cerr << "    ERROR - Unrecognized parameter:" << argv[npos] << std::endl;
      std::cerr << "            Skipping..." << std::endl;
      npos++;
    }
  }

  if (!haveMeshes) {
    meshFiles.resize(2);
    meshFiles[0] = std::string(TestDir + "/64bricks_1khex.h5m");
    meshFiles[1] = std::string(TestDir + "/64bricks_12ktet.h5m");
    std::cout << "Mesh files not entered; using default files "
              << meshFiles[0] << " and " << meshFiles[1] << std::endl;
  }

  if (!haveInterpTag) {
    interpTag = "vertex_field";
    std::cout << "Interpolation field name not given, using default of " << interpTag << std::endl;
  }

#ifdef HDF5_FILE
  if (1 == argc) {
    std::cout << "No arguments given; using output file dum.h5m." << std::endl;
    outFile = "dum.h5m";
  }
#endif

  return MB_SUCCESS;
}


// End new get_file_options()

// ErrorCode get_file_options(int argc, char **argv,
//                              std::vector<const char *> &filenames,
//                              std::string &interp_tag,
//                              std::string &out_fname,
//                              std::string &opts)
// {
//   int npos = 1;
//   int nfiles = atoi(argv[npos++]);

//     // get mesh filenames
//   filenames.resize(nfiles);
//   for (int i = 0; i < nfiles; i++) filenames[i] = argv[npos++];

//   if (npos < argc) interp_tag = argv[npos++];

//   if (npos < argc) out_fname = argv[npos++];

//     // get partition information
//   const char *tag_name = "GEOM_DIMENSION";
//   int tag_val = -1;
//   int distrib = 1;
//   int with_ghosts = 0;
//   if (npos < argc) tag_name = argv[npos++];
//   if (npos < argc) tag_val = strtol(argv[npos++], NULL, 0);
//   if (npos < argc) distrib = strtol(argv[npos++], NULL, 0);
//   if (npos < argc) with_ghosts = strtol(argv[npos++], NULL, 0);

//   if (-1 == tag_val && !strcmp(tag_name, "GEOM_DIMENSION"))
//     tag_val = 3;

//   std::ostringstream options;
//   options << "PARALLEL=READ_DELETE;PARTITION=" << tag_name;

//   if (-1 != tag_val)
//     options << ";PARTITION_VAL=" << tag_val;

//   if (1 == distrib)
//     options << ";PARTITION_DISTRIBUTE";

//   options << ";PARALLEL_RESOLVE_SHARED_ENTS";

//   if (1 == with_ghosts)
//     options << ";PARALLEL_GHOSTS=3.0.1";

//   options << ";CPUTIME";

//   opts = options.str();

//   return MB_SUCCESS;
// }

ErrorCode test_interpolation(Interface *mbImpl,
                             Coupler::Method method,
                             std::string &interpTag,
                             std::string &gNormTag,
                             std::string &ssNormTag,
                             std::vector<const char *> &ssTagNames,
                             std::vector<const char *> &ssTagValues,
                             iBase_EntitySetHandle *roots,
                             std::vector<ParallelComm *> &pcs,
                             double &instant_time,
                             double &pointloc_time,
                             double &interp_time,
                             double &gnorm_time,
                             double &ssnorm_time,
                             double & toler)
{
  assert(method >= Coupler::CONSTANT && method <= Coupler::SPECTRAL);

    // source is 1st mesh, target is 2nd
  Range src_elems, targ_elems, targ_verts;
  ErrorCode result = pcs[0]->get_part_entities(src_elems, 3);
  PRINT_LAST_ERROR;

  double start_time = MPI_Wtime();

    // instantiate a coupler, which also initializes the tree
  Coupler mbc(mbImpl, pcs[0], src_elems, 0);

  // initialize spectral elements, if they exist
  bool specSou=false, specTar = false;
//  result =  mbc.initialize_spectral_elements((EntityHandle)roots[0], (EntityHandle)roots[1], specSou, specTar);


  instant_time = MPI_Wtime();

    // get points from the target mesh to interpolate
  // we have to treat differently the case when the target is a spectral mesh
  // in that case, the points of interest are the GL points, not the vertex nodes
  std::vector<double> vpos; // this will have the positions we are interested in
  int numPointsOfInterest = 0;
  if (!specTar)
  {// usual case
    Range tmp_verts;

      // first get all vertices adj to partition entities in target mesh
    result = pcs[1]->get_part_entities(targ_elems, 3);
    PRINT_LAST_ERROR;
    if (Coupler::CONSTANT == method)
      targ_verts = targ_elems;
    else
      result = mbImpl->get_adjacencies(targ_elems, 0, false, targ_verts,
                                       Interface::UNION);
    PRINT_LAST_ERROR;

      // then get non-owned verts and subtract
    result = pcs[1]->get_pstatus_entities(0, PSTATUS_NOT_OWNED, tmp_verts);
    PRINT_LAST_ERROR;
    targ_verts = subtract( targ_verts, tmp_verts);
    // get position of these entities; these are the target points
    numPointsOfInterest = (int)targ_verts.size();
    vpos.resize(3*targ_verts.size());
    result = mbImpl->get_coords(targ_verts, &vpos[0]);
    PRINT_LAST_ERROR;
  }
  else
  {
    // in this case, the target mesh is spectral, we want values
    // interpolated on the GL positions; for each element, get the GL points, and construct CartVect!!!
    result = pcs[1]->get_part_entities(targ_elems, 3);
    PRINT_LAST_ERROR;
    result = mbc.get_gl_points_on_elements(targ_elems, vpos, numPointsOfInterest);
    PRINT_LAST_ERROR;
  }




    // locate those points in the source mesh
  std::cout<<"rank "<< pcs[0]->proc_config().proc_rank() << " points of interest: " << numPointsOfInterest << "\n";
  result = mbc.locate_points(&vpos[0], numPointsOfInterest, 0, toler);
  PRINT_LAST_ERROR;

  pointloc_time = MPI_Wtime();

    // now interpolate tag onto target points
  std::vector<double> field(numPointsOfInterest);

  result = mbc.interpolate(method, interpTag, &field[0]);
  PRINT_LAST_ERROR;

  interp_time = MPI_Wtime();

    // do global normalization if specified
  if (!gNormTag.empty()) {
    int err;

    // Normalize the source mesh
    err = mbc.normalize_mesh(roots[0], gNormTag.c_str(), Coupler::VOLUME, 4);
    PRINT_LAST_ERR;

    // Normalize the target mesh
    err = mbc.normalize_mesh(roots[1], gNormTag.c_str(), Coupler::VOLUME, 4);
    PRINT_LAST_ERR;
  }

  gnorm_time = MPI_Wtime();

    // do subset normalization if specified

  if (!ssNormTag.empty()) {
    int err;

    err = mbc.normalize_subset(roots[0],
                               ssNormTag.c_str(),
                               &ssTagNames[0],
                               ssTagNames.size(),
                               &ssTagValues[0],
                               Coupler::VOLUME,
                               4);
    PRINT_LAST_ERR;

    err = mbc.normalize_subset(roots[1],
                               ssNormTag.c_str(),
                               &ssTagNames[0],
                               ssTagNames.size(),
                               &ssTagValues[0],
                               Coupler::VOLUME,
                               4);
    PRINT_LAST_ERR;
  }

  ssnorm_time = MPI_Wtime();

  ssnorm_time -= gnorm_time;
  gnorm_time -= interp_time;
  interp_time -= pointloc_time;
  pointloc_time -= instant_time;
  instant_time -= start_time;

    // set field values as tag on target vertices
  if (specSou)
  {
    // create a new tag for the values on the target
    Tag tag;
    std::string newtag = interpTag +"_TAR";
    result = mbImpl->tag_get_handle(newtag.c_str(), 1, MB_TYPE_DOUBLE, tag, MB_TAG_CREAT|MB_TAG_DENSE);
    result = mbImpl->tag_set_data(tag, targ_verts, &field[0]);
    PRINT_LAST_ERROR;
  }
  else
  {
    if (!specTar)
    {
      // use original tag
      Tag tag;
      result = mbImpl->tag_get_handle(interpTag.c_str(), 1, MB_TYPE_DOUBLE, tag); PRINT_LAST_ERROR;
      result = mbImpl->tag_set_data(tag, targ_verts, &field[0]);
      PRINT_LAST_ERROR;
    }
    else
    {
      // we have the field values computed at each GL points, on each element
      // in the target mesh
      // we need to create a new tag, on elements, of size _ntot, to hold
      // all those values.
      // so it turns out we need ntot. maybe we can compute it from the
      // number of values computed, divided by number of elements
      int ntot = numPointsOfInterest/targ_elems.size();
      Tag tag;
      std::string newtag = interpTag +"_TAR";
      result = mbImpl->tag_get_handle(newtag.c_str(), ntot, MB_TYPE_DOUBLE, tag, MB_TAG_CREAT|MB_TAG_DENSE);
      result = mbImpl->tag_set_data(tag, targ_elems, &field[0]);
      PRINT_LAST_ERROR;

    }
  }
    // done
  return MB_SUCCESS;
}