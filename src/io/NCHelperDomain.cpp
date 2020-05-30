#include "NCHelperDomain.hpp"
#include "moab/FileOptions.hpp"

#include <cmath>
#include <sstream>

namespace moab {

bool NCHelperDomain::can_read_file(ReadNC* readNC, int fileId)
{
  std::vector<std::string>& dimNames = readNC->dimNames;

  // If dimension names "n" AND "ni" AND "nj" AND "nv" exist then it should be the Domain grid
  if ((std::find(dimNames.begin(), dimNames.end(), std::string("n")) != dimNames.end()) && (std::find(dimNames.begin(),
      dimNames.end(), std::string("ni")) != dimNames.end()) && (std::find(dimNames.begin(), dimNames.end(), std::string("nj"))
      != dimNames.end()) && (std::find(dimNames.begin(), dimNames.end(), std::string("nv")) != dimNames.end())) {
    // Make sure it is CAM grid
    std::map<std::string, ReadNC::AttData>::iterator attIt = readNC->globalAtts.find("source");
    if (attIt == readNC->globalAtts.end())
      return false;
    unsigned int sz = attIt->second.attLen;
    std::string att_data;
    att_data.resize(sz + 1);
    att_data[sz] = '\000';
    int success = NCFUNC(get_att_text)(fileId, attIt->second.attVarId, attIt->second.attName.c_str(), &att_data[0]);
    if (success)
      return false;
    /*if (att_data.find("CAM") == std::string::npos)
      return false;*/

    return true;
  }

  return false;
}

ErrorCode NCHelperDomain::init_mesh_vals()
{
  Interface*& mbImpl = _readNC->mbImpl;
  std::vector<std::string>& dimNames = _readNC->dimNames;
  std::vector<int>& dimLens = _readNC->dimLens;
  std::map<std::string, ReadNC::VarData>& varInfo = _readNC->varInfo;
  DebugOutput& dbgOut = _readNC->dbgOut;
  bool& isParallel = _readNC->isParallel;
  int& partMethod = _readNC->partMethod;
  ScdParData& parData = _readNC->parData;

  ErrorCode rval;

  // Look for names of i/j dimensions
  // First i
  std::vector<std::string>::iterator vit;
  unsigned int idx;
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "ni")) != dimNames.end())
    idx = vit - dimNames.begin();
  else {
    MB_SET_ERR(MB_FAILURE, "Couldn't find 'ni' variable");
  }
  iDim = idx;
  gDims[0] = 0;
  gDims[3] = dimLens[idx] - 1;

  // Then j
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "nj")) != dimNames.end())
    idx = vit - dimNames.begin();
  else {
    MB_SET_ERR(MB_FAILURE, "Couldn't find 'nj' variable");
  }
  jDim = idx;
  gDims[1] = 0;
  gDims[4] = dimLens[idx] - 1 ; // Add 2 for the pole points ? not needed



 // do not use gcdims ? or use only gcdims?

  // Try a truly 2D mesh
  gDims[2] = -1;
  gDims[5] = -1;

  // Get number of vertices per cell
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "nv")) != dimNames.end())
    idx = vit - dimNames.begin();
  else {
    MB_SET_ERR(MB_FAILURE, "Couldn't find 'nv' dimension");
  }
  nvDim = idx;
  nv = dimLens[idx];

  // Parse options to get subset
  int rank = 0, procs = 1;
#ifdef MOAB_HAVE_MPI
  if (isParallel) {
    ParallelComm*& myPcomm = _readNC->myPcomm;
    rank = myPcomm->proc_config().proc_rank();
    procs = myPcomm->proc_config().proc_size();
  }
#endif
  if (procs > 1) {
    for (int i = 0; i < 6; i++)
      parData.gDims[i] = gDims[i];
    parData.partMethod = partMethod;
    int pdims[3];

    locallyPeriodic[0] = locallyPeriodic[1] = locallyPeriodic[2] = 0;
    rval = ScdInterface::compute_partition(procs, rank, parData, lDims, locallyPeriodic, pdims);MB_CHK_ERR(rval);
    for (int i = 0; i < 3; i++)
      parData.pDims[i] = pdims[i];

    dbgOut.tprintf(1, "Partition: %dx%d (out of %dx%d)\n",
        lDims[3] - lDims[0] + 1, lDims[4] - lDims[1] + 1,
        gDims[3] - gDims[0] + 1, gDims[4] - gDims[1] + 1);
    if (0 == rank)
      dbgOut.tprintf(1, "Contiguous chunks of size %d bytes.\n", 8 * (lDims[3] - lDims[0] + 1) * (lDims[4] - lDims[1] + 1));
  }
  else {
    for (int i = 0; i < 6; i++)
      lDims[i] = gDims[i];
    locallyPeriodic[0] = globallyPeriodic[0];
  }


  // Now get actual coordinate values for vertices and cell centers
  lCDims[0] = lDims[0];

  lCDims[3] = lDims[3];


  // For FV models, will always be non-periodic in j
  lCDims[1] = lDims[1];
  lCDims[4] = lDims[4] - 1;

#if 0
  // Resize vectors to store values later
  if (-1 != lDims[0])
    ilVals.resize(lDims[3] - lDims[0] + 1);
  if (-1 != lCDims[0])
    ilCVals.resize(lCDims[3] - lCDims[0] + 1);
  if (-1 != lDims[1])
    jlVals.resize(lDims[4] - lDims[1] + 1);
  if (-1 != lCDims[1])
    jlCVals.resize(lCDims[4] - lCDims[1] + 1);
  if (nTimeSteps > 0)
    tVals.resize(nTimeSteps);
#endif


  dbgOut.tprintf(1, "I=%d-%d, J=%d-%d\n", lDims[0], lDims[3], lDims[1], lDims[4]);
  dbgOut.tprintf(1, "%d elements, %d vertices\n", (lDims[3] - lDims[0]) * (lDims[4] - lDims[1]), (lDims[3] - lDims[0]) * (lDims[4] - lDims[1]) * nv);

  // For each variable, determine the entity location type and number of levels
  std::map<std::string, ReadNC::VarData>::iterator mit;
  for (mit = varInfo.begin(); mit != varInfo.end(); ++mit) {
    ReadNC::VarData& vd = (*mit).second;

    // Default entLoc is ENTLOCSET
    if (std::find(vd.varDims.begin(), vd.varDims.end(), tDim) != vd.varDims.end()) {
      if ((std::find(vd.varDims.begin(), vd.varDims.end(), iCDim) != vd.varDims.end()) &&
          (std::find(vd.varDims.begin(), vd.varDims.end(), jCDim) != vd.varDims.end()))
        vd.entLoc = ReadNC::ENTLOCFACE;
      else if ((std::find(vd.varDims.begin(), vd.varDims.end(), jDim) != vd.varDims.end()) &&
          (std::find(vd.varDims.begin(), vd.varDims.end(), iCDim) != vd.varDims.end()))
        vd.entLoc = ReadNC::ENTLOCNSEDGE;
      else if ((std::find(vd.varDims.begin(), vd.varDims.end(), jCDim) != vd.varDims.end()) &&
          (std::find(vd.varDims.begin(), vd.varDims.end(), iDim) != vd.varDims.end()))
        vd.entLoc = ReadNC::ENTLOCEWEDGE;
    }

    // Default numLev is 0
    if (std::find(vd.varDims.begin(), vd.varDims.end(), levDim) != vd.varDims.end())
      vd.numLev = nLevels;
  }

  std::vector<std::string> ijdimNames(2);
  ijdimNames[0] = "__ni";
  ijdimNames[1] = "__nj";
  std::string tag_name;
  Tag tagh;

  // __<dim_name>_LOC_MINMAX (for slon, slat, lon and lat)
  for (unsigned int i = 0; i != ijdimNames.size(); i++) {
    std::vector<int> val(2, 0);
    if (ijdimNames[i] == "__ni") {
      val[0] = lDims[0];
      val[1] = lDims[3];
    }
    else if (ijdimNames[i] == "__nj") {
      val[0] = lDims[1];
      val[1] = lDims[4];
    }

    std::stringstream ss_tag_name;
    ss_tag_name << ijdimNames[i] << "_LOC_MINMAX";
    tag_name = ss_tag_name.str();
    rval = mbImpl->tag_get_handle(tag_name.c_str(), 2, MB_TYPE_INTEGER, tagh,
                                  MB_TAG_SPARSE | MB_TAG_CREAT);MB_CHK_SET_ERR(rval, "Trouble creating conventional tag " << tag_name);
    rval = mbImpl->tag_set_data(tagh, &_fileSet, 1, &val[0]);MB_CHK_SET_ERR(rval, "Trouble setting data to conventional tag " << tag_name);
    if (MB_SUCCESS == rval)
      dbgOut.tprintf(2, "Conventional tag %s is created.\n", tag_name.c_str());
  }

  // __<dim_name>_LOC_VALS (for slon, slat, lon and lat)
  // Assume all have the same data type as lon (expected type is float or double)
  switch (varInfo["xc"].varDataType) {
    case NC_FLOAT:
    case NC_DOUBLE:
      break;
    default:
      MB_SET_ERR(MB_FAILURE, "Unexpected variable data type for 'lon'");
  }

  // do not need conventional tags
  Tag convTagsCreated = 0;
  int def_val = 0;
  rval = mbImpl->tag_get_handle("__CONV_TAGS_CREATED", 1, MB_TYPE_INTEGER, convTagsCreated,
                                MB_TAG_SPARSE | MB_TAG_CREAT, &def_val);MB_CHK_SET_ERR(rval, "Trouble getting _CONV_TAGS_CREATED tag");
  int create_conv_tags_flag = 1;
  rval = mbImpl->tag_set_data(convTagsCreated, &_fileSet, 1, &create_conv_tags_flag);MB_CHK_SET_ERR(rval, "Trouble setting _CONV_TAGS_CREATED tag");

  return MB_SUCCESS;
}

ErrorCode NCHelperDomain::create_mesh(Range& faces)
{
  Interface*& mbImpl = _readNC->mbImpl;
  std::string& fileName = _readNC->fileName;
  Tag& mGlobalIdTag = _readNC->mGlobalIdTag;
  const Tag*& mpFileIdTag = _readNC->mpFileIdTag;
  DebugOutput& dbgOut = _readNC->dbgOut;
  int& gatherSetRank = _readNC->gatherSetRank;
  int& trivialPartitionShift = _readNC->trivialPartitionShift;

  int rank = 0;
  int procs = 1;
#ifdef MOAB_HAVE_MPI
  bool& isParallel = _readNC->isParallel;
  if (isParallel) {
    ParallelComm*& myPcomm = _readNC->myPcomm;
    rank = myPcomm->proc_config().proc_rank();
    procs = myPcomm->proc_config().proc_size();
  }
#endif

  ErrorCode rval;
  int success = 0;


  bool create_gathers = false;
  if (rank == gatherSetRank)
    create_gathers = true;

  // Shift rank to obtain a rotated trivial partition
  int shifted_rank = rank;
  if (procs >= 2 && trivialPartitionShift > 0)
    shifted_rank = (rank + trivialPartitionShift) % procs;

#if 0

  int local_elems =
  // num_coarse_quads is the number of elems instantiated in MOAB; if !spectralMesh, num_coarse_quads = num_fine_quads
  num_elems = int(std::floor(1.0 * num_quads / (spectral_unit * procs)));
  // start_idx is the starting index in the HommeMapping connectivity list for this proc, before converting to coarse quad representation
  start_idx = 4 * shifted_rank * num_coarse_quads * spectral_unit;
  // iextra = # coarse quads extra after equal split over procs
  int iextra = num_quads % (procs * spectral_unit);
  if (shifted_rank < iextra)
    num_coarse_quads++;
  start_idx += 4 * spectral_unit * std::min(shifted_rank, iextra);
  // num_fine_quads is the number of quads in the connectivity list in HommeMapping file assigned to this proc
  num_fine_quads = spectral_unit * num_coarse_quads;

  // Now create num_coarse_quads
  EntityHandle* conn_arr;
  EntityHandle start_vertex;
  Range tmp_range;

  // Read connectivity into that space
  EntityHandle* sv_ptr = NULL;
  EntityHandle start_quad;
  EntityType mdb_type = MBTRI;
  if (nv==3)
    mdb_type = MBTRI;
  else if (nv==4)
    mdb_type = MBQUAD;
  else // (nv > 4)
    mdb_type = MBPOLYGON;

  if (!spectralMesh) {
    rval = _readNC->readMeshIface->get_element_connect(num_coarse_quads, 4,
                                                      MBQUAD, 0, start_quad, conn_arr,
                                                      // Might have to create gather mesh later
                                                      (create_gathers ? num_coarse_quads + num_quads : num_coarse_quads));MB_CHK_SET_ERR(rval, "Failed to create local quads");
    tmp_range.insert(start_quad, start_quad + num_coarse_quads - 1);
    int* tmp_conn_end = (&tmp_conn[start_idx + 4 * num_fine_quads-1])+1;
    std::copy(&tmp_conn[start_idx], tmp_conn_end, conn_arr);
    std::copy(conn_arr, conn_arr + 4 * num_fine_quads, range_inserter(localGidVerts));
  }
  else {
    rval = smt.create_spectral_elems(&tmp_conn[0], num_fine_quads, 2, tmp_range, start_idx, &localGidVerts);MB_CHK_SET_ERR(rval, "Failed to create spectral elements");
    int count, v_per_e;
    rval = mbImpl->connect_iterate(tmp_range.begin(), tmp_range.end(), conn_arr, v_per_e, count);MB_CHK_SET_ERR(rval, "Failed to get connectivity of spectral elements");
    rval = mbImpl->tag_iterate(smt.spectral_vertices_tag(true), tmp_range.begin(), tmp_range.end(),
                               count, (void*&)sv_ptr);MB_CHK_SET_ERR(rval, "Failed to get fine connectivity of spectral elements");
  }

  // Create vertices
  nLocalVertices = localGidVerts.size();
  std::vector<double*> arrays;
  rval = _readNC->readMeshIface->get_node_coords(3, nLocalVertices, 0, start_vertex, arrays,
                                                // Might have to create gather mesh later
                                                (create_gathers ? nLocalVertices + nVertices : nLocalVertices));MB_CHK_SET_ERR(rval, "Failed to create local vertices");

  // Set vertex coordinates
  Range::iterator rit;
  double* xptr = arrays[0];
  double* yptr = arrays[1];
  double* zptr = arrays[2];
  int i;
  for (i = 0, rit = localGidVerts.begin(); i < nLocalVertices; i++, ++rit) {
    assert(*rit < xVertVals.size() + 1);
    xptr[i] = xVertVals[(*rit) - 1]; // lon
    yptr[i] = yVertVals[(*rit) - 1]; // lat
  }

  // Convert lon/lat/rad to x/y/z
  const double pideg = acos(-1.0) / 180.0;
  double rad = (isConnFile) ? 8000.0 : 8000.0 + levVals[0];
  for (i = 0; i < nLocalVertices; i++) {
    double cosphi = cos(pideg * yptr[i]);
    double zmult = sin(pideg * yptr[i]);
    double xmult = cosphi * cos(xptr[i] * pideg);
    double ymult = cosphi * sin(xptr[i] * pideg);
    xptr[i] = rad * xmult;
    yptr[i] = rad * ymult;
    zptr[i] = rad * zmult;
  }

  // Get ptr to gid memory for vertices
  Range vert_range(start_vertex, start_vertex + nLocalVertices - 1);
  void* data;
  int count;
  rval = mbImpl->tag_iterate(mGlobalIdTag, vert_range.begin(), vert_range.end(),
                             count, data);MB_CHK_SET_ERR(rval, "Failed to iterate global id tag on local vertices");
  assert(count == nLocalVertices);
  int* gid_data = (int*) data;
  std::copy(localGidVerts.begin(), localGidVerts.end(), gid_data);

  // Duplicate global id data, which will be used to resolve sharing
  if (mpFileIdTag) {
    rval = mbImpl->tag_iterate(*mpFileIdTag, vert_range.begin(), vert_range.end(),
                               count, data);MB_CHK_SET_ERR(rval, "Failed to iterate file id tag on local vertices");
    assert(count == nLocalVertices);
    int bytes_per_tag = 4;
    rval = mbImpl->tag_get_bytes(*mpFileIdTag, bytes_per_tag);MB_CHK_SET_ERR(rval, "Can't get number of bytes for file id tag");
    if (4 == bytes_per_tag) {
      gid_data = (int*) data;
      std::copy(localGidVerts.begin(), localGidVerts.end(), gid_data);
    }
    else if (8 == bytes_per_tag) { // Should be a handle tag on 64 bit machine?
      long* handle_tag_data = (long*)data;
      std::copy(localGidVerts.begin(), localGidVerts.end(), handle_tag_data);
    }
  }

  // Create map from file ids to vertex handles, used later to set connectivity
  std::map<EntityHandle, EntityHandle> vert_handles;
  for (rit = localGidVerts.begin(), i = 0; rit != localGidVerts.end(); ++rit, i++)
    vert_handles[*rit] = start_vertex + i;

  // Compute proper handles in connectivity using offset
  for (int q = 0; q < 4 * num_coarse_quads; q++) {
    conn_arr[q] = vert_handles[conn_arr[q]];
    assert(conn_arr[q]);
  }
  if (spectralMesh) {
    int verts_per_quad = (_spectralOrder + 1) * (_spectralOrder + 1);
    for (int q = 0; q < verts_per_quad * num_coarse_quads; q++) {
      sv_ptr[q] = vert_handles[sv_ptr[q]];
      assert(sv_ptr[q]);
    }
  }

  // Add new vertices and quads to current file set
  faces.merge(tmp_range);
  tmp_range.insert(start_vertex, start_vertex + nLocalVertices - 1);
  rval = mbImpl->add_entities(_fileSet, tmp_range);MB_CHK_SET_ERR(rval, "Failed to add new vertices and quads to current file set");

  // Mark the set with the spectral order
  Tag sporder;
  rval = mbImpl->tag_get_handle("SPECTRAL_ORDER", 1, MB_TYPE_INTEGER, sporder,
                                MB_TAG_SPARSE | MB_TAG_CREAT);MB_CHK_SET_ERR(rval, "Trouble creating SPECTRAL_ORDER tag");
  rval = mbImpl->tag_set_data(sporder, &_fileSet, 1, &_spectralOrder);MB_CHK_SET_ERR(rval, "Trouble setting data to SPECTRAL_ORDER tag");

  if (create_gathers) {
    EntityHandle gather_set;
    rval = _readNC->readMeshIface->create_gather_set(gather_set);MB_CHK_SET_ERR(rval, "Failed to create gather set");

    // Create vertices
    arrays.clear();
    // Don't need to specify allocation number here, because we know enough verts were created before
    rval = _readNC->readMeshIface->get_node_coords(3, nVertices, 0, start_vertex, arrays);MB_CHK_SET_ERR(rval, "Failed to create gather set vertices");

    xptr = arrays[0];
    yptr = arrays[1];
    zptr = arrays[2];
    for (i = 0; i < nVertices; i++) {
      double cosphi = cos(pideg * yVertVals[i]);
      double zmult = sin(pideg * yVertVals[i]);
      double xmult = cosphi * cos(xVertVals[i] * pideg);
      double ymult = cosphi * sin(xVertVals[i] * pideg);
      xptr[i] = rad * xmult;
      yptr[i] = rad * ymult;
      zptr[i] = rad * zmult;
    }

    // Get ptr to gid memory for vertices
    Range gather_set_verts_range(start_vertex, start_vertex + nVertices - 1);
    rval = mbImpl->tag_iterate(mGlobalIdTag, gather_set_verts_range.begin(), gather_set_verts_range.end(),
                               count, data);MB_CHK_SET_ERR(rval, "Failed to iterate global id tag on gather set vertices");
    assert(count == nVertices);
    gid_data = (int*) data;
    for (int j = 1; j <= nVertices; j++)
      gid_data[j - 1] = j;
    // Set the file id tag too, it should be bigger something not interfering with global id
    if (mpFileIdTag) {
      rval = mbImpl->tag_iterate(*mpFileIdTag, gather_set_verts_range.begin(), gather_set_verts_range.end(),
                                 count, data);MB_CHK_SET_ERR(rval, "Failed to iterate file id tag on gather set vertices");
      assert(count == nVertices);
      int bytes_per_tag = 4;
      rval = mbImpl->tag_get_bytes(*mpFileIdTag, bytes_per_tag);MB_CHK_SET_ERR(rval, "Can't get number of bytes for file id tag");
      if (4 == bytes_per_tag) {
        gid_data = (int*)data;
        for (int j = 1; j <= nVertices; j++)
          gid_data[j - 1] = nVertices + j; // Bigger than global id tag
      }
      else if (8 == bytes_per_tag) { // Should be a handle tag on 64 bit machine?
        long* handle_tag_data = (long*)data;
        for (int j = 1; j <= nVertices; j++)
          handle_tag_data[j - 1] = nVertices + j; // Bigger than global id tag
      }
    }

    rval = mbImpl->add_entities(gather_set, gather_set_verts_range);MB_CHK_SET_ERR(rval, "Failed to add vertices to the gather set");

    // Create quads
    Range gather_set_quads_range;
    // Don't need to specify allocation number here, because we know enough quads were created before
    rval = _readNC->readMeshIface->get_element_connect(num_quads, 4, MBQUAD, 0,
                                                       start_quad, conn_arr);MB_CHK_SET_ERR(rval, "Failed to create gather set quads");
    gather_set_quads_range.insert(start_quad, start_quad + num_quads - 1);
    int* tmp_conn_end = (&tmp_conn[4 * num_quads-1]) + 1;
    std::copy(&tmp_conn[0], tmp_conn_end, conn_arr);
    for (i = 0; i != 4 * num_quads; i++)
      conn_arr[i] += start_vertex - 1; // Connectivity array is shifted by where the gather verts start
    rval = mbImpl->add_entities(gather_set, gather_set_quads_range);MB_CHK_SET_ERR(rval, "Failed to add quads to the gather set");
  }

#endif
  return MB_SUCCESS;

}
} // namespace moab
