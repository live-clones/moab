// Usage:
// tools/h5mtoscrip -w map_atm_to_ocn.h5m -s map_atm_to_ocn.nc --coords
// 
#include <iostream>
#include <exception>
#include <cmath>
#include <cassert>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>

#include "moab/ProgOptions.hpp"
#include "moab/Core.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#endif

#ifndef MOAB_HAVE_TEMPESTREMAP
#error Tool requires compilation with TempestRemap dependency
#endif

// TempestRemap includes
#include "OfflineMap.h"
#include "netcdfcpp.h"
#include "NetCDFUtilities.h"
#include "DataArray2D.h"

using namespace moab;

template<typename T>
ErrorCode get_vartag_data(moab::Interface* mbCore, Tag tag, moab::Range& sets, int& out_data_size, std::vector<T>& data)
{
  int* tag_sizes = new int [ sets.size() ];
  const void** tag_data = (const void**) new void* [sets.size()];

  ErrorCode rval = mbCore->tag_get_by_ptr(tag, sets, tag_data, tag_sizes );MB_CHK_SET_ERR(rval, "Getting matrix rows failed");

  out_data_size = 0;
  for (unsigned is=0; is < sets.size(); ++is)
    out_data_size += tag_sizes[is];

  data.resize(out_data_size);
  int ioffset = 0;
  for (unsigned index=0; index < sets.size(); index++)
  {
    T* m_vals = (T*)tag_data[index];
    for (int k=0; k < tag_sizes[index]; k++)
    {
      data[ioffset++] = m_vals[k];
    }
  }

  return moab::MB_SUCCESS;
}

void ReadFileMetaData(std::string& metaFilename, std::map<std::string, std::string>& metadataVals)
{
  std::ifstream metafile;
  std::string line;

  metafile.open(metaFilename.c_str());
  metadataVals["Title"] = "MOAB-TempestRemap (MBTR) Offline Regridding Weight Converter (h5mtoscrip)";
  std::string key, value;
  while (std::getline(metafile, line))
  {
    size_t lastindex = line.find_last_of("=");
    key = line.substr(0, lastindex-1);
    value = line.substr(lastindex+2, line.length());

    metadataVals[std::string(key)] = std::string(value);
  }
  metafile.close();
}

int main(int argc, char* argv[])
{
  moab::ErrorCode rval;
  int dimension=2;
  NcError error2 ( NcError::verbose_nonfatal );
  std::stringstream sstr;
  ProgOptions opts;
  std::string h5mfilename, scripfile;
  bool noMap = false;
  bool writeXYCoords = false;

#ifdef MOAB_HAVE_MPI
    MPI_Init ( &argc, &argv );
#endif

  opts.addOpt<std::string>("weights,w", "h5m remapping weights filename", &h5mfilename);
  opts.addOpt<std::string>("scrip,s", "Output SCRIP map filename", &scripfile);
  opts.addOpt<int>("dim,d", "Dimension of entities to use for partitioning", &dimension);
  opts.addOpt<void>("mesh,m", "Only convert the mesh and exclude the remap weight details", &noMap);
  opts.addOpt<void>("coords,c", "Write the center and vertex coordinates in lat/lon format", &writeXYCoords);
  
  opts.parseCommandLine(argc, argv);

  if (h5mfilename.empty() || scripfile.empty())
  {
    opts.printHelp();
    exit(1);
  }

  moab::Interface* mbCore = new ( std::nothrow ) moab::Core;

  if ( NULL == mbCore ) { return 1; }

  //Set the read options for parallel file loading
  const std::string partition_set_name = "PARALLEL_PARTITION";
  const std::string global_id_name = "GLOBAL_ID";

  //Load file
  rval = mbCore->load_mesh(h5mfilename.c_str());MB_CHK_ERR(rval);

  try {

    // Temporarily change rval reporting
    NcError error_temp(NcError::verbose_fatal);

    // Open an output file
    NcFile ncMap(scripfile.c_str(), NcFile::Replace, NULL, 0, NcFile::Offset64Bits);
    if (!ncMap.is_valid()) {
      _EXCEPTION1("Unable to open output map file \"%s\"",
        scripfile.c_str());
    }

    {
      // NetCDF-SCRIP Global Attributes
      std::map<std::string, std::string> mapAttributes;
      size_t lastindex = h5mfilename.find_last_of(".");
      std::stringstream sstr;
      sstr << h5mfilename.substr(0, lastindex) << ".meta";
      std::string metaFilename = sstr.str();
      ReadFileMetaData(metaFilename, mapAttributes);
      mapAttributes["Command"] = "Converted with MOAB:h5mtoscrip with --w=" + h5mfilename + " and --s=" + scripfile;

      // Add global attributes
      std::map<std::string, std::string>::const_iterator iterAttributes =
        mapAttributes.begin();
      for (; iterAttributes != mapAttributes.end(); iterAttributes++) {

        std::cout << iterAttributes->first << " -- " << iterAttributes->second << std::endl;
        ncMap.add_att(
          iterAttributes->first.c_str(),
          iterAttributes->second.c_str());
      }
      std::cout << "\n";
    }

    Tag globalIDTag, materialSetTag;
    globalIDTag = mbCore->globalId_tag();
    // materialSetTag = mbCore->material_tag();
    rval = mbCore->tag_get_handle( "MATERIAL_SET" , 1, MB_TYPE_INTEGER, materialSetTag, MB_TAG_SPARSE);MB_CHK_ERR(rval);

    // Get sets entities, by type
    moab::Range meshsets;
    rval = mbCore->get_entities_by_type_and_tag(0, MBENTITYSET, &globalIDTag, NULL, 1, meshsets, moab::Interface::UNION, true);MB_CHK_ERR(rval);

    moab::EntityHandle rootset = 0;
    ///////////////////////////////////////////////////////////////////////////
    // The metadata in H5M file contains the following data:
    //
    //   1. n_a: Total source entities: (number of elements in source mesh)
    //   2. n_b: Total target entities: (number of elements in target mesh)
    //   3. nv_a: Max edge size of elements in source mesh
    //   4. nv_b: Max edge size of elements in target mesh
    //   5. maxrows: Number of rows in remap weight matrix
    //   6. maxcols: Number of cols in remap weight matrix
    //   7. nnz: Number of total nnz in sparse remap weight matrix
    //   8. np_a: The order of the field description on the source mesh: >= 1 
    //   9. np_b: The order of the field description on the target mesh: >= 1
    //   10. method_a: The type of discretization for field on source mesh: [0 = FV, 1 = cGLL, 2 = dGLL]
    //   11. method_b: The type of discretization for field on target mesh: [0 = FV, 1 = cGLL, 2 = dGLL]
    //   12. conserved: Flag to specify whether the remap operator has conservation constraints: [0, 1]
    //   13. monotonicity: Flags to specify whether the remap operator has monotonicity constraints: [0, 1, 2]
    //
    ///////////////////////////////////////////////////////////////////////////    
    Tag smatMetadataTag;
    int smat_metadata_glb[13];
    rval = mbCore->tag_get_handle( "SMAT_DATA" , 13, MB_TYPE_INTEGER, smatMetadataTag, MB_TAG_SPARSE);MB_CHK_ERR(rval);
    rval = mbCore->tag_get_data(smatMetadataTag, &rootset, 1, smat_metadata_glb);MB_CHK_ERR(rval);
    // std::cout << "Number of mesh sets is " << meshsets.size() << std::endl;

#define DTYPE(a) { ((a == 0) ? "FV" : ((a == 1) ? "cGLL" : "dGLL")) }
    // Map dimensions
    int nA = smat_metadata_glb[0];
    int nB = smat_metadata_glb[1];
    int nVA = smat_metadata_glb[2];
    int nVB = smat_metadata_glb[3];
    int nDofB = smat_metadata_glb[4];
    int nDofA = smat_metadata_glb[5];
    int NNZ = smat_metadata_glb[6];
    int nOrdA = smat_metadata_glb[7];
    int nOrdB = smat_metadata_glb[8];
    int nBasA = smat_metadata_glb[9];
    std::string methodA = DTYPE(nBasA);
    int nBasB = smat_metadata_glb[10];
    std::string methodB = DTYPE(nBasB);
    int bConserved = smat_metadata_glb[11];
    int bMonotonicity = smat_metadata_glb[12];

    EntityHandle source_mesh=0, target_mesh=0, overlap_mesh=0;
    for (unsigned im=0; im < meshsets.size(); ++im) {
      moab::Range elems;
      rval = mbCore->get_entities_by_dimension(meshsets[im], 2, elems);MB_CHK_ERR(rval);
      if (elems.size() - nA == 0 && source_mesh == 0)
        source_mesh = meshsets[im];
      else if (elems.size() - nB == 0 && target_mesh == 0)
        target_mesh = meshsets[im];
      else if (overlap_mesh == 0)
        overlap_mesh = meshsets[im];
      else continue;
    }

    Tag srcIDTag, srcAreaTag, tgtIDTag, tgtAreaTag;
    rval = mbCore->tag_get_handle( "SourceGIDS" , srcIDTag);MB_CHK_ERR(rval);
    rval = mbCore->tag_get_handle( "SourceAreas" , srcAreaTag);MB_CHK_ERR(rval);
    rval = mbCore->tag_get_handle( "TargetGIDS" , tgtIDTag);MB_CHK_ERR(rval);
    rval = mbCore->tag_get_handle( "TargetAreas" , tgtAreaTag);MB_CHK_ERR(rval);
    Tag smatRowdataTag, smatColdataTag, smatValsdataTag;
    rval = mbCore->tag_get_handle( "SMAT_ROWS" , smatRowdataTag);MB_CHK_ERR(rval);
    rval = mbCore->tag_get_handle( "SMAT_COLS" , smatColdataTag);MB_CHK_ERR(rval);
    rval = mbCore->tag_get_handle( "SMAT_VALS" , smatValsdataTag);MB_CHK_ERR(rval);
    Tag srcCenterLon, srcCenterLat, tgtCenterLon, tgtCenterLat;
    rval = mbCore->tag_get_handle( "SourceCoordCenterLon" , srcCenterLon);MB_CHK_ERR(rval);
    rval = mbCore->tag_get_handle( "SourceCoordCenterLat" , srcCenterLat);MB_CHK_ERR(rval);
    rval = mbCore->tag_get_handle( "TargetCoordCenterLon" , tgtCenterLon);MB_CHK_ERR(rval);
    rval = mbCore->tag_get_handle( "TargetCoordCenterLat" , tgtCenterLat);MB_CHK_ERR(rval);
    Tag srcVertexLon, srcVertexLat, tgtVertexLon, tgtVertexLat;
    rval = mbCore->tag_get_handle( "SourceCoordVertexLon" , srcVertexLon);MB_CHK_ERR(rval);
    rval = mbCore->tag_get_handle( "SourceCoordVertexLat" , srcVertexLat);MB_CHK_ERR(rval);
    rval = mbCore->tag_get_handle( "TargetCoordVertexLon" , tgtVertexLon);MB_CHK_ERR(rval);
    rval = mbCore->tag_get_handle( "TargetCoordVertexLat" , tgtVertexLat);MB_CHK_ERR(rval);
    
    // Get sets entities, by type
    moab::Range sets;
    // rval = mbCore->get_entities_by_type(0, MBENTITYSET, sets);MB_CHK_ERR(rval);
    rval = mbCore->get_entities_by_type_and_tag(0, MBENTITYSET, &smatRowdataTag, NULL, 1, sets, moab::Interface::UNION, true);MB_CHK_ERR(rval);

    std::vector<int> src_gids, tgt_gids;
    std::vector<double> src_areas, tgt_areas;
    int srcID_size, tgtID_size, srcArea_size, tgtArea_size;
    rval = get_vartag_data(mbCore, srcIDTag, sets, srcID_size, src_gids);MB_CHK_SET_ERR(rval, "Getting source mesh IDs failed");
    rval = get_vartag_data(mbCore, tgtIDTag, sets, tgtID_size, tgt_gids);MB_CHK_SET_ERR(rval, "Getting target mesh IDs failed");
    rval = get_vartag_data(mbCore, srcAreaTag, sets, srcArea_size, src_areas);MB_CHK_SET_ERR(rval, "Getting source mesh areas failed");
    rval = get_vartag_data(mbCore, tgtAreaTag, sets, tgtArea_size, tgt_areas);MB_CHK_SET_ERR(rval, "Getting target mesh areas failed");

    assert(srcArea_size == srcID_size);
    assert(tgtArea_size == tgtID_size);

    std::vector<double> src_glob_areas(nDofA, 0.0), tgt_glob_areas(nDofB, 0.0);
    for (int i=0; i < srcArea_size; ++i) {
        // printf("%d/%d: %d = Found ID %d and area %5.6e\n", i, srcArea_size, nDofA, src_gids[i], src_areas[i]);
        assert(i < srcID_size);
        assert(src_gids[i] < nDofA);
        if (src_areas[i] > src_glob_areas[src_gids[i]])
          src_glob_areas[src_gids[i]] = src_areas[i];
    }
    for (int i=0; i < tgtArea_size; ++i) {
        // printf("%d/%d: %d = Found ID %d and area %5.6e\n", i, tgtArea_size, nDofB, tgt_gids[i], tgt_areas[i]);
        assert(i < tgtID_size);
        assert(tgt_gids[i] < nDofB);
        if (tgt_areas[i] > tgt_glob_areas[tgt_gids[i]])
          tgt_glob_areas[tgt_gids[i]] = tgt_areas[i];
    }

    // Write output dimensions entries
    int nSrcGridDims = 1;
    int nDstGridDims = 1;

    NcDim * dimSrcGridRank = ncMap.add_dim("src_grid_rank", nSrcGridDims);
    NcDim * dimDstGridRank = ncMap.add_dim("dst_grid_rank", nDstGridDims);

    NcVar * varSrcGridDims =
      ncMap.add_var("src_grid_dims", ncInt, dimSrcGridRank);
    NcVar * varDstGridDims =
      ncMap.add_var("dst_grid_dims", ncInt, dimDstGridRank);

    if (nA == nDofA) {
      varSrcGridDims->put(&nA, 1);
      varSrcGridDims->add_att("name0", "num_elem");
    }
    else {
      varSrcGridDims->put(&nDofA, 1);
      varSrcGridDims->add_att("name1", "num_dof");
    }

    if (nB == nDofB) {
      varDstGridDims->put(&nB, 1);
      varDstGridDims->add_att("name0", "num_elem");
    }
    else {
      varDstGridDims->put(&nDofB, 1);
      varDstGridDims->add_att("name1", "num_dof");
    }

    // Source and Target mesh resolutions
    NcDim * dimNA = ncMap.add_dim("n_a", nDofA);
    NcDim * dimNB = ncMap.add_dim("n_b", nDofB);

    // Source and Target verticecs per elements
    const int nva = (nA == nDofA ? nVA : 1);
    const int nvb = (nB == nDofB ? nVB : 1);
    NcDim * dimNVA = ncMap.add_dim("nv_a", nva);
    NcDim * dimNVB = ncMap.add_dim("nv_b", nvb);

    // Source and Target verticecs per elements
    // NcDim * dimNEA = ncMap.add_dim("ne_a", nA);
    // NcDim * dimNEB = ncMap.add_dim("ne_b", nB);

    if (writeXYCoords)
    {
      // Write coordinates
      NcVar * varYCA = ncMap.add_var("yc_a", ncDouble, dimNA/*dimNA*/);
      NcVar * varYCB = ncMap.add_var("yc_b", ncDouble, dimNB/*dimNB*/);

      NcVar * varXCA = ncMap.add_var("xc_a", ncDouble, dimNA/*dimNA*/);
      NcVar * varXCB = ncMap.add_var("xc_b", ncDouble, dimNB/*dimNB*/);

      NcVar * varYVA = ncMap.add_var("yv_a", ncDouble, dimNA/*dimNA*/, dimNVA);
      NcVar * varYVB = ncMap.add_var("yv_b", ncDouble, dimNB/*dimNB*/, dimNVB);

      NcVar * varXVA = ncMap.add_var("xv_a", ncDouble, dimNA/*dimNA*/, dimNVA);
      NcVar * varXVB = ncMap.add_var("xv_b", ncDouble, dimNB/*dimNB*/, dimNVB);

      varYCA->add_att("units", "degrees");
      varYCB->add_att("units", "degrees");

      varXCA->add_att("units", "degrees");
      varXCB->add_att("units", "degrees");

      varYVA->add_att("units", "degrees");
      varYVB->add_att("units", "degrees");

      varXVA->add_att("units", "degrees");
      varXVB->add_att("units", "degrees");

      std::vector<double> src_centerlat, src_centerlon;
      int srccenter_size;
      rval = get_vartag_data(mbCore, srcCenterLat, sets, srccenter_size, src_centerlat);MB_CHK_SET_ERR(rval, "Getting source mesh areas failed");
      rval = get_vartag_data(mbCore, srcCenterLon, sets, srccenter_size, src_centerlon);MB_CHK_SET_ERR(rval, "Getting target mesh areas failed");
      std::vector<double> src_glob_centerlat(nDofA, 0.0), src_glob_centerlon(nDofA, 0.0);

      for (int i=0; i < srccenter_size; ++i) {
          assert(i < srcID_size);
          assert(src_gids[i] < nDofA);

          src_glob_centerlat[src_gids[i]] = src_centerlat[i];
          src_glob_centerlon[src_gids[i]] = src_centerlon[i];
      }

      std::vector<double> tgt_centerlat, tgt_centerlon;
      int tgtcenter_size;
      rval = get_vartag_data(mbCore, tgtCenterLat, sets, tgtcenter_size, tgt_centerlat);MB_CHK_SET_ERR(rval, "Getting source mesh areas failed");
      rval = get_vartag_data(mbCore, tgtCenterLon, sets, tgtcenter_size, tgt_centerlon);MB_CHK_SET_ERR(rval, "Getting target mesh areas failed");
      std::vector<double> tgt_glob_centerlat(nDofB, 0.0), tgt_glob_centerlon(nDofB, 0.0);
      for (int i=0; i < tgtcenter_size; ++i) {
          assert(i < tgtID_size);
          assert(tgt_gids[i] < nDofB);

          tgt_glob_centerlat[tgt_gids[i]] = tgt_centerlat[i];
          tgt_glob_centerlon[tgt_gids[i]] = tgt_centerlon[i];
      }

      varYCA->put(&(src_glob_centerlat[0]), nDofA);
      varYCB->put(&(tgt_glob_centerlat[0]), nDofB);
      varXCA->put(&(src_glob_centerlon[0]), nDofA);
      varXCB->put(&(tgt_glob_centerlon[0]), nDofB);

      src_centerlat.clear();
      src_centerlon.clear();
      tgt_centerlat.clear();
      tgt_centerlon.clear();

      DataArray2D<double> src_glob_vertexlat(nDofA,nva), src_glob_vertexlon(nDofA,nva);
      if (nva > 1)
      {
        std::vector<double> src_vertexlat, src_vertexlon;
        int srcvertex_size1, srcvertex_size2;
        rval = get_vartag_data(mbCore, srcVertexLat, sets, srcvertex_size1, src_vertexlat);MB_CHK_SET_ERR(rval, "Getting source mesh areas failed");
        rval = get_vartag_data(mbCore, srcVertexLon, sets, srcvertex_size2, src_vertexlon);MB_CHK_SET_ERR(rval, "Getting target mesh areas failed");
        int offset = 0;
        printf("Source: %d, %d, %d, %d, %d, %d, %d, %d\n", nva, nDofA, srcvertex_size1, srcvertex_size2, src_gids.size()*nva, nDofA*nva, src_gids.size(), src_vertexlat.size());
        for (unsigned vIndex = 0; vIndex < src_gids.size(); ++vIndex)
        {
            for (int vNV = 0; vNV < nva; ++vNV)
            {
              // assert(offset < srcvertex_size1);
              if (offset < srcvertex_size1) src_glob_vertexlat[src_gids[vIndex]][vNV] = src_vertexlat[offset];
              else { printf("Offset = %d, and srcvertex_size1 = %d\n", offset, srcvertex_size1); }
              // assert(offset < srcvertex_size2);
              if (offset < srcvertex_size2) src_glob_vertexlon[src_gids[vIndex]][vNV] = src_vertexlon[offset];
              else { printf("Offset = %d, and srcvertex_size2 = %d\n", offset, srcvertex_size2); }
              offset++;
            }
        }
      }

      DataArray2D<double> tgt_glob_vertexlat(nDofB,nvb), tgt_glob_vertexlon(nDofB,nvb);
      if (nvb > 1)
      {
        std::vector<double> tgt_vertexlat, tgt_vertexlon;
        int tgtvertex_size;
        rval = get_vartag_data(mbCore, tgtVertexLat, sets, tgtvertex_size, tgt_vertexlat);MB_CHK_SET_ERR(rval, "Getting source mesh areas failed");
        rval = get_vartag_data(mbCore, tgtVertexLon, sets, tgtvertex_size, tgt_vertexlon);MB_CHK_SET_ERR(rval, "Getting target mesh areas failed");
        int offset = 0;
        for (unsigned vIndex = 0; vIndex < tgt_gids.size(); ++vIndex)
        {
          for (int vNV = 0; vNV < nvb; ++vNV)
          {
            assert(offset < tgtvertex_size);
            tgt_glob_vertexlat[tgt_gids[vIndex]][vNV] = tgt_vertexlat[offset];
            tgt_glob_vertexlon[tgt_gids[vIndex]][vNV] = tgt_vertexlon[offset];
            offset++;
          }
        }
      }

      varYVA->put(&(src_glob_vertexlat[0][0]), nDofA, nva);
      varYVB->put(&(tgt_glob_vertexlat[0][0]), nDofB, nvb);

      varXVA->put(&(src_glob_vertexlon[0][0]), nDofA, nva);
      varXVB->put(&(tgt_glob_vertexlon[0][0]), nDofB, nvb);
    }

    // Write areas
    NcVar * varAreaA = ncMap.add_var("area_a", ncDouble, dimNA);
    varAreaA->put(&(src_glob_areas[0]), nDofA);
    // varAreaA->add_att("units", "steradians");

    NcVar * varAreaB = ncMap.add_var("area_b", ncDouble, dimNB);
    varAreaB->put(&(tgt_glob_areas[0]), nDofB);
    // varAreaB->add_att("units", "steradians");

    std::vector<int> mat_rows, mat_cols;
    std::vector<double> mat_vals;
    int row_sizes, col_sizes, val_sizes;
    rval = get_vartag_data(mbCore, smatRowdataTag, sets, row_sizes, mat_rows);MB_CHK_SET_ERR(rval, "Getting matrix row data failed");
    assert(row_sizes == NNZ);
    rval = get_vartag_data(mbCore, smatColdataTag, sets, col_sizes, mat_cols);MB_CHK_SET_ERR(rval, "Getting matrix col data failed");
    assert(col_sizes == NNZ);
    rval = get_vartag_data(mbCore, smatValsdataTag, sets, val_sizes, mat_vals);MB_CHK_SET_ERR(rval, "Getting matrix values failed");
    assert(val_sizes == NNZ);

    // Let us form the matrix in-memory and consolidate shared DoF rows from shared-process contributions
    SparseMatrix<double> mapMatrix;

    for (int innz = 0; innz < NNZ; ++innz)
    {
#ifdef VERBOSE
      if (fabs(mapMatrix(mat_rows[innz], mat_cols[innz])) > 1e-12) 
      {
        printf("Adding to existing loc: (%d, %d) = %12.8f\n", mat_rows[innz], mat_cols[innz], mapMatrix(mat_rows[innz], mat_cols[innz]));
      }
#endif
      mapMatrix(mat_rows[innz], mat_cols[innz]) += mat_vals[innz];
    }

    // Write SparseMatrix entries
    DataArray1D<int> vecRow;
    DataArray1D<int> vecCol;
    DataArray1D<double> vecS;

    mapMatrix.GetEntries(vecRow, vecCol, vecS);

    int nS = vecS.GetRows();

    // Print more information about what we are converting:
    // Source elements/vertices/type (Discretization ?)
    // Target elements/vertices/type (Discretization ?)
    // Overlap elements/types
    // Rmeapping weights matrix: rows/cols/NNZ
    // Output the number of sets
    printf("Primary sets: %15zu\n", sets.size());
    printf("Original NNZ: %18d\n", NNZ);
    printf("Consolidated Total NNZ: %8d\n", nS);
    printf("Conservative weights ? %6d\n", (bConserved > 0));
    printf("Monotone weights ? %10d\n", (bMonotonicity > 0));

    printf("\n--------------------------------------------------------------\n");
    printf("%20s %21s %15s\n", "Description", "Source", "Target");
    printf("--------------------------------------------------------------\n");

    printf("%25s %15d %15d\n", "Number of elements:", nA, nB);
    printf("%25s %15d %15d\n", "Number of DoFs:", nDofA, nDofB);
    printf("%25s %15d %15d\n", "Maximum vertex/element:", nVA, nVB);
    printf("%25s %15s %15s\n", "Discretization type:", methodA.c_str(), methodB.c_str());
    printf("%25s %15d %15d\n", "Discretization order:", nOrdA, nOrdB);
    
    // Calculate and write fractional coverage arrays
    {
      DataArray1D<double> dFracA(nDofA);
      DataArray1D<double> dFracB(nDofB);
    
      for (int i = 0; i < nS; i++) {
        // std::cout << i << " - mat_vals = " << mat_vals[i] << " dFracA = " << mat_vals[i] / src_glob_areas[vecCol[i]] * tgt_glob_areas[vecRow[i]] << std::endl;
        dFracA[vecCol[i]] += vecS[i] / src_glob_areas[vecCol[i]] * tgt_glob_areas[vecRow[i]];
        dFracB[vecRow[i]] += vecS[i];
      }

      NcVar * varFracA = ncMap.add_var("frac_a", ncDouble, dimNA);
      varFracA->put(&(dFracA[0]), nDofA);
      varFracA->add_att("name", "fraction of target coverage of source dof");
      varFracA->add_att("units", "unitless");

      NcVar * varFracB = ncMap.add_var("frac_b", ncDouble, dimNB);
      varFracB->put(&(dFracB[0]), nDofB);
      varFracB->add_att("name", "fraction of source coverage of target dof");
      varFracB->add_att("units", "unitless");
    }

    // Write out data
    NcDim * dimNS = ncMap.add_dim("n_s", nS);

    NcVar * varRow = ncMap.add_var("row", ncInt, dimNS);
    varRow->add_att("name", "sparse matrix target dof index");
    varRow->add_att("first_index", "1");

    NcVar * varCol = ncMap.add_var("col", ncInt, dimNS);
    varCol->add_att("name", "sparse matrix source dof index");
    varCol->add_att("first_index", "1");

    NcVar * varS = ncMap.add_var("S", ncDouble, dimNS);
    varS->add_att("name", "sparse matrix coefficient");

    // Increment vecRow and vecCol: make it 1-based
    for (int i = 0; i < nS; i++) {
      vecRow[i]++;
      vecCol[i]++;
    }

    varRow->set_cur((long)0);
    varRow->put(&(vecRow[0]), nS);

    varCol->set_cur((long)0);
    varCol->put(&(vecCol[0]), nS);

    varS->set_cur((long)0);
    varS->put(&(vecS[0]), nS);

    ncMap.close();

    // rval = mbCore->write_file(scripfile.c_str());MB_CHK_ERR(rval);
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
