#include "moab/Core.hpp"
#include "moab/CartVect.hpp"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <netcdf.h>
#include <math.h>

using namespace moab;
using namespace std;

string test_file_name = string("dual_100kisland_2.h5m");
string mpas_file=  string("mpas.nc");


int main(int argc, char **argv)
{
  // Get MOAB instance
  Interface* mb = new (std::nothrow) Core;
  if (NULL == mb)
    return 1;

  // Need option handling here for input filename
  if (argc > 1) {
    // User has input a mesh file
    test_file_name = argv[1];
  }

  // Load the mesh from h5m  file
  ErrorCode rval = mb->load_mesh(test_file_name.c_str());MB_CHK_ERR(rval);
  cout<< " read " << test_file_name << endl;
  if (argc > 2) {
    mpas_file= argv[2];
  }

  // Try to open the file after gather mesh info succeeds
  int ncFile;
  int fail = nc_create(mpas_file.c_str(),  NC_CLOBBER, &ncFile);
  if (NC_NOERR != fail) {
    cout<< "can't open file\n";
    return 1;
  }
  cout << "open " << mpas_file << endl;

  Range polys;
  Range edges;
  Range verts;

  rval = mb->get_entities_by_dimension(0, 2, polys); MB_CHK_ERR(rval);
  rval = mb->get_adjacencies(polys,1,true,edges, Interface::UNION) ;  MB_CHK_ERR(rval);
  rval = mb->get_connectivity(polys, verts); MB_CHK_ERR(rval);

  cout << " cells:" << polys.size() << " edges: " << edges.size() << " vertices: " << verts.size() << "\n";



  const char * working_title ="MPAS mesh";

  if (NC_NOERR != nc_put_att_text(ncFile, NC_GLOBAL, "title", 9, working_title)) {
     MB_SET_ERR(MB_FAILURE, "WriteNCDF: failed to define title attribute");
   }
  const char * on_sphere="YES ";
  double sphere_radius=1.0;
  int nCells = (int)polys.size();
  int nEdges = (int)edges.size();
  int nVertices = (int)verts.size();
  fail = nc_put_att_text(ncFile, NC_GLOBAL, "on_a_sphere", 4, on_sphere);
  fail = nc_put_att_double(ncFile, NC_GLOBAL, "sphere_radius", NC_DOUBLE, 1, &sphere_radius);
  int np=nCells, n_scvt_iterations=10;
  double eps = 1.0e-9;
  fail = nc_put_att_int(ncFile, NC_GLOBAL, "np", NC_INT, 1, &np);
  fail = nc_put_att_int(ncFile, NC_GLOBAL, "n_scvt_iterations", NC_INT, 1, &n_scvt_iterations);
  fail = nc_put_att_double(ncFile, NC_GLOBAL, "eps", NC_DOUBLE, 1, &eps);

  fail = nc_put_att_text(ncFile, NC_GLOBAL, "Convergence",2,"L2");

  // Define dimensions
  int wrDimIDnCells, wrDimIDnEdges, wrDimIDnVertices, wrDimIDmaxEdges,
    wrDimIDmaxEdges2, wrDimIDTWO, wrDimIDvertexDegree, wrDimIDTime;

  if (nc_def_dim(ncFile, "nCells", nCells, &wrDimIDnCells) != NC_NOERR) {
    MB_SET_ERR(MB_FAILURE, "fail to create nCells");
  }
  if (nc_def_dim(ncFile, "nEdges", nEdges, &wrDimIDnEdges) != NC_NOERR) {
    MB_SET_ERR(MB_FAILURE, "fail to create nEdges dim");
  }
  if (nc_def_dim(ncFile, "nVertices", nVertices, &wrDimIDnVertices) != NC_NOERR) {
     MB_SET_ERR(MB_FAILURE, "fail to create nVertices dim");
   }
  // find max edges from polys
  int maxEdges = -1;
  for (Range::iterator it=polys.begin(); it!=polys.end(); it++)
  {
    const EntityHandle * conn=NULL;
    int nnodes=0;
    rval = mb->get_connectivity(*it, conn, nnodes); MB_CHK_ERR(rval);
    if (nnodes>maxEdges) maxEdges = nnodes;
  }
  if (nc_def_dim(ncFile, "maxEdges",maxEdges , &wrDimIDmaxEdges) != NC_NOERR) {
     MB_SET_ERR(MB_FAILURE, "fail to create maxEdges dim");
   }
  if (nc_def_dim(ncFile, "maxEdges2",  2*maxEdges, &wrDimIDmaxEdges2) != NC_NOERR) {
     MB_SET_ERR(MB_FAILURE, "fail to create maxEdges2 dim");
   }
  if (nc_def_dim(ncFile, "TWO",  2, &wrDimIDTWO) != NC_NOERR) {
       MB_SET_ERR(MB_FAILURE, "fail to create TWO dim");
     }
  if (nc_def_dim(ncFile, "vertexDegree",  3, &wrDimIDvertexDegree) != NC_NOERR) {
       MB_SET_ERR(MB_FAILURE, "fail to create vertexDegree dim");
     }
  if (nc_def_dim(ncFile, "Time",  NC_UNLIMITED, &wrDimIDTime) != NC_NOERR) {
         MB_SET_ERR(MB_FAILURE, "fail to create Time dim");
       }
/*
  !
  ! Define variables
  !
  dimlist( 1) = wrDimIDnCells
  nferr = nf_def_var(wr_ncid, 'latCell', NF_DOUBLE,  1, dimlist, wrVarIDlatCell) */

  int dimlist[3];
  int wrVarIDlatCell;
  dimlist[0] = wrDimIDnCells;
  if (NC_NOERR != nc_def_var(ncFile, "latCell", NC_DOUBLE, 1, dimlist, &wrVarIDlatCell)) {
    MB_SET_ERR(MB_FAILURE, "fail to define latCell var");
  }

  /*
  dimlist( 1) = wrDimIDnCells
  nferr = nf_def_var(wr_ncid, 'lonCell', NF_DOUBLE,  1, dimlist, wrVarIDlonCell)*/
  int wrVarIDlonCell;
  dimlist[0] = wrDimIDnCells;
  if (NC_NOERR != nc_def_var(ncFile, "lonCell", NC_DOUBLE, 1, dimlist, &wrVarIDlonCell)) {
    MB_SET_ERR(MB_FAILURE, "fail to define lonCell var");
  }

  int wrVarIDmeshDensity;
  dimlist[0] = wrDimIDnCells;
  if (NC_NOERR != nc_def_var(ncFile, "meshDensity", NC_DOUBLE, 1, dimlist, &wrVarIDmeshDensity)) {
    MB_SET_ERR(MB_FAILURE, "fail to define wrVarIDmeshDensity var");
  }

  int wrVarIDxCell, wrVarIDyCell, wrVarIDzCell;
  if (NC_NOERR != nc_def_var(ncFile, "xCell", NC_DOUBLE, 1, dimlist, &wrVarIDxCell)) {
      MB_SET_ERR(MB_FAILURE, "fail to define wrVarIDmeshDensity var");
    }
  if (NC_NOERR != nc_def_var(ncFile, "yCell", NC_DOUBLE, 1, dimlist, &wrVarIDyCell)) {
      MB_SET_ERR(MB_FAILURE, "fail to define wrVarIDmeshDensity var");
    }
  if (NC_NOERR != nc_def_var(ncFile, "zCell", NC_DOUBLE, 1, dimlist, &wrVarIDzCell)) {
      MB_SET_ERR(MB_FAILURE, "fail to define wrVarIDmeshDensity var");
    }


  int wrVarIDindexToCellID;
  if (NC_NOERR != nc_def_var(ncFile, "indexToCellID", NC_INT, 1, dimlist, &wrVarIDindexToCellID)) {
    MB_SET_ERR(MB_FAILURE, "fail to define indexToCellID var");
  }

  // edges:
  dimlist[0] = wrDimIDnEdges;
  int wrVarIDlatEdge, wrVarIDlonEdge, wrVarIDxEdge, wrVarIDyEdge, wrVarIDzEdge, wrVarIDindexToEdgeID;
  if (NC_NOERR != nc_def_var(ncFile, "latEdge", NC_DOUBLE, 1, dimlist, &wrVarIDlatEdge)) {
    MB_SET_ERR(MB_FAILURE, "fail to define latEdge var");
  }
  if (NC_NOERR != nc_def_var(ncFile, "lonEdge", NC_DOUBLE, 1, dimlist, &wrVarIDlonEdge)) {
    MB_SET_ERR(MB_FAILURE, "fail to define lonEdge var");
  }
  if (NC_NOERR != nc_def_var(ncFile, "xEdge", NC_DOUBLE, 1, dimlist, &wrVarIDxEdge)) {
    MB_SET_ERR(MB_FAILURE, "fail to define xEdge var");
  }
  if (NC_NOERR != nc_def_var(ncFile, "yEdge", NC_DOUBLE, 1, dimlist, &wrVarIDyEdge)) {
    MB_SET_ERR(MB_FAILURE, "fail to define yEdge var");
  }
  if (NC_NOERR != nc_def_var(ncFile, "zEdge", NC_DOUBLE, 1, dimlist, &wrVarIDzEdge)) {
    MB_SET_ERR(MB_FAILURE, "fail to define zEdge var");
  }
  if (NC_NOERR != nc_def_var(ncFile, "indexToEdgeID", NC_INT, 1, dimlist, &wrVarIDindexToEdgeID)) {
    MB_SET_ERR(MB_FAILURE, "fail to define indexToEdgeID var");
  }
  dimlist[0]= wrDimIDTWO;
  dimlist[1] = wrDimIDnEdges;
  int wrVarIDcellsOnEdge;
  if (NC_NOERR != nc_def_var(ncFile, "cellsOnEdge", NC_INT, 2, dimlist, &wrVarIDcellsOnEdge)) {
    MB_SET_ERR(MB_FAILURE, "fail to define cellsOnEdge var");
  }

  int wrVarIDlatVertex, wrVarIDlonVertex, wrVarIDxVertex, wrVarIDyVertex, wrVarIDzVertex;
  int wrVarIDindexToVertexID;
  dimlist[0] = wrDimIDnVertices;
  if (NC_NOERR != nc_def_var(ncFile, "latVertex", NC_DOUBLE,  1, dimlist, &wrVarIDlatVertex) )
  {
    MB_SET_ERR(MB_FAILURE, "fail to define latVertex var");
  }
  if (NC_NOERR != nc_def_var(ncFile, "lonVertex", NC_DOUBLE,  1, dimlist, &wrVarIDlonVertex) )
  {
    MB_SET_ERR(MB_FAILURE, "fail to define lonVertex var");
  }
  if (NC_NOERR != nc_def_var(ncFile, "xVertex", NC_DOUBLE,  1, dimlist, &wrVarIDxVertex))
  {
    MB_SET_ERR(MB_FAILURE, "fail to define xVertex var");
  }
  if (NC_NOERR != nc_def_var(ncFile, "yVertex", NC_DOUBLE,  1, dimlist, &wrVarIDyVertex))
  {
    MB_SET_ERR(MB_FAILURE, "fail to define yVertex var");
  }
  if (NC_NOERR != nc_def_var(ncFile, "zVertex", NC_DOUBLE,  1, dimlist, &wrVarIDzVertex))
  {
    MB_SET_ERR(MB_FAILURE, "fail to define zVertex var");
  }
  if (NC_NOERR != nc_def_var(ncFile, "indexToVertexID", NC_INT,  1, dimlist, &wrVarIDindexToVertexID))
  {
    MB_SET_ERR(MB_FAILURE, "fail to define indexToVertexID var");
  }

  /* */


  /*
  dimlist( 1) = wrDimIDnCells
  nferr = nf_def_var(wr_ncid, 'nEdgesOnCell', NF_INT,  1, dimlist, wrVarIDnEdgesOnCell)
  dimlist( 1) = wrDimIDnEdges
  nferr = nf_def_var(wr_ncid, 'nEdgesOnEdge', NF_INT,  1, dimlist, wrVarIDnEdgesOnEdge)
  dimlist( 1) = wrDimIDmaxEdges
  dimlist( 2) = wrDimIDnCells
  nferr = nf_def_var(wr_ncid, 'edgesOnCell', NF_INT,  2, dimlist, wrVarIDedgesOnCell)
  dimlist( 1) = wrDimIDmaxEdges2
  dimlist( 2) = wrDimIDnEdges
  nferr = nf_def_var(wr_ncid, 'edgesOnEdge', NF_INT,  2, dimlist, wrVarIDedgesOnEdge)
  dimlist( 1) = wrDimIDmaxEdges2
  dimlist( 2) = wrDimIDnEdges
  nferr = nf_def_var(wr_ncid, 'weightsOnEdge', NF_DOUBLE,  2, dimlist, wrVarIDweightsOnEdge)
  dimlist( 1) = wrDimIDnEdges
  nferr = nf_def_var(wr_ncid, 'dvEdge', NF_DOUBLE,  1, dimlist, wrVarIDdvEdge)
  dimlist( 1) = wrDimIDnEdges
  nferr = nf_def_var(wr_ncid, 'dv1Edge', NF_DOUBLE,  1, dimlist, wrVarIDdv1Edge)
  dimlist( 1) = wrDimIDnEdges
  nferr = nf_def_var(wr_ncid, 'dv2Edge', NF_DOUBLE,  1, dimlist, wrVarIDdv2Edge)
  dimlist( 1) = wrDimIDnEdges
  nferr = nf_def_var(wr_ncid, 'dcEdge', NF_DOUBLE,  1, dimlist, wrVarIDdcEdge)
  dimlist( 1) = wrDimIDnEdges
  nferr = nf_def_var(wr_ncid, 'angleEdge', NF_DOUBLE,  1, dimlist, wrVarIDangleEdge)
  dimlist( 1) = wrDimIDnCells
  nferr = nf_def_var(wr_ncid, 'areaCell', NF_DOUBLE,  1, dimlist, wrVarIDareaCell)
  dimlist( 1) = wrDimIDnVertices
  nferr = nf_def_var(wr_ncid, 'areaTriangle', NF_DOUBLE,  1, dimlist, wrVarIDareaTriangle)
  dimlist( 1) = wrDimIDmaxEdges
  dimlist( 2) = wrDimIDnCells
  nferr = nf_def_var(wr_ncid, 'cellsOnCell', NF_INT,  2, dimlist, wrVarIDcellsOnCell)
  dimlist( 1) = wrDimIDmaxEdges
  dimlist( 2) = wrDimIDnCells
  nferr = nf_def_var(wr_ncid, 'verticesOnCell', NF_INT,  2, dimlist, wrVarIDverticesOnCell)
  dimlist( 1) = wrDimIDTWO
  dimlist( 2) = wrDimIDnEdges
  nferr = nf_def_var(wr_ncid, 'verticesOnEdge', NF_INT,  2, dimlist, wrVarIDverticesOnEdge)
  dimlist( 1) = wrDimIDvertexDegree
  dimlist( 2) = wrDimIDnVertices
  nferr = nf_def_var(wr_ncid, 'edgesOnVertex', NF_INT,  2, dimlist, wrVarIDedgesOnVertex)
  dimlist( 1) = wrDimIDvertexDegree
  dimlist( 2) = wrDimIDnVertices
  nferr = nf_def_var(wr_ncid, 'cellsOnVertex', NF_INT,  2, dimlist, wrVarIDcellsOnVertex)
  dimlist( 1) = wrDimIDvertexDegree
  dimlist( 2) = wrDimIDnVertices
  nferr = nf_def_var(wr_ncid, 'kiteAreasOnVertex', NF_DOUBLE,  2, dimlist, wrVarIDkiteAreasOnVertex)
  dimlist( 1) = wrDimIDnEdges
  nferr = nf_def_var(wr_ncid, 'fEdge', NF_DOUBLE,  1, dimlist, wrVarIDfEdge)
  dimlist( 1) = wrDimIDnVertices
  nferr = nf_def_var(wr_ncid, 'fVertex', NF_DOUBLE,  1, dimlist, wrVarIDfVertex)

  nferr = nf_enddef(wr_ncid)
*/

  // Take it out of define mode
  if (NC_NOERR != nc_enddef(ncFile)) {
    MB_SET_ERR(MB_FAILURE, " Trouble leaving define mode");
  }

  // write data
  /*start1(1) = 1
    count1( 1) = wrLocalnCells
    nferr = nf_put_vara_double(wr_ncid, wrVarIDlatCell, start1, count1, latCell)

    start1(1) = 1
    count1( 1) = wrLocalnCells
    nferr = nf_put_vara_double(wr_ncid, wrVarIDlonCell, start1, count1, lonCell)
    */
  // get the coordinates of centers of cells
  std::vector<double> centers(nCells*3);
  // TODO: retrieve actual generators
  rval = mb->get_coords(polys, &centers[0]);  MB_CHK_ERR(rval);
  std::vector<double> wkarr1(nCells), wkarr2(nCells);
  // convert coords to lat lon
  for (int i=0; i<nCells; i++)
  {
    /*
     *   lat(i) = (pii/2.0 - acos(z(i)))
         lon(i) = atan2(y(i),x(i))
     */
    wkarr1[i] = M_PI/2- acos(centers[i*3+2]);
    wkarr2[i] = atan2(centers[i*3+1], centers[i*3]); //
  }
  size_t start, count;
  start = 0;
  count = nCells;
  fail = nc_put_vara_double(ncFile, wrVarIDlatCell, &start, &count, &wkarr1[0]);
  if (NC_NOERR != fail) {
    MB_SET_ERR(MB_FAILURE, "Failed writing latCell");
  }
  fail = nc_put_vara_double(ncFile, wrVarIDlonCell, &start, &count, &wkarr2[0]);
  if (NC_NOERR != fail) {
    MB_SET_ERR(MB_FAILURE, "Failed writing lonCell");
  }
  /*start1(1) = 1
       count1( 1) = wrLocalnCells
       nferr = nf_put_vara_double(wr_ncid, wrVarIDmeshDensity, start1, count1, meshDensity)*/
  for (int i=0; i<nCells; i++)
    wkarr1[i] = 1; // dummy mesh density TODO
  fail = nc_put_vara_double(ncFile, wrVarIDmeshDensity, &start, &count, &wkarr1[0]);
  if (NC_NOERR != fail) {
    MB_SET_ERR(MB_FAILURE, "Failed writing mesh density");
  }
  /*
   * nferr = nf_put_vara_double(wr_ncid, wrVarIDxCell, start1, count1, xCell)

      start1(1) = 1
      count1( 1) = wrLocalnCells
      nferr = nf_put_vara_double(wr_ncid, wrVarIDyCell, start1, count1, yCell)

      start1(1) = 1
      count1( 1) = wrLocalnCells
      nferr = nf_put_vara_double(wr_ncid, wrVarIDzCell, start1, count1, zCell)
   */
   for (int i=0; i<nCells; i++)
      wkarr1[i] = centers[3*i]; //
  fail = nc_put_vara_double(ncFile, wrVarIDxCell, &start, &count, &wkarr1[0]);
  if (NC_NOERR != fail) {
    MB_SET_ERR(MB_FAILURE, "Failed writing xCell");
  }
  for (int i = 0; i < nCells; i++)
    wkarr1[i] = centers[3 * i + 1]; //
  fail = nc_put_vara_double(ncFile, wrVarIDyCell, &start, &count, &wkarr1[0]);
  if (NC_NOERR != fail) {
    MB_SET_ERR(MB_FAILURE, "Failed writing yCell");
  }

  for (int i = 0; i < nCells; i++)
    wkarr1[i] = centers[3 * i + 2]; //
  fail = nc_put_vara_double(ncFile, wrVarIDzCell, &start, &count, &wkarr1[0]);
  if (NC_NOERR != fail) {
    MB_SET_ERR(MB_FAILURE, "Failed writing zCell");
  }
/*

      start1(1) = 1
      count1( 1) = wrLocalnCells
      nferr = nf_put_vara_int(wr_ncid, wrVarIDindexToCellID, start1, count1, indexToCellID)
 */
  std::vector<int> wkint(nCells);
  for (int i=0; i<nCells; i++)
    wkint[i] = i+1;
  fail = nc_put_vara_int(ncFile, wrVarIDindexToCellID, &start, &count, &wkint[0]);
  if (NC_NOERR != fail) {
    MB_SET_ERR(MB_FAILURE, "Failed writing indexToCellID");
  }

  /*
      start1(1) = 1
      count1( 1) = wrLocalnEdges
      nferr = nf_put_vara_double(wr_ncid, wrVarIDlatEdge, start1, count1, latEdge)

      start1(1) = 1
      count1( 1) = wrLocalnEdges
      nferr = nf_put_vara_double(wr_ncid, wrVarIDlonEdge, start1, count1, lonEdge)

      start1(1) = 1
      count1( 1) = wrLocalnEdges
      nferr = nf_put_vara_double(wr_ncid, wrVarIDxEdge, start1, count1, xEdge)

      start1(1) = 1
      count1( 1) = wrLocalnEdges
      nferr = nf_put_vara_double(wr_ncid, wrVarIDyEdge, start1, count1, yEdge)

      start1(1) = 1
      count1( 1) = wrLocalnEdges
      nferr = nf_put_vara_double(wr_ncid, wrVarIDzEdge, start1, count1, zEdge)

      start1(1) = 1
      count1( 1) = wrLocalnEdges
      nferr = nf_put_vara_int(wr_ncid, wrVarIDindexToEdgeID, start1, count1, indexToEdgeID) */
  /*start2(2) = 1
        count2( 1) = 2
        count2( 2) = wrLocalnEdges
        nferr = nf_put_vara_int(wr_ncid, wrVarIDcellsOnEdge, start2, count2, cellsOnEdge)
    start1(1) = 1
      count1( 1) = wrLocalnEdges
      nferr = nf_put_vara_int(wr_ncid, wrVarIDnEdgesOnEdge, start1, count1, nEdgesOnEdge)
    start2(2) = 1
        count2( 1) = 2
        count2( 2) = wrLocalnEdges
        nferr = nf_put_vara_int(wr_ncid, wrVarIDverticesOnEdge, start2, count2, verticesOnEdge)

      start2(2) = 1
      count2( 1) = 2*wrLocalmaxEdges
      count2( 2) = wrLocalnEdges
      nferr = nf_put_vara_int(wr_ncid, wrVarIDedgesOnEdge, start2, count2, edgesOnEdge)

      start2(2) = 1
      count2( 1) = 2*wrLocalmaxEdges
      count2( 2) = wrLocalnEdges
      nferr = nf_put_vara_double(wr_ncid, wrVarIDweightsOnEdge, start2, count2, weightsOnEdge)

      start1(1) = 1
      count1( 1) = wrLocalnEdges
      nferr = nf_put_vara_double(wr_ncid, wrVarIDdvEdge, start1, count1, dvEdge)

      start1(1) = 1
      count1( 1) = wrLocalnEdges
      nferr = nf_put_vara_double(wr_ncid, wrVarIDdv1Edge, start1, count1, dv1Edge)

      start1(1) = 1
      count1( 1) = wrLocalnEdges
      nferr = nf_put_vara_double(wr_ncid, wrVarIDdv2Edge, start1, count1, dv2Edge)

      start1(1) = 1
      count1( 1) = wrLocalnEdges
      nferr = nf_put_vara_double(wr_ncid, wrVarIDdcEdge, start1, count1, dcEdge)

      start1(1) = 1
      count1( 1) = wrLocalnEdges
      nferr = nf_put_vara_double(wr_ncid, wrVarIDangleEdge, start1, count1, angleEdge)


 */
  /* edges */
  std::vector<double> latEdge(nEdges), lonEdge(nEdges);
  std::vector<double> xEdge(nEdges), yEdge(nEdges), zEdge(nEdges);
  std::vector<int>  indexToEdgeID(nEdges);
  std::vector<int> cellsOnEdge(2*nEdges); // 2d array;
  std::vector<int> verticesOnEdge(2*nEdges);
  std::vector<int> nEdgesOnEdge(nEdges);
  std::vector<int> edgesOnEdge(2*maxEdges*nEdges); // will have to fill up to maxEdges
  for (int i=0; i<nEdges; i++)
  {
    EntityHandle edge=edges[i];
    const EntityHandle * conn2=NULL;
    int nnodes=2;
    rval = mb->get_connectivity(edge, conn2, nnodes);  MB_CHK_ERR(rval);
    if (nnodes!=2)
      return 1;
    verticesOnEdge[2*i]=(int)conn2[0];
    verticesOnEdge[2*i+1]=(int)conn2[1];
    CartVect v[2];
    rval = mb->get_coords(conn2, 2, &(v[0][0])); MB_CHK_ERR(rval);
    CartVect mid=0.5*(v[0]+v[1]);
    mid.normalize();
    xEdge[i]=mid[0]; yEdge[i]=mid[1]; zEdge[i] = mid[2];
    latEdge[i] = M_PI/2- acos(mid[2]);
    latEdge[i] = atan2(mid[1], mid[0]); //
    indexToEdgeID[i] = i+1;
    Range adjCells;
    rval = mb->get_adjacencies(&edge, 1, 2, false, adjCells); MB_CHK_ERR(rval);
    // it should be at least one
    if (adjCells.size()<1)
      return 1; // error
    EntityHandle cell0=adjCells[0], cell1=0;
    if (adjCells.size()==2)
      cell1=adjCells[1];
    cellsOnEdge[i*2] = polys.index(cell0)+1;
    cellsOnEdge[i*2+1] = polys.index(cell1)+1; // could be 0
  }
  start = 0;
  count = nEdges;
  fail = nc_put_vara_double(ncFile, wrVarIDlatEdge, &start, &count, &latEdge[0]);
  if (NC_NOERR != fail) {
    MB_SET_ERR(MB_FAILURE, "Failed writing latEdge");
  }

  fail = nc_put_vara_double(ncFile, wrVarIDlonEdge, &start, &count, &lonEdge[0]);
  if (NC_NOERR != fail) {
    MB_SET_ERR(MB_FAILURE, "Failed writing lonEdge");
  }

  fail = nc_put_vara_double(ncFile, wrVarIDxEdge, &start, &count, &xEdge[0]);
  if (NC_NOERR != fail) {
    MB_SET_ERR(MB_FAILURE, "Failed writing xEdge");
  }

  fail = nc_put_vara_double(ncFile, wrVarIDyEdge, &start, &count, &yEdge[0]);
  if (NC_NOERR != fail) {
    MB_SET_ERR(MB_FAILURE, "Failed writing yEdge");
  }

  fail = nc_put_vara_double(ncFile, wrVarIDzEdge,  &start, &count, &zEdge[0]);
  if (NC_NOERR != fail) {
    MB_SET_ERR(MB_FAILURE, "Failed writing zEdge");
  }

  fail = nc_put_vara_int(ncFile, wrVarIDindexToEdgeID, &start, &count, &indexToEdgeID[0]);
  if (NC_NOERR != fail) {
    MB_SET_ERR(MB_FAILURE, "Failed writing indexToEdgeID");
  }
  size_t start2[2], count2[2];
  start2[0]=start2[1]=0;
  count2[0]=2; count2[1]=nEdges;

  fail = nc_put_vara_int(ncFile, wrVarIDcellsOnEdge, start2, count2, &cellsOnEdge[0]);
  if (NC_NOERR != fail) {
    MB_SET_ERR(MB_FAILURE, "Failed writing indexToEdgeID");
  }

  std::vector<double> coordsv(3*verts.size());
  rval = mb->get_coords(verts, &coordsv[0]);  MB_CHK_ERR(rval);
  std::vector<double> wkx(nVertices);

  count = nVertices;
  std::vector<double> latVertex(nVertices), lonVertex(nVertices);

  for (int i=0; i<nVertices; i++)
  {
    latVertex[i] = M_PI/2- acos(coordsv[i*3+2]);
    lonVertex[i] = atan2(coordsv[i*3+1], coordsv[i*3]); //
  }
  fail = nc_put_vara_double(ncFile, wrVarIDlatVertex, &start, &count, &latVertex[0]);
  if (NC_NOERR != fail) {
    MB_SET_ERR(MB_FAILURE, "Failed writing latVertex");
  }

  fail = nc_put_vara_double(ncFile,  wrVarIDlonVertex, &start, &count, &lonVertex[0]);

  for (int i=0; i<nVertices; i++)
    wkx[i]=coordsv[3*i];

  fail = nc_put_vara_double(ncFile, wrVarIDxVertex, &start, &count, &wkx[0]);
  if (NC_NOERR != fail) {
    MB_SET_ERR(MB_FAILURE, "Failed writing xVertex");
  }

  for (int i=0; i<nVertices; i++)
    wkx[i]=coordsv[3*i+1];

  fail = nc_put_vara_double(ncFile, wrVarIDyVertex, &start, &count, &wkx[0]);
  if (NC_NOERR != fail) {
    MB_SET_ERR(MB_FAILURE, "Failed writing yVertex");
  }
  for (int i=0; i<nVertices; i++)
    wkx[i]=coordsv[3*i+2];

  fail = nc_put_vara_double(ncFile, wrVarIDzVertex, &start, &count, &wkx[0]);
  if (NC_NOERR != fail) {
    MB_SET_ERR(MB_FAILURE, "Failed writing zVertex");
  }

  wkint.resize(nVertices);
  for (int i=0;i<nVertices; i++)
    wkint[i] = i+1;

  fail = nc_put_vara_int(ncFile, wrVarIDindexToVertexID, &start, &count, &wkint[0]);
  if (NC_NOERR != fail) {
    MB_SET_ERR(MB_FAILURE, "Failed writing indexToVertexID");
  }
  /*

      start1(1) = 1
      count1( 1) = wrLocalnCells
      nferr = nf_put_vara_int(wr_ncid, wrVarIDnEdgesOnCell, start1, count1, nEdgesOnCell)

      start2(2) = 1
      count2( 1) = wrLocalmaxEdges
      count2( 2) = wrLocalnCells
      nferr = nf_put_vara_int(wr_ncid, wrVarIDedgesOnCell, start2, count2, edgesOnCell)

      start1(1) = 1
      count1( 1) = wrLocalnCells
      nferr = nf_put_vara_double(wr_ncid, wrVarIDareaCell, start1, count1, areaCell)

      start1(1) = 1
      count1( 1) = wrLocalnVertices
      nferr = nf_put_vara_double(wr_ncid, wrVarIDareaTriangle, start1, count1, areaTriangle)

      start2(2) = 1
      count2( 1) = wrLocalmaxEdges
      count2( 2) = wrLocalnCells
      nferr = nf_put_vara_int(wr_ncid, wrVarIDcellsOnCell, start2, count2, cellsOnCell)

      start2(2) = 1
      count2( 1) = wrLocalmaxEdges
      count2( 2) = wrLocalnCells
      nferr = nf_put_vara_int(wr_ncid, wrVarIDverticesOnCell, start2, count2, verticesOnCell)*/



  /*

      start2(2) = 1
      count2( 1) = 3
      count2( 2) = wrLocalnVertices
      nferr = nf_put_vara_int(wr_ncid, wrVarIDedgesOnVertex, start2, count2, edgesOnVertex)

      start2(2) = 1
      count2( 1) = 3
      count2( 2) = wrLocalnVertices
      nferr = nf_put_vara_int(wr_ncid, wrVarIDcellsOnVertex, start2, count2, cellsOnVertex)

      start2(2) = 1
      count2( 1) = 3
      count2( 2) = wrLocalnVertices
      nferr = nf_put_vara_double(wr_ncid, wrVarIDkiteAreasOnVertex, start2, count2, kiteAreasOnVertex)

      start1(1) = 1
      count1( 1) = wrLocalnEdges
      nferr = nf_put_vara_double(wr_ncid, wrVarIDfEdge, start1, count1, fEdge)

      start1(1) = 1
      count1( 1) = wrLocalnVertices
      nferr = nf_put_vara_double(wr_ncid, wrVarIDfVertex, start1, count1, fVertex)
 */

  fail = nc_close(ncFile);
  if (NC_NOERR != fail) {
    MB_SET_ERR(MB_FAILURE, "Trouble closing file");
  }
  cout << " close " << mpas_file << endl;



}
