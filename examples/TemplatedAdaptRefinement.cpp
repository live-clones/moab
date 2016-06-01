/*This unit test is for the uniform refinement capability based on AHF datastructures*/
#include <iostream>
#include <math.h>
#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"


using namespace moab;

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

#define MAX_VERTS 16
#define MAX_CHILDS 10
struct refTemplates {
  int num_verts;
  int num_ents;
  double vertex_natcoords[MAX_VERTS][3];
  int ents_conn[MAX_CHILDS][4];
};

static const refTemplates quadTemp = {
  4, 5, {{-1/3.0, -1/3.0,0},{1/3.0, -1/3.0,0},{1/3.0,1/3.0,0},{-1/3.0, 1/3.0,0}}, {{0,1,5,4},{1,2,6,5},{3,7,6,2},{0,4,7,3},{4,5,6,7}}
};

bool Errfunc(double point[3], double tol)
{
  double center[3] = {0,0,0};
  double val = pow(point[0]-center[0],2)+pow(point[1]-center[1],2)+pow(point[2]-center[2],2);
  if (exp(-val) > tol)
    return true;
  else
    return false;
}

void compute_coords(double *corner_coords, int num_coords, double xi, double eta, double verts[3])
{
  double N[4];
   N[0] = (1-xi)*(1-eta)/4; N[1] = (1+xi)*(1-eta)/4; N[2] = (1+xi)*(1+eta)/4, N[3] = (1-xi)*(1+eta)/4;

   double x=0, y=0, z=0;
   for (int j=0; j<num_coords; j++)
     {
       x += N[j]*corner_coords[3*j];
       y += N[j]*corner_coords[3*j+1];
       z += N[j]*corner_coords[3*j+2];
     }
   verts[0] = x; verts[1] = y; verts[2] = z;

}

int main(int argc, char *argv[])
{

  Core mb;
  Interface* mbImpl = &mb;
  ErrorCode error;


  if (argc==1)
    {
      std::cerr << "Usage: " << argv[0] << " [filename]" << std::endl;
      return 1;
    }

  //Load Mesh
  const char *filename = argv[1];
  EntityHandle inmesh;
  error = mbImpl->create_meshset(MESHSET_SET, inmesh);MB_CHK_ERR(error);
  error = mbImpl->load_file(filename, &inmesh);MB_CHK_ERR(error);
  Range allverts, allquads;
  error = mbImpl->get_entities_by_dimension(inmesh, 0, allverts);MB_CHK_ERR(error);
  error = mbImpl->get_entities_by_dimension(inmesh, 2, allquads);MB_CHK_ERR(error);

  //Create refinement set from error indicators
  EntityHandle refset;
  error = mbImpl->create_meshset(MESHSET_SET, refset);MB_CHK_ERR(error);
  double centroid[3];
  double tol = 0.5;
  for (Range::iterator it = allquads.begin(); it != allquads.end(); it++)
    {
      error = mbImpl->get_coords(&(*it), 1, &centroid[0]);MB_CHK_ERR(error);
      bool select = Errfunc(centroid, tol);
      if (select)
        {
          error = mbImpl->add_entities(refset, &(*it), 1);MB_CHK_ERR(error);
        }
    }

  //Create storage for the refined entities
  //ReadUtilIface *iface;
  //error = mbImpl->query_interface(iface);MB_CHK_ERR(error);



  //Refine the entities in the refinement set
  EntityHandle finemesh;
  error = mbImpl->create_meshset(MESHSET_SET, finemesh);MB_CHK_ERR(error);

  Range coarse_quads;
  error = mbImpl->get_entities_by_dimension(refset, 2, coarse_quads);MB_CHK_ERR(error);

  Range remquads = subtract(allquads, coarse_quads);
  error = mbImpl->add_entities(finemesh, remquads);MB_CHK_ERR(error);

  for (Range::iterator it = coarse_quads.begin(); it != coarse_quads.end(); it++)
    {
      int nconn = 0;
      const EntityHandle *conn;
      error = mbImpl->get_connectivity(*it, conn, nconn);MB_CHK_ERR(error);

      double *coords = new double[3*nconn];
      error = mbImpl->get_coords(conn, nconn, coords);MB_CHK_ERR(error);

      //Create new vertices
      EntityHandle nwvert;  double x[3];
      EntityHandle *vbuffer = new EntityHandle[nconn+quadTemp.num_verts];

      for (int i=0; i<nconn; i++)
        vbuffer[i] = conn[i];

      for (int i = 0; i < quadTemp.num_verts; i++)
        {
          double xi = quadTemp.vertex_natcoords[i][0];
          double eta = quadTemp.vertex_natcoords[i][1];
          compute_coords(coords, nconn, xi, eta, x);

          error = mbImpl->create_vertex(x, nwvert);MB_CHK_ERR(error);
          vbuffer[i+nconn] = nwvert;
        }

      //Create new quads
      EntityHandle nwquad;
      EntityHandle *nwconn = new EntityHandle[nconn];
      for (int i=0; i<quadTemp.num_ents; i++)
        {
          for (int k=0; k<nconn; k++)
            {
              int idx = quadTemp.ents_conn[i][k];
              nwconn[k] = vbuffer[idx];
            }
          error = mbImpl->create_element(MBQUAD, nwconn, nconn, nwquad);MB_CHK_ERR(error);
          error = mbImpl->add_entities(finemesh, &nwquad, 1);MB_CHK_ERR(error);
        }

      delete [] vbuffer;
      delete [] coords;
      delete [] nwconn;
    }


  //Write out
  std::stringstream file;
  file <<"refined_mesh.vtk";
  error = mbImpl->write_file(file.str().c_str(), 0, NULL, &finemesh, 1);MB_CHK_ERR(error);

  return 0;
}

