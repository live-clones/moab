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
  4, 5, {{sqrt(2)-1, sqrt(2)-1,0},{1-sqrt(2), sqrt(2)-1,0},{3*sqrt(2)-1, 3*sqrt(2)-1,0},{1-3*sqrt(2), 3*sqrt(2)-1,0}}, {{0,1,5,4},{1,2,6,5},{3,7,6,2},{0,4,7,3},{4,5,6,7}}
};

bool Errfunc(double coords[3], double tol)
{
  double center[3] = {0,0,0};
  double val = pow(coords[0]-center[0],2)+pow(coords[1]-center[1],2)+pow(coords[2]-center[2],2);
  if (exp(-val) > tol)
    return true;
  else
    return false;
}

void compute_coords(double *corner_coords, double xi, double eta, double vert[3])
{
  double N[4];
   N[0] = (1-xi)*(1-eta)/4; N[1] = (1+xi)*(1-eta)/4; N[2] = (1+xi)*(1+eta)/4, N[3] = (1-xi)*(1+eta)/4;

   for (int j=0; j<4; j++)
     {
       verts[0] += N[j]*corner_coords[3*j];
       verts[1] += N[j]*corner_coords[3*j+1];
       verts[2] += N[j]*corner_coords[3*j+2];
     }

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
  error = mbImpl->load_file(filename, inmesh);MB_CHK_ERR(error);
  Range allverts, allquads;
  error = mbImpl->get_entities_by_dimension(inmesh, 0, allverts);MB_CHK_ERR(error);
  error = mbImpl->get_entities_by_dimension(inmesh, 2, allquads);MB_CHK_ERR(error);

  //Create refinement set from error indicators
  EntityHandle refset;
  error = mbImpl->create_meshset(MESHSET_SET, refset);MB_CHK_ERR(error);
  double qcoords[3];
  double tol = 0.5;
  for (Range::iterator it = allquads.begin(); it != allquads.end(); it++)
    {
      error = mbImpl->get_coords(&(*it), 1, &qcoords[0]);MB_CHK_ERR(error);
      bool select = Errfunc(qcoords, tol);
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


  EntityHandle conn[4];
  double coords[12];
  for (Range::iterator it = coarse_quads.begin(); it != coarse_quads.end(); it++)
    {
      error = mbImpl->get_connectivity(*it, 1, &conn[0]);MB_CHK_ERR(error);
      error = mbImpl->get_coords(&conn[0], 4, &coords[0]);MB_CHK_ERR(error);

      //Create new vertices
      EntityHandle nwvert;
      double x[3], vbuffer[8];

      for (int i=0; i<4; i++)
        vbuffer[i] = conn[i];

      for (int i = 0; i< quadTemp.num_verts; i++)
        {
          double xi = quadTemp.vertex_natcoords[i][0];
          double eta = quadTemp.vertex_natcoords[i][1];
          compute_coords(coords, xi, eta, x);

          error = mbImpl->create_vertex(x, nwvert);MB_CHK_ERR(error);
          vbuffer[i+4] = nwvert;
        }

      //Create new quads
      EntityHandle nwquad, nwconn[4];
      for (int i=0; i<quadTemp.num_ents; i++)
        {
          for (int k=0; k<4; k++)
            {
              int idx = quadTemp.ents_conn[i][k];
              nwconn[k] = vbuffer[idx];
            }
          error = mbImpl->create_element(MBQUAD, nwconn, 4, nwquad);MB_CHK_ERR(error);
          error = mbImpl->add_entities(finemesh, &nwquad, 1);MB_CHK_ERR(error);
        }
    }


  //Write out
  std::stringstream file;
  file <<"refined_mesh.vtk";
  error = mbImpl->write_file(file.str().c_str(), 0,0,finemesh, 1);MB_CHK_ERR(error);

  return 0;
}

