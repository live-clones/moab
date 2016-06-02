/*This unit test is for the uniform refinement capability based on AHF datastructures*/
#include <iostream>
#include <math.h>
#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"

#include "moab/ElemEvaluator.hpp"
#include "moab/LinearHex.hpp"
#include "moab/LinearQuad.hpp"
#include "moab/QuadraticQuad.hpp"

using namespace moab;

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

#define MAX_VERTS 16
#define MAX_CHILDS 10
#define MAX_CONN   27
struct refTemplates {
  int num_verts;
  int num_ents;
  double vertex_natcoords[MAX_VERTS][3];
  int ents_conn[MAX_CHILDS][MAX_CONN];
};

struct refTemplates quadTemp = {
  4, 
  5,
   {{-1/3.0, -1/3.0,0},{1/3.0, -1/3.0,0},{1/3.0,1/3.0,0},{-1/3.0, 1/3.0,0}},
    {{0,1,5,4},{1,2,6,5},{3,7,6,2},{0,4,7,3},{4,5,6,7}}
};

#define alfa 1./3

struct refTemplates hexTemp = {
  8,
  7,
   {
     { -alfa, -alfa, -alfa },
     {  alfa, -alfa, -alfa },
     {  alfa,  alfa, -alfa },
     { -alfa,  alfa, -alfa },
     { -alfa, -alfa,  alfa },
     {  alfa, -alfa,  alfa },
     {  alfa,  alfa,  alfa },
     { -alfa,  alfa,  alfa }
   },
   {
     {0,1,5,4, 8,9,13,12}, 
     {1,2,6,5, 9,10,14,13}, 
     {2,3,7,6, 10,11,15,14}, 
     {3,0,4,7, 11,8,12,15}, 
     {0,3,2,1, 8,11,10,9}, 
     {4,5,6,7, 12,13,14,15}, 
     {8,9,10,11, 12,13,14,15} 
   }
 };


struct refTemplates quadTempQuadratic = {
  16, // 16 new nodes; 4 new corners, 4 on internal edges, 4 on radial edges, 4 interior
  5,  // 5 quads, as in the linear case 
        /* original connectivity for quad 9 in MOAB
         { -1, -1, 0},  // 0                                  
         {  1, -1, 0},  // 1                                 3           6          2
         {  1,  1, 0},  // 2                                   20       23      19
         { -1,  1, 0},  // 3                                      12    15   11  
         {  0, -1, 0},  // 4                                   

         {  1,  0, 0},  // 5                                    
         {  0,  1, 0},  // 6                                 7 24 16     8   14 22  5 
         { -1,  0, 0},  // 7                                                      
         {  0,  0, 0}   // 8  */                                
   {                                                                
     { -alfa, -alfa, 0 },  // 9                                    9    13   10
     {  alfa, -alfa, 0 },  // 10                               17       21       18
     {  alfa,  alfa, 0 },  // 11                             0           4          1    
     { -alfa,  alfa, 0 },  // 12                               
     {     0, -alfa, 0 },  // 13
     {  alfa,     0, 0 },  // 14
     {     0,  alfa, 0 },  // 15
     { -alfa,     0, 0 },  // 16
     { (-1-alfa)/2, (-1-alfa)/2, 0 }, // 17  = (0+9)/2
     { ( 1+alfa)/2, (-1-alfa)/2, 0 }, // 18  = (1+10)/2
     { ( 1+alfa)/2, ( 1+alfa)/2, 0 }, // 19  = (2+11)/2
     { (-1-alfa)/2, ( 1+alfa)/2, 0 }, // 20  = (3+12)/2
     {           0, (-1-alfa)/2, 0 }, // 21  = (4+13)/2
     { ( 1+alfa)/2,           0, 0 }, // 22  = (5+14)/2
     {           0, ( 1+alfa)/2, 0 }, // 23  = (6+15)/2
     { (-1-alfa)/2,           0, 0 }, // 24  = (7+16)/2

   },

    {
     {0,1,10, 9, 4, 18, 13, 17, 21 },
     {1,2,11,10, 5, 19, 14, 18, 22 },
     {2,3, 12, 11, 6, 20, 15, 19, 23},
     {3,0, 9,12, 7, 17, 16, 20, 24 },
     {9, 10, 11, 12, 13, 14, 15, 16, 8} 
    }
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

int main(int argc, char *argv[])
{

  Core mb;
  Interface* mbImpl = &mb;
  ErrorCode error;


  if (argc<=2)
    {
      std::cerr << "Usage: " << argv[0] << " [filename] [dim]" << std::endl;
      return 1;
    }

  //Load Mesh
  const char *filename = argv[1];
  // 
  int dim = atoi(argv[2]);
  if (dim!=2 && dim!=3) 
  {
     std::cerr<< " dimension is not 2 or 3 :" << dim << "\n";
  }
  EntityHandle inmesh;
  error = mbImpl->create_meshset(MESHSET_SET, inmesh);MB_CHK_ERR(error);
  error = mbImpl->load_file(filename, &inmesh);MB_CHK_ERR(error);
  Range ents;

  error = mbImpl->get_entities_by_dimension(inmesh, dim, ents);MB_CHK_ERR(error);
  const struct refTemplates * entTemplate;
  EntityType entType = MBQUAD;
  // number of nodes of first entity quad 4 or 9?
  EntityHandle firstEnt = ents[0];
  const EntityHandle * co;
  int numNodes = 4;
  error = mbImpl->get_connectivity(firstEnt, co, numNodes); MB_CHK_ERR(error);
  if (2==dim)
  { 
    if (numNodes == 4)
      entTemplate = &quadTemp;
    else if (numNodes ==9)
      entTemplate = &quadTempQuadratic;
    else
      std::cerr << " not supported yet\n";
  }
  else
  {
    entTemplate = &hexTemp;
    entType = MBHEX;
  }
  //Create refinement set from error indicators
  EntityHandle refset;
  error = mbImpl->create_meshset(MESHSET_SET, refset);MB_CHK_ERR(error);
  double centroid[3];
  double tol = 0.5;
  for (Range::iterator it = ents.begin(); it != ents.end(); it++)
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

  Range coarseEnts;
  error = mbImpl->get_entities_by_dimension(refset, dim, coarseEnts);MB_CHK_ERR(error);

  Range remEnts = subtract(ents, coarseEnts);
  error = mbImpl->add_entities(finemesh, remEnts);MB_CHK_ERR(error);

  EntityHandle vbuffer[25]; // max number of nodes in template // quad 9 has 9 + 16 = 25 nodes
  for (Range::iterator it = coarseEnts.begin(); it != coarseEnts.end(); it++)
    {
      int nconn = 0;
      const EntityHandle *conn;
      EntityHandle ent=*it;
      error = mbImpl->get_connectivity(ent, conn, nconn);MB_CHK_ERR(error);

      ElemEvaluator ee(mbImpl, ent, 0);
      ee.set_tag_handle(0, 0);
      if (dim == 2)
      {
        if (nconn==4)
           ee.set_eval_set(MBQUAD, LinearQuad::eval_set());
        else if (nconn == 9)
           ee.set_eval_set(MBQUAD, QuadraticQuad::eval_set());
        else
           std::cerr << " not supported yet\n";
      }
      else
        ee.set_eval_set(MBHEX, LinearHex::eval_set());
      
      
      for (int k=0; k<entTemplate->num_verts; k++) // 4 or 8
      {
        double coords[3];
        error = ee.eval(entTemplate->vertex_natcoords[k], coords);MB_CHK_ERR(error);
        error = mbImpl->create_vertex(coords, vbuffer[nconn+k]);MB_CHK_ERR(error);
      } 

      for (int i=0; i<nconn; i++)
        vbuffer[i] = conn[i];

      //Create new ents
      EntityHandle newEnt;
      EntityHandle nwconn[27]; // max could be 27 
      for (int i=0; i<entTemplate->num_ents; i++)
        {
          for (int k=0; k<nconn; k++)
            {
              int idx = entTemplate->ents_conn[i][k];
              nwconn[k] = vbuffer[idx];
            }
          error = mbImpl->create_element(entType, nwconn, nconn, newEnt);MB_CHK_ERR(error);
          error = mbImpl->add_entities(finemesh, &newEnt, 1);MB_CHK_ERR(error);
        }
    }


  //Write out
  std::stringstream file;
  file <<"refined_mesh.h5m";
  error = mbImpl->write_file(file.str().c_str(), 0, NULL, &finemesh, 1);MB_CHK_ERR(error);

  return 0;
}

