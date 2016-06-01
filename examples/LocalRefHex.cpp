/** @example LocalRefHex.cpp  \n
 * \brief Read a hexa mesh, and refine locally some hexas, with a 5 hex template
 *
 */

#include "moab/Core.hpp"
#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#endif
#include "MBParallelConventions.h"
#include <iostream>
#include "moab/CartVect.hpp"
#include "moab/ElemEvaluator.hpp"
#include "moab/LinearHex.hpp"

using namespace moab;
using namespace std;

string test_file_name = string(MESH_DIR) + string("/64bricks_512hex_256part.h5m");

ErrorCode SimpleRefine(Interface * mb, Range & hexas, Range & outHexas)
{
   double alfa=1./3; // 
  CartVect newNodesPar[8] = { CartVect ( -alfa, -alfa, -alfa ),
                              CartVect (  alfa, -alfa, -alfa ),
                              CartVect (  alfa,  alfa, -alfa ),
                              CartVect ( -alfa,  alfa, -alfa ),
                              CartVect ( -alfa, -alfa,  alfa ),
                              CartVect (  alfa, -alfa,  alfa ),
                              CartVect (  alfa,  alfa,  alfa ),
                              CartVect ( -alfa,  alfa,  alfa ) };
  ErrorCode rval;
  for (Range::iterator eit= hexas.begin(); eit!=hexas.end(); eit++)
  {
    EntityHandle hex= *eit;  
    ElemEvaluator ee(mb, hex, 0);
    ee.set_tag_handle(0, 0);
    ee.set_eval_set(MBHEX, LinearHex::eval_set());
    const EntityHandle * C=NULL; // old connectivity
    int nnodes;
    rval = mb->get_connectivity(hex, C, nnodes); MB_CHK_ERR(rval);
    EntityHandle N[8];
    for (int k=0; k<8; k++)
    {
      double coords[3];
      rval = ee.eval(newNodesPar[k].array(), coords);MB_CHK_ERR(rval);
      rval = mb->create_vertex(coords, N[k]);
    } 
    // build the new hexas; the original entity handle will have new nodes
     // hex-face from MBCNArrays.hpp
    //  { {0,1,5,4}, {1,2,6,5}, {2,3,7,6}, {3,0,4,7}, {0,3,2,1}, {4,5,6,7} } },
    EntityHandle  conn[48] =  { N[0], N[1], N[5], N[4],
                                C[0], C[1], C[5], C[4],  // first face

                                N[1], N[2], N[6], N[5],
                                C[1], C[2], C[6], C[5],  // second face

                                N[2], N[3], N[7], N[6],
                                C[2], C[3], C[7], C[6],  // third  face {2,3,7,6}

                                N[3], N[0], N[4], N[7],
                                C[3], C[0], C[4], C[7],  // fourth  face {3,0.4,7}

                                N[0], N[3], N[2], N[1],
                                C[0], C[3], C[2], C[1],  // fifth   face {0,3,2,1}

                                N[4], N[5], N[6], N[7],
                                C[4], C[5], C[6], C[7]  // sixth   face {4,5,6,7}
                               
                               }; // 6 new hexas, one on each face  
     for (int i=0; i<6; i++)
     {
       EntityHandle newHex;
       rval = mb->create_element(MBHEX, &conn[8*i], 8, newHex); MB_CHK_ERR(rval);
       outHexas.insert(newHex);
     }
     // reset connectivity of the old hex:
     rval = mb->set_connectivity(hex, N, 8); MB_CHK_ERR(rval);
    
  }
  return MB_SUCCESS;
}
int main(int argc, char **argv)
{

  string options;

  // Need option handling here for input filename
  if (argc > 1) {
    // User has input a mesh file
    test_file_name = argv[1];
  }  

  Interface* mb = new (std::nothrow) Core;
  if (NULL == mb)
    return 1;

    cout << "Reading file " << test_file_name << "\n with options: " << options << "\n";

  // Read the file with the specified options
  ErrorCode rval = mb->load_file(test_file_name.c_str(), 0, options.c_str());MB_CHK_ERR(rval);
  Range hexas;
  rval = mb->get_entities_by_dimension(0, 3, hexas);MB_CHK_ERR(rval);
  
  Range toRefine(hexas[0], hexas[0]+3);// only 4 hexas
  Range newHexas;
  rval = SimpleRefine(mb, toRefine, newHexas);  MB_CHK_ERR(rval);
  
  rval = mb->write_file("new_file.h5m"); MB_CHK_ERR(rval);

  delete mb;

  return 0;
}
