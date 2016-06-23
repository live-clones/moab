/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    (2006) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file ElementSolIndQM.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "ElementSolIndQM.hpp"
#include "MsqError.hpp"
#include "PatchData.hpp"
#include <limits>

namespace MBMesquite {

std::string ElementSolIndQM::get_name() const
{
  return "ElementSolIndQM(UserIndicator)";
}


int ElementSolIndQM::get_negate_flag() const
{
  return 1;
}

/*! \file LocalSizeQualityMetric.cpp
  \author Michael Brewer
  \date April 9, 2003
  Evaluates the corner volume (or areas) of the element corners
  attached to a given vertiex and then averages those values
  together.
*/


   //!Calculate the area of the triangle formed by the three vertices.
static inline double compute_corner_area( PatchData &pd,
                                          size_t vert_1,
                                          size_t vert_2,
                                          size_t vert_3,
                                          MsqError &err)
{
  const MsqVertex* verts = pd.get_vertex_array(err);
  Vector3D vec_1=verts[vert_2]-verts[vert_1];
  Vector3D vec_2=verts[vert_3]-verts[vert_1];
  Vector3D cross_vec=vec_1*vec_2;
  return (cross_vec.length()/2.0);
}
   
   //!Calculate the volume of the tetrahedron formed by the four vertices.
static inline double compute_corner_volume( PatchData &pd,
                                            size_t vert_1,
                                            size_t vert_2,
                                            size_t vert_3,
                                            size_t vert_4,
                                            MsqError &err)
{
  const MsqVertex* verts = pd.get_vertex_array(err);
  Vector3D vec_1=verts[vert_2]-verts[vert_1];
  Vector3D vec_2=verts[vert_3]-verts[vert_1];
  Vector3D vec_3=verts[vert_4]-verts[vert_1];
  return fabs((vec_3%(vec_1*vec_2))/6.0);

}  


/*!For the given vertex, vert, with connected elements, e_i for i=1...K,
  the ElementSolIndQM computes the corner volumes (or areas) of
  each e_i at the corner defined by vert.  The corner volume is defined
  as the volume of the tet defined by the edges of an element which contain
  the common vertex, vert.  That volume is then diveded by the average corner
  volume of all the element corners connected to this vertex.  For
  vertices attached to pyramid elements, this metric is undefined.
*/
bool ElementSolIndQM::evaluate( PatchData &pd, size_t this_vert,
                                       double &fval, MsqError &err )
{
  fval=0.0;
    //get the element array
  MsqMeshEntity* elems = pd.get_element_array(err);  MSQ_ERRZERO(err);
    //get the vertex to element array and the offset array
  //const size_t* elem_offset = pd.get_vertex_to_elem_offset(err);  MSQ_ERRZERO(err);
  //const size_t* v_to_e_array = pd.get_vertex_to_elem_array(err);  MSQ_ERRZERO(err);
    //find the offset for this vertex
  //size_t this_offset = elem_offset[this_vert];
    //get the number of elements attached to this vertex (given by the
    //first entry in the vertex to element array)
  //size_t num_elems = v_to_e_array[this_offset];
    //PRINT_INFO("\nIN LOCAL SIZE CPP, num_elements = %i",num_elems);
  size_t num_elems;
  const size_t *v_to_e_array = pd.get_vertex_element_adjacencies( this_vert, num_elems, err ); 
  MSQ_ERRZERO(err);
  
  if(num_elems <= 0){
    return true;
  }
  
    //create an array to store the local metric values before averaging
    //Can we remove this dynamic allocatio?
  double* met_vals = new double[num_elems];
    //vector to hold the other verts which form a corner.
  std::vector<size_t> other_vertices;
  other_vertices.reserve(4);
  double total_val=0.0;
  size_t i=0;
    //loop over the elements attached to this vertex
  for(i=0;i<num_elems;++i){
      //get the vertices which (with this_vert) form the corner of
      //the ith element.
    elems[v_to_e_array[i]].get_connected_vertices(this_vert,
                                                              other_vertices,
                                                              err);  MSQ_ERRZERO(err);
      ////PRINT_INFO("\nINSIDE LOCAL SIZE CPP other_vertices size = %i",other_vertices.size());
    
    switch(other_vertices.size()){
        //if a surface element, compute the corner area
      case 2:
        met_vals[i] = compute_corner_area(pd, this_vert, other_vertices[0],
                                          other_vertices[1], err);  MSQ_ERRZERO(err);
        break;
          //if a volume element, compute the corner volume 
      case 3:
        met_vals[i] = compute_corner_volume(pd, this_vert, other_vertices[0],
                                            other_vertices[1],
                                            other_vertices[2], err);  MSQ_ERRZERO(err);
        break;
      default:
          //otherwise, there is was an error.  Either the wrong number
          //of vertices were returned fom get_connected_vertices or
          //the element does not have the correct number of edges
          //connected to this vertex (possibly a pyramid element).
        met_vals[i]=0.0;
        MSQ_SETERR(err)("Incorrect number of vertices returned from "
                        "get_connected_vertices.", MsqError::UNSUPPORTED_ELEMENT);
        return false;
    };
    // Scale by the user solution indicator [0, 1]
    met_vals[i]*=m_solution_indicator[i];
      //keep track of total so that we can compute the linear average
    total_val+=met_vals[i];
      //clear the vector of other_vertices for re-use.
    other_vertices.clear();
    
  }
  //calculate fval by averaging the corner values
  fval = average_metrics(met_vals, num_elems, err);  MSQ_ERRZERO(err);
  //clean up the dynamically allocated array
  delete []met_vals;
  //always return true... the vertex position is never invalid
  return true;
}

     
bool ElementSolIndQM::evaluate_with_indices( PatchData& pd,
                                                    size_t vertex,
                                                    double& value,
                                                    std::vector<size_t>& indices,
                                                    MsqError& err )
{
  indices.clear();
  pd.get_adjacent_vertex_indices( vertex, indices, err ); MSQ_ERRZERO(err);
  
  std::vector<size_t>::iterator r, w;
  for (r = w = indices.begin(); r != indices.end(); ++r) {
    if (*r < pd.num_free_vertices()) {
      *w = *r;
      ++w;
    }
  }
  indices.erase( w, indices.end() );
  if (vertex < pd.num_free_vertices())
    indices.push_back( vertex );
  
  bool rval = evaluate( pd, vertex, value, err );
  return !MSQ_CHKERR(err) && rval;
}


} // namespace Mesquite
