#include <iostream>
#include <limits>
#include <assert.h>

#include "ElemUtil.hpp"
#include "Matrix3.hpp"
#include "types.h"

namespace moab { 
namespace ElemUtil {

/**\brief Class representing a 3-D mapping function (e.g. shape function for volume element) */
class VolMap {
  public:
      /**\brief Return $\vec \xi$ corresponding to logical center of element */
    virtual CartVect center_xi() const = 0;
      /**\brief Evaluate mapping function (calculate $\vec x = F($\vec \xi)$ )*/
    virtual CartVect evaluate( const CartVect& xi ) const = 0;
      /**\brief Evaluate Jacobian of mapping function */
    virtual Matrix3 jacobian( const CartVect& xi ) const = 0;
      /**\brief Evaluate inverse of mapping function (calculate $\vec \xi = F^-1($\vec x)$ )*/
    bool solve_inverse( const CartVect& x, CartVect& xi, double tol ) const ;
};

bool VolMap::solve_inverse( const CartVect& x, CartVect& xi, double tol ) const
{
  const double error_tol_sqr = tol*tol;
  double det;
  xi = center_xi();
  CartVect delta = evaluate(xi) - x;
  Matrix3 J;
  while (delta % delta > error_tol_sqr) {
    J = jacobian(xi);
    det = J.determinant();
    if (det < std::numeric_limits<double>::epsilon())
      return false;
    xi -= J.inverse(1.0/det) * delta;
    delta = evaluate( xi ) - x;
  }
  return true;
}

/**\brief Shape function for trilinear hexahedron */
class LinearHexMap : public VolMap {
  public:
    LinearHexMap( const CartVect* corner_coords ) : corners(corner_coords) {}
    virtual CartVect center_xi() const;
    virtual CartVect evaluate( const CartVect& xi ) const;
    virtual Matrix3 jacobian( const CartVect& xi ) const;
  private:
    const CartVect* corners;
    static const double corner_xi[8][3];
};

const double LinearHexMap::corner_xi[8][3] = { { -1, -1, -1 },
                                               {  1, -1, -1 },
                                               {  1,  1, -1 },
                                               { -1,  1, -1 },
                                               { -1, -1,  1 },
                                               {  1, -1,  1 },
                                               {  1,  1,  1 },
                                               { -1,  1,  1 } };
CartVect LinearHexMap::center_xi() const
  { return CartVect(0.0); }

CartVect LinearHexMap::evaluate( const CartVect& xi ) const
{
  CartVect x(0.0);
  for (unsigned i = 0; i < 8; ++i) {
    const double N_i = (1 + xi[0]*corner_xi[i][0])
                     * (1 + xi[1]*corner_xi[i][1])
                     * (1 + xi[2]*corner_xi[i][2]);
    x += N_i * corners[i];
  }
  x *= 0.125;
  return x;
}

Matrix3 LinearHexMap::jacobian( const CartVect& xi ) const
{
  Matrix3 J(0.0);
  for (unsigned i = 0; i < 8; ++i) {
    const double   xi_p = 1 + xi[0]*corner_xi[i][0];
    const double  eta_p = 1 + xi[1]*corner_xi[i][1];
    const double zeta_p = 1 + xi[2]*corner_xi[i][2];
    const double dNi_dxi   = corner_xi[i][0] * eta_p * zeta_p;
    const double dNi_deta  = corner_xi[i][1] *  xi_p * zeta_p;
    const double dNi_dzeta = corner_xi[i][2] *  xi_p *  eta_p;
    J(0,0) += dNi_dxi   * corners[i][0];
    J(1,0) += dNi_dxi   * corners[i][1];
    J(2,0) += dNi_dxi   * corners[i][2];
    J(0,1) += dNi_deta  * corners[i][0];
    J(1,1) += dNi_deta  * corners[i][1];
    J(2,1) += dNi_deta  * corners[i][2];
    J(0,2) += dNi_dzeta * corners[i][0];
    J(1,2) += dNi_dzeta * corners[i][1];
    J(2,2) += dNi_dzeta * corners[i][2];
  }
  return J *= 0.125;
}

bool nat_coords_trilinear_hex( const CartVect* corner_coords,
                               const CartVect& x,
                               CartVect& xi,
                               double tol )
{
  return LinearHexMap( corner_coords ).solve_inverse( x, xi, tol );
}


//
// nat_coords_trilinear_hex2
//  Duplicate functionality of nat_coords_trilinear_hex using hex_findpt
// 
void nat_coords_trilinear_hex2(const CartVect hex[8], 
                               const CartVect& xyz,
                               CartVect &ncoords,
                               double etol)       

{
  const int ndim = 3;
  const int nverts = 8;
  const int vertMap[nverts] = {0,1,3,2, 4,5,7,6}; //Map from nat to lex ordering

  const int n = 2; //linear
  real coords[ndim*nverts]; //buffer

  real *xm[ndim];
  for(int i=0; i<ndim; i++)
    xm[i] = coords + i*nverts;
    
  //stuff hex into coords
  for(int i=0; i<nverts; i++){
    real vcoord[ndim];
    hex[i].get(vcoord);
   
    for(int d=0; d<ndim; d++)
      coords[d*nverts + vertMap[i]] = vcoord[d];
    
  }

  double dist = 0.0;
  ElemUtil::hex_findpt(xm, n, xyz, ncoords, dist);
  if (3*EPS < dist) {
      // outside element, set extremal values to something outside range
    for (int j = 0; j < 3; j++) {
      if (ncoords[j] < (-1.0-etol) || ncoords[j] > (1.0+etol))
        ncoords[j] *= 10;
    }
  }
  
}

bool point_in_trilinear_hex(const CartVect *hex, 
                            const CartVect& xyz,
                            double etol) 
{
  CartVect xi;
  return nat_coords_trilinear_hex( hex, xyz, xi, etol )
      && fabs(xi[0])-1 < etol 
      && fabs(xi[1])-1 < etol 
      && fabs(xi[2])-1 < etol;
}


bool point_in_trilinear_hex(const CartVect *hex, 
                            const CartVect& xyz, 
                            const CartVect& box_min, 
                            const CartVect& box_max,
                            double etol) 
{
    // all values scaled by 2 (eliminates 3 flops)
  const CartVect mid = box_max + box_min;
  const CartVect dim = box_max - box_min;
  const CartVect pt = 2*xyz - mid;
  return fabs(pt[0]) - dim[0] < etol &&
         fabs(pt[1]) - dim[1] < etol &&
         fabs(pt[2]) - dim[2] < etol &&
         point_in_trilinear_hex( hex, xyz, etol );
}



// Wrapper to James Lottes' findpt routines
// hex_findpt
// Find the parametric coordinates of an xyz point inside
// a 3d hex spectral element with n nodes per dimension
// xm: coordinates fields, value of x,y,z for each of then n*n*n gauss-lobatto nodes. Nodes are in lexicographical order (x is fastest-changing)
// n: number of nodes per dimension -- n=2 for a linear element
// xyz: input, point to find
// rst: output: parametric coords of xyz inside the element. If xyz is outside the element, rst will be the coords of the closest point
// dist: output: distance between xyz and the point with parametric coords rst
extern "C"{
#include "types.h"
#include "poly.h"
#include "tensor.h"
#include "findpt.h"
#include "extrafindpt.h"
#include "errmem.h"
}

void hex_findpt(real *xm[3],
                int n,
                CartVect xyz,
                CartVect &rst,
                double &dist)       
{

  //compute stuff that only depends on the order -- could be cached
  real *z[3];
  lagrange_data ld[3];
  opt_data_3 data;

  //triplicates
  for(int d=0; d<3; d++){
    MALLOC(z[d], n, real);
    lobatto_nodes(z[d], n); 
    lagrange_setup(&ld[d], z[d], n);
  }

  opt_alloc_3(&data, ld);

  //find nearest point
  real x_star[3];
  xyz.get(x_star);

  real r[3] = {0, 0, 0 }; // initial guess for parametric coords
  unsigned c = opt_no_constraints_3;
  dist = opt_findpt_3(&data, (const real **)xm, x_star, r, &c);
  //c tells us if we landed inside the element or exactly on a face, edge, or node

  //copy parametric coords back
  rst = r;

  //Clean-up (move to destructor if we decide to cache)
  opt_free_3(&data);  
  for(int d=0; d<3; ++d) 
    lagrange_free(&ld[d]);
  for(int d=0; d<3; ++d) 
    free(z[d]);
}




// hex_eval
// Evaluate a field in a 3d hex spectral element with n nodes per dimension, at some given parametric coordinates
// field: field values for each of then n*n*n gauss-lobatto nodes. Nodes are in lexicographical order (x is fastest-changing)
// n: number of nodes per dimension -- n=2 for a linear element
// rst: input: parametric coords of the point where we want to evaluate the field
// value: output: value of field at rst

void hex_eval(real *field,
	      int n,
	      CartVect rstCartVec,
	      double &value)       
{
  int d;
  real rst[3];
  rstCartVec.get(rst);

  //can cache stuff below
  lagrange_data ld[3]; 
  real *z[3];
  for(d=0;d<3;++d){
    MALLOC(z[d], n, real);
    lobatto_nodes(z[d], n);
    lagrange_setup(&ld[d], z[d], n);
  } 

  //cut and paste -- see findpt.c
  const unsigned 
    nf = n*n,
    ne = n,
    nw = 2*n*n + 3*n;
  real *od_work;
  MALLOC(od_work, 6*nf + 9*ne + nw, real);

  //piece that we shouldn't want to cache
  for(d=0; d<3; d++){
    lagrange_0(&ld[d], rst[d]);
  }
  
  value = tensor_i3(ld[0].J,ld[0].n,
		    ld[1].J,ld[1].n,
		    ld[2].J,ld[2].n,
		    field,
		    od_work);

  //all this could be cached
  for(d=0; d<3; d++){
    free(z[d]);
    lagrange_free(&ld[d]); 
  }
  free(od_work);
}

} // namespace ElemUtil

} // namespace moab
