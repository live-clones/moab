#ifndef MOAB_LINEAR_HEX_HPP
#define MOAB_LINEAR_HEX_HPP

#include "moab/Matrix3.hpp"
#include "moab/CartVect.hpp"
#include <sstream>
#include <iomanip>
#include <iostream>

namespace moab { 

namespace element_utility {

namespace {

template< typename Vector>
double normsq( const Vector & v) { 
	const double n= (v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2]); 
	return n;
}

template< typename Vector>
void vec_subtract( Vector & v, const Vector & u){ 
	for(int i = 0; i < 3; ++i){ v[ i] -= u[ i]; }
}

} //non-exported functionality

template< typename _Matrix>
class Linear_hex_map {
  public:
	typedef _Matrix Matrix;
  private:
	typedef Linear_hex_map< Matrix> Self;
  public: 
    //Constructor
    Linear_hex_map() {}
    //Copy constructor
    Linear_hex_map( const Self & /* f */ ) {}

 public:
    //Natural coordinates
    template< typename Moab, typename Entity_handle, 
	      typename Points, typename Point>
    std::pair< bool, Point> operator()( const Moab & /* moab */,
					const Entity_handle & /* h */,
					const Points & v, 
					const Point & p, 
					const double tol = 1.e-6) const {
      Point result(3, 0.0);
      solve_inverse( p, result, v);
      bool point_found = solve_inverse( p, result, v, tol) &&
                is_contained( result, tol);
      return std::make_pair( point_found, result);
    }

  private:
    //This is a hack to avoid a .cpp file and C++11
    //reference_points(i,j) will be a 1 or -1;
    //This should unroll..
    inline double reference_points( const std::size_t& i,
          				  const std::size_t& j) const{
    const double rpts[8][3] = { { -1, -1, -1 },
                                {  1, -1, -1 },
                                {  1,  1, -1 },
                                { -1,  1, -1 },
                                { -1, -1,  1 },
                                {  1, -1,  1 },
                                {  1,  1,  1 },
                                { -1,  1,  1 } };
	  return rpts[ i][ j];
    }

    template< typename Point>
    bool is_contained( const Point & p, const double tol) const{
     //just look at the box+tol here
     return ( p[0]>=-1.-tol) && (p[0]<=1.+tol) &&
            ( p[1]>=-1.-tol) && (p[1]<=1.+tol) &&
            ( p[2]>=-1.-tol) && (p[2]<=1.+tol);
    }

    template< typename Point, typename Points>
    bool solve_inverse( const Point & x, 
			Point & xi,
			const Points & points, 
			const double tol=1.e-6) const {
      const double error_tol_sqr = tol*tol;
      Point delta(3,0.0);
      xi = delta;
      evaluate( xi, points, delta);
      vec_subtract( delta, x);
      std::size_t num_iterations=0;
      #ifdef LINEAR_HEX_DEBUG
 	std::stringstream ss;
	ss << "Point: "; 
       ss << x[ 0 ] << ", " << x[ 1] 
          << ", " << x [ 2] << std::endl;
	ss << "Hex: ";
	for(int i = 0; i < 8; ++i){
 	      	ss << points[ i][ 0] << ", " << points[ i][ 1] << ", "
		   << points[ i][ 2] << std::endl;
	}
	ss << std::endl;
      #endif
      while ( normsq( delta) > error_tol_sqr) {
	#ifdef LINEAR_HEX_DEBUG
	ss << "Iter #: "  << num_iterations 
	   << " Err: " << sqrt( normsq( delta)) << " Iterate: ";
	ss << xi[ 0 ] << ", " << xi[ 1] 
		<< ", " << xi[ 2] << std::endl;
	#endif
	if( ++num_iterations >= 5){ return false; }
        Matrix J;
	jacobian( xi, points, J);
        double det = moab::Matrix::determinant3( J);
        if (fabs(det) < 1.e-10){
		#ifdef LINEAR_HEX_DEBUG
			std::cerr << ss.str();
		#endif
		#ifndef LINEAR_HEX_DEBUG
		std::cerr << x[ 0 ] << ", " << x[ 1] 
			  << ", " << x [ 2] << std::endl;
		#endif
		std::cerr << "inverse solve failure: det: " << det << std::endl;
		exit( -1);
	}
        vec_subtract( xi, moab::Matrix::inverse(J, 1.0/det) * delta);
        evaluate( xi, points, delta);
	vec_subtract( delta, x);
      }
       return true;
    }

    template< typename Point, typename Points>
    Point& evaluate( const Point & p, const Points & points, Point & f) const{ 
	typedef typename Points::value_type Vector;
	Vector result;
	for(int i = 0; i < 3; ++i){ result[ i] = 0; }
	for (unsigned i = 0; i < 8; ++i) {
	  const double N_i = (1 + p[0]*reference_points(i,0))
	                   * (1 + p[1]*reference_points(i,1))
	                   * (1 + p[2]*reference_points(i,2));
	    result += N_i * points[ i];
	}
	result *= 0.125;
	for (int i = 0; i < 3; ++i){ f[ i] = result[ i]; }
	return f;
    }

    template< typename Point, typename Field>
    double evaluate_scalar_field( const Point & p, 
				  const Field & field_values) const{
	double f = 0.0; 
        for (unsigned i = 0; i < 8; ++i) {
          const double N_i = (1 + p[0]*reference_points(i,0))
                           * (1 + p[1]*reference_points(i,1))
                           * (1 + p[2]*reference_points(i,2));
            f += N_i * field_values[ i];
        }
        f *= 0.125;
        return f;
    }

    template< typename Points, typename Field> 
    double integrate_scalar_field( const Points & points, 
				   const Field & field_values) const {
      typedef typename Points::value_type Point;
      const double gauss[1][2] = { {  2.0,           0.0          } };
      const std::size_t gauss_count = 1;
      double I(0.0);
      for(std::size_t j1 = 0; j1 < gauss_count; ++j1) {
        double x1 = gauss[j1][1];
        double w1 = gauss[j1][0];
        for(std::size_t j2 = 0; j2 < gauss_count; ++j2) {
          double x2 = gauss[j2][1];
          double w2 = gauss[j2][0];
          for(std::size_t j3 = 0; j3 < gauss_count; ++j3) {
            double x3 = gauss[j3][1];
            double w3 = gauss[j3][0];
            Point x(x1,x2,x3);
	    Matrix J;
	    jacobian( x, points, J); 
            I += evaluate_scalar_field(x,field_values)*w1*w2*w3*J.determinant();
          }
        }
      }
      return I;
    }

    template< typename Point, typename Points>
    Matrix& jacobian( const Point & p, const Points & points, Matrix & J) const{
	  J = Matrix(0.0);
	  for (std::size_t i = 0; i < 8; ++i) {
	    const double   xi_p = 1 + p[0]*reference_points(i,0);
	    const double  eta_p = 1 + p[1]*reference_points(i,1);
	    const double zeta_p = 1 + p[2]*reference_points(i,2);
	    const double dNi_dxi   = reference_points(i, 0) * eta_p * zeta_p;
	    const double dNi_deta  = reference_points(i, 1) *  xi_p * zeta_p;
	    const double dNi_dzeta = reference_points(i, 2) *  xi_p *  eta_p;
	    J(0,0) += dNi_dxi   * points[i][0];
	    J(1,0) += dNi_dxi   * points[i][1];
	    J(2,0) += dNi_dxi   * points[i][2];
	    J(0,1) += dNi_deta  * points[i][0];
	    J(1,1) += dNi_deta  * points[i][1];
	    J(2,1) += dNi_deta  * points[i][2];
	    J(0,2) += dNi_dzeta * points[i][0];
	    J(1,2) += dNi_dzeta * points[i][1];
	    J(2,2) += dNi_dzeta * points[i][2];
	  }
	  return J *= 0.125;
   }
  private:
}; //Class Linear_hex_map

}// namespace element_utility

}// namespace moab
#endif //MOAB_LINEAR_HEX_nPP
