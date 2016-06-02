#include "moab/QuadraticQuad.hpp"
#include "moab/Forward.hpp"

namespace moab 
{
    
      // those are not just the corners, but for simplicity, keep this name
      //
    const int QuadraticQuad::corner[9][2] = {
        { -1, -1 },
        {  1, -1 },
        {  1,  1 },  // corner nodes: 0-3
        { -1,  1 },  // mid-edge nodes: 4-7

        {  0, -1 }, //
        {  1,  0 }, //
        {  0,  1 }, //
        { -1,  0 }, //
        {  0,  0 } // center ndoe
    };

    double QuadraticQuad::SH(const int i, const double params)
    {
      switch (i)
      {
        case -1: return (params*params-params)/2;
        case 0: return 1-params*params;
        case 1: return (params*params+params)/2;
        default: return 0.;
      }
    }
    double QuadraticQuad::DSH(const int i, const double params)
    {
      switch (i)
      {
        case -1: return params-0.5;
        case 0: return -2*params;
        case 1: return params+0.5;
        default: return 0.;
      }
    }

    ErrorCode QuadraticQuad::evalFcn(const double *params, const double *field, const int /*ndim*/, const int num_tuples,
                                    double */*work*/, double *result)
    {
      assert(params && field && num_tuples > 0);
      std::fill(result, result+num_tuples, 0.0);
      for (int i=0; i<9; i++)
      {
        const double sh = SH(corner[i][0], params[0]) * SH(corner[i][1], params[1]) ;
        for (int j = 0; j < num_tuples; j++) 
          result[j] += sh * field[num_tuples*i+j];
      }

      return MB_SUCCESS;
    }

    ErrorCode QuadraticQuad::jacobianFcn(const double *params, const double *verts, const int nverts, const int ndim,
                                        double */*work*/, double *result)
    {
      assert(9 == nverts && params && verts);
      if (9 != nverts) return MB_FAILURE;
      Matrix3 *J = reinterpret_cast<Matrix3*>(result);
      *J = Matrix3(0.0);
      for (int i=0; i<9; i++)
      {
        const double sh[2]={ SH(corner[i][0], params[0]),
                             SH(corner[i][1], params[1])};
        const double dsh[2]={ DSH(corner[i][0], params[0]),
                              DSH(corner[i][1], params[1])};


        for (int j=0; j<2; j++)
        {
          (*J)(j,0)+=dsh[0]*sh[1]*verts[ndim*i+j]; // dxj/dr first column
          (*J)(j,1)+=sh[0]*dsh[1]*verts[ndim*i+j]; // dxj/ds
        }
      }
      (*J)(2,2) = 1; //
      
      return MB_SUCCESS;
    }

    ErrorCode QuadraticQuad::integrateFcn(const double */*field*/, const double */*verts*/, const int /*nverts*/, const int /*ndim*/, const int /*num_tuples*/,
                                         double */*work*/, double */*result*/)
    {
      return MB_NOT_IMPLEMENTED;
    }

    ErrorCode QuadraticQuad::reverseEvalFcn(EvalFcn eval, JacobianFcn jacob, InsideFcn ins,
                                           const double *posn, const double *verts, const int nverts, const int ndim,
                                           const double iter_tol, const double inside_tol, double *work, 
                                           double *params, int *is_inside) 
    {
      assert(posn && verts);
      return EvalSet::evaluate_reverse(eval, jacob, ins, posn, verts, nverts, ndim, iter_tol, inside_tol, 
                                       work, params, is_inside);
    } 

    int QuadraticQuad::insideFcn(const double *params, const int ndim, const double tol)
    {
      return EvalSet::inside_function(params, ndim, tol);
    }

    ErrorCode QuadraticQuad::normalFcn(const int /*ientDim*/, const int /*facet*/, const int /*nverts*/, const double */*verts*/,  double * /*normal[3]*/)
    {
      return MB_NOT_IMPLEMENTED;
    }
} // namespace moab
