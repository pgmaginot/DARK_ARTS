/** @file   Fem_Quadrature.cc
  *   @author pmaginot
  *   @brief Implement the Fem_quadrature class, hold all of the various FEM bits and pieces (quadrature, integration, etc)
 */

#include "Fem_Quadrature.h"
#include "Quadrule_New.h"

Fem_Quadrature::Fem_Quadrature(Input_Reader&  input_reader)
{
  /// Initialize a Quadrule object to be able to get all of the quadrature we need
  Quadrule_New quad_fun;
  
  /// Get the DFEM interpolation points
  const int n_interp_point = input_reader.get_dfem_degree() + 1;
  
  m_dfem_interpolation_points.resize(n_interp_point,0.);
  m_dfem_interpolation_weights.resize(n_interp_point,0.);
  const QUADRATURE_TYPE dfem_point_type = input_reader.get_dfem_interpolation_point_type();
  switch(dfem_point_type)
  {
    case GAUSS:
    {
      quad_fun.legendre_dr_compute( n_interp_point , m_dfem_interpolation_points, m_dfem_interpolation_weights);
      break;
    }
    case LOBATTO:
    {
      quad_fun.lobatto_compute(n_interp_point , m_dfem_interpolation_points, m_dfem_interpolation_weights);
      break;
    }
    case EQUAL_SPACED:
    {
      quad_fun.ncc_compute(n_interp_point , m_dfem_interpolation_points, m_dfem_interpolation_weights);
      break;
    }
    case INVALID_QUADRATURE_TYPE:
    {
      std::cout << "Invalid Interpolation point in quadrature calculation" << std::endl;
      exit(EXIT_FAILURE);
    }
  }  
  
  /** Get the points where we will evaluate material properties at
      These points are different than the points which will be integrated
  */
  int n_opacity_eval_points = 0;
  const OPACITY_TREATMENT xs_treatment = input_reader.get_opacity_treatment();
  switch( xs_treatment)
  {
    case MOMENT_PRESERVING:
    {
      n_opacity_eval_points = m_xs_extra_points + input_reader.get_opacity_degree();
      break;
    }
    case INTERPOLATING:
    {
      n_opacity_eval_points = 1 + input_reader.get_opacity_degree();
      break;
    }
    case INVALID_OPACITY_TREATMENT:
    {
      std::cout << "Invalid Opacity treatement in quadrature" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  
  m_xs_eval_points.resize(n_opacity_eval_points,0.);
  m_xs_eval_weights.resize(n_opacity_eval_points,0.);
  if(xs_treatment == INTERPOLATING)
  {
    switch(input_reader.get_opacity_interpolation_point_type())
    {
      case GAUSS:
      {
        quad_fun.legendre_dr_compute( n_opacity_eval_points , m_xs_eval_points, m_xs_eval_weights);
        break;
      }
      case LOBATTO:
      {
        quad_fun.lobatto_compute(n_opacity_eval_points , m_xs_eval_points, m_xs_eval_weights);
        break;
      }
      case EQUAL_SPACED:
      {
        quad_fun.ncc_compute(n_opacity_eval_points , m_xs_eval_points, m_xs_eval_weights);
        break;
      }
      case INVALID_QUADRATURE_TYPE:
      {
        std::cout << "Bad Opacity Interpolation Point Type in Fem_Quadrature" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  }
  else{
    /// moment preserving, use Gauss quad to maximize accuracy
    quad_fun.legendre_dr_compute( n_opacity_eval_points , m_xs_eval_points, m_xs_eval_weights);
  }
  
  /// Get the quadrature points we are going to use to form the matrices
  int n_integration_points = 0;
  MATRIX_INTEGRATION int_method = input_reader.get_integration_method();
  switch(int_method)
  {
    case SELF_LUMPING:
    {
      n_integration_points = n_interp_point;
      break;
    }
    case TRAD_LUMPING:
    {
      n_integration_points = input_reader.get_dfem_degree() + 1 + input_reader.get_opacity_degree();
      break;
    }
    case EXACT:
    {
      /// going to use Gauss quadrature point (accurate for P <= 2*Np - 1)
      /// integrating P_dfem + P_dfem + P_xs functions
      /// P_dfem + 1 + P_xs will more than cover it
      n_integration_points = input_reader.get_dfem_degree() + 1 + input_reader.get_opacity_degree();
      break;
    }
    case INVALID_MATRIX_INTEGRATION:
    {
      std::cout << "Error.  Invalid matrix integration strategy in Fem_Quadrature" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  
  m_integration_points.resize(n_integration_points,0.);
  m_integration_weights.resize(n_integration_points,0.);
  
  if(int_method == SELF_LUMPING)
  {
    switch(dfem_point_type)
    {
      case GAUSS:
      {
        quad_fun.legendre_dr_compute( n_integration_points , m_integration_points, m_integration_weights);
        break;
      }
      case LOBATTO:
      {
        quad_fun.lobatto_compute(n_integration_points , m_integration_points, m_integration_weights);
        break;
      }
      case EQUAL_SPACED:
      {
        quad_fun.ncc_compute(n_integration_points , m_integration_points, m_integration_weights);
        break;
      }
      case INVALID_QUADRATURE_TYPE:
      {
        std::cout << "Invalid integration points in Fem_Quadrature.cc " << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  }
  else
  {  
    quad_fun.legendre_dr_compute( n_integration_points , m_integration_points, m_integration_weights);
  }
}
