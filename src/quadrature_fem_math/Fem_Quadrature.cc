/** @file   Fem_Quadrature.cc
  *   @author pmaginot
  *   @brief Implement the Fem_quadrature class, hold all of the various FEM bits and pieces (quadrature, integration, etc)
 */

#include "Fem_Quadrature.h"

Fem_Quadrature::Fem_Quadrature(const Input_Reader& input_reader, const Quadrule_New& quad_fun)
{  
  /// Get the DFEM interpolation points
  m_n_interpolation_points = input_reader.get_dfem_degree() + 1;
  
  m_dfem_interpolation_points.resize(m_n_interpolation_points,0.);
  m_dfem_interpolation_weights.resize(m_n_interpolation_points,0.);
  const QUADRATURE_TYPE dfem_point_type = input_reader.get_dfem_interpolation_point_type();
  switch(dfem_point_type)
  {
    case GAUSS:
    {
      quad_fun.legendre_dr_compute( m_n_interpolation_points , m_dfem_interpolation_points, m_dfem_interpolation_weights);
      break;
    }
    case LOBATTO:
    {
      quad_fun.lobatto_compute(m_n_interpolation_points , m_dfem_interpolation_points, m_dfem_interpolation_weights);
      break;
    }
    case EQUAL_SPACED:
    {
      quad_fun.ncc_compute(m_n_interpolation_points , m_dfem_interpolation_points, m_dfem_interpolation_weights);
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
    case SLXS:
    {
      n_opacity_eval_points = m_n_interpolation_points;
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
    /// evaluate xs interpolatory polynomials at dfem integration points
    evaluate_lagrange_func(m_xs_eval_points, m_integration_points, m_xs_poly_at_integration_points);
    
  }
  else if(xs_treatment == SLXS)
  {
    switch(input_reader.get_dfem_interpolation_point_type())
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
  m_int_method = input_reader.get_integration_method();
  switch(m_int_method)
  {
    case SELF_LUMPING:
    {
      m_n_integration_points = m_n_interpolation_points;
      break;
    }
    case TRAD_LUMPING:
    {
      m_n_integration_points = input_reader.get_dfem_degree() + 1 + input_reader.get_opacity_degree();
      break;
    }
    case EXACT:
    {
      /// going to use Gauss quadrature point (accurate for P <= 2*Np - 1)
      /// integrating P_dfem + P_dfem + P_xs functions
      /// P_dfem + 1 + P_xs will more than cover it
      m_n_integration_points = input_reader.get_dfem_degree() + 1 + input_reader.get_opacity_degree();
      break;
    }
    case INVALID_MATRIX_INTEGRATION:
    {
      std::cout << "Error.  Invalid matrix integration strategy in Fem_Quadrature" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  
  m_integration_points.resize(m_n_integration_points,0.);
  m_integration_weights.resize(m_n_integration_points,0.);
  
  if(m_int_method == SELF_LUMPING)
  {
    switch(dfem_point_type)
    {
      case GAUSS:
      {
        quad_fun.legendre_dr_compute( m_n_integration_points , m_integration_points, m_integration_weights);
        break;
      }
      case LOBATTO:
      {
        quad_fun.lobatto_compute(m_n_integration_points , m_integration_points, m_integration_weights);
        break;
      }
      case EQUAL_SPACED:
      {
        quad_fun.ncc_compute(m_n_integration_points , m_integration_points, m_integration_weights);
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
    quad_fun.legendre_dr_compute( m_n_integration_points , m_integration_points, m_integration_weights);
  }
  
  /// Evaluate DFEM basis functions at matrix integration points
  evaluate_lagrange_func(m_dfem_interpolation_points, m_integration_points,
    m_basis_at_integration_points);
  
  /// Evaluate derivatives of basis functions on the reference element
  evaluate_lagrange_func_derivatives(m_dfem_interpolation_points, m_integration_points,
    m_d_basis_d_s_at_integration_points);
    
  /// Evaluate DFEM basis functions at xs evaluation points
  evaluate_lagrange_func(m_dfem_interpolation_points, m_xs_eval_points,
    m_basis_at_xs_points);
    
  /// make a dummy vector to use in filling out the edge values
  std::vector<double> edge_vec(1,-1.);
  evaluate_lagrange_func(m_dfem_interpolation_points, edge_vec,m_dfem_at_left_edge);  
  evaluate_lagrange_func_derivatives(m_dfem_interpolation_points, edge_vec,m_d_dfem_d_s_at_left_edge);
  
  edge_vec[0] = 1.;
  evaluate_lagrange_func(m_dfem_interpolation_points, edge_vec,m_dfem_at_right_edge);  
  evaluate_lagrange_func_derivatives(m_dfem_interpolation_points, edge_vec,m_d_dfem_d_s_at_right_edge);
  
}

int Fem_Quadrature::get_number_of_integration_points(void) const
{
  return m_n_integration_points; 
}

int Fem_Quadrature::get_number_of_interpolation_points(void) const
{
  return m_n_interpolation_points; 
}

void Fem_Quadrature::get_xs_eval_points(std::vector<double>& xs_eval_pts) const
{
  xs_eval_pts = m_xs_eval_points;
  return;
}

 void Fem_Quadrature::get_xs_eval_weights(std::vector<double>& xs_eval_wts) const
{
  xs_eval_wts = m_xs_eval_weights;
  return;
}

void Fem_Quadrature::get_dfem_at_xs_eval_points(std::vector<double>& dfem_at_xs_pts) const
{
  dfem_at_xs_pts = m_basis_at_xs_points;
  return;
}

void Fem_Quadrature::get_dfem_integration_points(std::vector<double>& dfem_integration_pts) const
{
  dfem_integration_pts = m_integration_points;
  return;
}

void Fem_Quadrature::get_dfem_at_edges(std::vector<double>& dfem_at_left_edge,
  std::vector<double>& dfem_at_right_edge) const
{
  dfem_at_left_edge = m_dfem_at_left_edge;
  dfem_at_right_edge = m_dfem_at_right_edge;
  return;
}

int Fem_Quadrature::get_number_of_xs_point(void) const
{
  return m_n_xs_evaluation_points;
}

void Fem_Quadrature::get_xs_at_dfem_integration_points(
  std::vector<double>& xs_at_dfem_integration_pts) const
{
  xs_at_dfem_integration_pts = m_xs_poly_at_integration_points;
  return;
}

MATRIX_INTEGRATION Fem_Quadrature::get_integration_type(void) const
{
  return m_int_method;
}

void Fem_Quadrature::get_dfem_at_integration_points(std::vector<double>& dfem_at_int_pts) const
{
  dfem_at_int_pts = m_basis_at_integration_points;
  return;
}  
void Fem_Quadrature::get_dfem_derivatives_at_integration_points(std::vector<double>& dfem_deriv_at_int_pts) const
{
  dfem_deriv_at_int_pts = m_d_basis_d_s_at_integration_points;
  return;
}

void Fem_Quadrature::get_integration_weights(std::vector<double>& weights) const
{
  weights = m_integration_weights;
  return;
}


void Fem_Quadrature::evaluate_lagrange_func(const std::vector<double>& interp_points, 
  const std::vector<double>& eval_points, std::vector<double>& func_evals)
{
  const int n_eval_p = eval_points.size();
  const int n_interp_p = interp_points.size();
  
  /// allocate space for func_evals vector
  /// vector laid out as [B_1(p_1) B_1(p_2) ... B_2(p_1) ... B_N(p_N) ]
  func_evals.resize(n_eval_p*n_interp_p,0.);
  
  int pos = 0;
  /// loop over basis functions
  for(int j=0; j< n_interp_p ; j++)
  {
    /// loop over points we want to evaluate at
    for(int pnt = 0; pnt < n_eval_p ; pnt++)
    {
      double val = 1.;
      for(int k = 0; k< n_interp_p ; k++)
      {
        if(j == k)
          continue;
          
        val *= (eval_points[pnt] - interp_points[k]) / 
          (interp_points[j] - interp_points[k]);
      }
      func_evals[pos] = val;
      pos++;
    }    
  }
  
  return;
}

void Fem_Quadrature::evaluate_lagrange_func_derivatives(const std::vector<double>& interp_points, 
  const std::vector<double>& eval_points, std::vector<double>& deriv_evals)
{
  /// taken from http://www.phys.ufl.edu/~coldwell/wsteve/FDERS/The%20Lagrange%20Polynomial.htm
  const int n_eval_p = eval_points.size();
  const int n_interp_p = interp_points.size();
  
  deriv_evals.resize(n_eval_p*n_interp_p,0.);
  int cnt = 0;
  for( int p=0;p<n_eval_p;p++)
  {
    for(int m=0;m<n_interp_p; m++)
    {
      double sum_val = 0;
      for(int l=0; l<n_interp_p; l++)
      {
        if(m==l)
          continue;
          
        double inner_mult = 1.;
        for(int j = 0; j < n_interp_p; j++)
        {
          if( (j==l) || (j==m) )
            continue;
            
          inner_mult *= (eval_points[p] - interp_points[j])/(interp_points[m] - interp_points[j]);
        }
        sum_val += inner_mult/(interp_points[m] - interp_points[l]);
     
      }  
      deriv_evals[cnt] = sum_val;
      cnt++;
    }
  }
  
  return;
}
    
    