/** @file   V_Matrix_Construction.cc
  *   @author pmaginot
  *   @brief Implement the (partially) virtual Matrix_Construction class, 
  *     get all data from Fem_Quadrature object, calculate gradient matrix, and calculate upwinding vectors 
  *       does not change amongst the different matrix integration schemes
 */

#include "V_Matrix_Construction.h"

V_Matrix_Construction::V_Matrix_Construction(const Fem_Quadrature& fem_quadrature)
:
  m_n_quad_pts{ fem_quadrature.get_number_of_interpolation_points() }, 
  m_n_basis_pts{ fem_quadrature.get_number_of_integration_points()} 
{    
  fem_quadrature.get_integration_weights(m_integration_weights);
  
  fem_quadrature.get_dfem_at_edges(m_basis_at_left_edge,m_basis_at_right_edge);

  fem_quadrature.get_dfem_at_integration_points(m_basis_at_quad);
  fem_quadrature.get_dfem_derivatives_at_integration_points(m_basis_deriv_at_quad);
}
    
void V_Matrix_Construction::construct_pos_gradient_matrix(Eigen::MatrixXd& l_mat)
{
  double temp_sum = 0.;
  for(int i=0; i<m_n_basis_pts;i++)
  {
    for(int j=0; j<m_n_basis_pts;j++)
    {
      temp_sum =0.;
      l_mat(i,j) = m_basis_at_right_edge[i]*m_basis_at_right_edge[j];
      for(int q=0;q<m_n_quad_pts;q++)
      {  
        temp_sum += m_integration_weights[q]*
          m_basis_deriv_at_quad[q+i*m_n_quad_pts]*m_basis_at_quad[q+j*m_n_quad_pts];
      }
      l_mat(i,j) -= temp_sum;
    }
  }
  return;
}

void V_Matrix_Construction::construct_neg_gradient_matrix(Eigen::MatrixXd& l_mat)
{
  double temp_sum = 0.;
  for(int i=0; i<m_n_basis_pts;i++)
  {
    for(int j=0; j<m_n_basis_pts;j++)
    {
      temp_sum =0.;
      l_mat(i,j) = -m_basis_at_left_edge[i]*m_basis_at_left_edge[j];
      for(int q=0;q<m_n_quad_pts;q++)
      {  
        temp_sum += m_integration_weights[q]*
          m_basis_deriv_at_quad[q+i*m_n_quad_pts]*m_basis_at_quad[q+j*m_n_quad_pts];
      }
      l_mat(i,j) -= temp_sum;
    }
  }
  return;
}

void V_Matrix_Construction::construct_left_upwind_vector(Eigen::VectorXd& f_mu_pos)
{
  for(int j=0; j<m_n_basis_pts;j++)
  {
    f_mu_pos(j) = m_basis_at_left_edge[j];
  }
  return;
}
  
void V_Matrix_Construction::construct_right_upwind_vector(Eigen::VectorXd& f_mu_neg)
{
  for(int j=0; j<m_n_basis_pts;j++)
  {
    f_mu_neg(j) = m_basis_at_right_edge[j];
  }
  return;
}   

/** Use integration points already stored.  In the future, may want to have seperate quadrature that is more exact,
  because NSE article showed that exact integration of moments is more robust (and most correct)
*/
void V_Matrix_Construction::construct_source_moments(Eigen::VectorXd& source_mom, 
  std::vector<double>& source_evals)
{
  for(int j=0; j<m_n_basis_pts;j++)
  {
    source_mom(j) = 0.;
    for(int q=0;q<m_n_quad_pts;q++)
      source_mom(j) += m_integration_weights[q]*m_basis_at_quad[q+j*m_n_quad_pts]*source_evals[q];
  }
  
  return;
}