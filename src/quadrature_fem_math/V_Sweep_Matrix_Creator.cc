/** @file   V_Sweep_Matrix_Creator.cc
  *   @author pmaginot
  *   @brief Implement the (partially) virtual Matrix_Construction class, 
  *     get all data from Fem_Quadrature object, calculate gradient matrix, and calculate upwinding vectors 
  *       does not change amongst the different matrix integration schemes
 */

#include "V_Sweep_Matrix_Creator.h"

V_Sweep_Matrix_Creator::V_Sweep_Matrix_Creator(const Fem_Quadrature& fem_quadrature)
:
  m_matrix_type{fem_quadrature.get_integration_type() },
  m_np{fem_quadrature.get_number_of_interpolation_points()},
  m_no_mu_pos_l_matrix{ Eigen::MatrixXd::Zero(m_np,m_np) },
  m_no_mu_neg_l_matrix{ Eigen::MatrixXd::Zero(m_np,m_np) },
  m_no_mu_pos_f_vector{ Eigen::VectorXd::Zero(m_np) },
  m_no_mu_neg_f_vector{ Eigen::VectorXd::Zero(m_np) }
{  
  /// initialize matrix constructor
  if(m_matrix_type ==  EXACT)
  {
    m_mtrx_builder = std::shared_ptr<V_Matrix_Construction> (new    
      Matrix_Construction_Exact(fem_quadrature) );
  }
  else if(m_matrix_type == TRAD_LUMPING)
  {
    m_mtrx_builder = std::shared_ptr<V_Matrix_Construction> (new 
      Matrix_Construction_Trad_Lumping(fem_quadrature) );
  }
  else if(m_matrix_type == SELF_LUMPING)
  {
    m_mtrx_builder = std::shared_ptr<V_Matrix_Construction> (new 
      Matrix_Construction_Self_Lumping(fem_quadrature) );
  }
  
  /// construct gradient terms once
  m_mtrx_builder->construct_pos_gradient_matrix(m_no_mu_pos_l_matrix);
  m_mtrx_builder->construct_neg_gradient_matrix(m_no_mu_neg_l_matrix);
  
  m_mtrx_builder->construct_left_upwind_vector(m_no_mu_pos_f_vector);
  m_mtrx_builder->construct_right_upwind_vector(m_no_mu_neg_f_vector);
}

void V_Sweep_Matrix_Creator::construct_l_matrix(const double mu, Eigen::MatrixXd& l_matrix)
{
  if(mu > 0.)
  {
    l_matrix = mu*m_no_mu_pos_l_matrix;
  }
  else
  {
    l_matrix = mu*m_no_mu_neg_l_matrix;
  }
  
  return;
}

void V_Sweep_Matrix_Creator::construct_f_vector(const double mu, Eigen::VectorXd& f_vector)
{
  if(mu > 0.)
  {
    f_vector = mu*m_no_mu_pos_f_vector;
  }
  else
  {
    f_vector = mu*m_no_mu_neg_f_vector;
  }
  
  return;
}

