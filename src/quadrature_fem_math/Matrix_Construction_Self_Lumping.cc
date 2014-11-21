/** @file   Matrix_Construction_Self_Lumping.cc
  *   @author pmaginot
  *   @brief Implement the Matrix_Construction_Self_Lumping class, 
  *     just need to specify reaction matrix and mass matrix construction
  *   Important: will need to multiply by dx/2 outside of Matrix_Construction
 */

#include "Matrix_Construction_Self_Lumping.h"
Matrix_Construction_Self_Lumping::Matrix_Construction_Self_Lumping(
  const Fem_Quadrature& fem_quadrature, Materials& materials)
    : V_Matrix_Construction(fem_quadrature,materials)
    {}

void Matrix_Construction_Self_Lumping::construct_dimensionless_mass_matrix(Eigen::MatrixXd& mass) 
{
  mass = Eigen::MatrixXd::Zero(m_n_basis_pts,m_n_basis_pts);
  for(int i=0;i<m_n_basis_pts;i++)
    mass(i,i) = m_integration_weights[i];
    
  return;
}
  
void Matrix_Construction_Self_Lumping::construct_reaction_matrix(
  Eigen::MatrixXd& rx_mat, std::vector<double>& xs) 
{
  rx_mat = Eigen::MatrixXd::Zero(m_n_basis_pts,m_n_basis_pts);
  for(int i=0;i<m_n_basis_pts;i++)
    rx_mat(i,i) = m_integration_weights[i]*xs[i];
    
  return;
}

