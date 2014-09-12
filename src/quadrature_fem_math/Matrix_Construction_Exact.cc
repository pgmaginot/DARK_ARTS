/** @file   Matrix_Construction_Exact.cc
  *   @author pmaginot
  *   @brief Implement the Matrix_Construction_Exact class, 
  *     just need to specify reaction matrix and mass matrix construction
  *   Important: will need to multiply by dx/2 outside of Matrix_Construction
  *    Matrices will be dense!
 */

#include "Matrix_Construction_Exact.h"

Matrix_Construction_Exact::Matrix_Construction_Exact(const Fem_Quadrature& fem_quadrature,  Materials* const materials_ptr)
  : 
  V_Matrix_Construction(fem_quadrature, materials_ptr)
{

}

void Matrix_Construction_Exact::construct_dimensionless_mass_matrix(
  Eigen::MatrixXd& mass_mat) 
{
  for(int i=0;i<m_n_basis_pts;i++)
  {
    for(int j=0;j<m_n_basis_pts;j++)
    {
      double sum = 0.;
      for(int q=0; q < m_n_quad_pts; q++)
      {
        sum += m_integration_weights[q]*m_basis_at_quad[q+i*m_n_quad_pts]*m_basis_at_quad[q+j*m_n_quad_pts];
      }
      mass_mat(i,j) = sum;
    }
  }
  return;
}
  
void Matrix_Construction_Exact::construct_reaction_matrix(
  Eigen::MatrixXd& rx_mat, std::vector<double>& xs) 
{
  for(int i=0;i<m_n_basis_pts;i++)
  {
    for(int j=0;j<m_n_basis_pts;j++)
    {
      double sum = 0.;
      for(int q=0; q < m_n_quad_pts; q++)
      {
        sum += m_integration_weights[q]*xs[q]*m_basis_at_quad[q+i*m_n_quad_pts]*m_basis_at_quad[q+j*m_n_quad_pts];
      }
      rx_mat(i,j) = sum;
    }
  }
    
  return;
}

