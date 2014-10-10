/** @file   Matrix_Construction_Trad_Lumping.cc
  *   @author pmaginot
  *   @brief Implement the Matrix_Construction_Trad_Lumping class, 
  *     just need to specify reaction matrix and mass matrix construction
  *   Important: will need to multiply by dx/2 outside of Matrix_Construction
  *    Matrices will be diagonal, after "exact" quadrature integration and then row collapse
 */

#include "Matrix_Construction_Trad_Lumping.h"
Matrix_Construction_Trad_Lumping::Matrix_Construction_Trad_Lumping(
  const Fem_Quadrature& fem_quadrature, Materials& materials)
  : V_Matrix_Construction(fem_quadrature,materials)
  {}


void Matrix_Construction_Trad_Lumping::construct_dimensionless_mass_matrix(Eigen::MatrixXd& mass) 
{  
  for(int i=0;i<m_n_basis_pts;i++)
  {
    double sum = 0.;
    for(int j=0;j<m_n_basis_pts;j++)
    {
      
      for(int q=0; q < m_n_quad_pts; q++)
      {
        sum += m_integration_weights[q]*m_basis_at_quad[q+i*m_n_quad_pts]*m_basis_at_quad[q+j*m_n_quad_pts];
      }      
    }
    mass(i,i) = sum;
  }
  return;
}
  
void Matrix_Construction_Trad_Lumping::construct_reaction_matrix(
  Eigen::MatrixXd& rx_mat, std::vector<double>& xs) 
{
  for(int i=0;i<m_n_basis_pts;i++)
  {
    double sum = 0.;
    for(int j=0;j<m_n_basis_pts;j++)
    {      
      for(int q=0; q < m_n_quad_pts; q++)
      {
        sum += m_integration_weights[q]*xs[q]*m_basis_at_quad[q+i*m_n_quad_pts]*m_basis_at_quad[q+j*m_n_quad_pts];
      }      
    }
    rx_mat(i,i) = sum;
  }
    
  return;
}

