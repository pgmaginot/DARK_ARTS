/** @file   V_Matrix_Construction.cc
  *   @author pmaginot
  *   @brief Implement the (partially) virtual Matrix_Construction class, 
  *     get all data from Fem_Quadrature object, calculate gradient matrix, and calculate upwinding vectors 
  *       does not change amongst the different matrix integration schemes
 */

#include "V_Matrix_Construction.h"

V_Matrix_Construction::V_Matrix_Construction(const Fem_Quadrature& fem_quadrature, Materials* const materials_ptr)
:
  m_materials_ptr{materials_ptr},
  m_n_quad_pts{ fem_quadrature.get_number_of_interpolation_points() }, 
  m_n_basis_pts{ fem_quadrature.get_number_of_integration_points() } ,
  m_n_source_quad_pts{ fem_quadrature.get_number_of_source_points() }
{    
  m_xs_evals.resize(m_n_quad_pts,0.);
  m_source_evals.resize(m_n_source_quad_pts,0.);
  
  fem_quadrature.get_dfem_at_source_points(m_dfem_at_source_quad);
  fem_quadrature.get_source_weights(m_source_weights);
  
  fem_quadrature.get_integration_weights(m_integration_weights);
  
  fem_quadrature.get_dfem_at_edges(m_basis_at_left_edge,m_basis_at_right_edge);

  fem_quadrature.get_dfem_at_integration_points(m_basis_at_quad);
  fem_quadrature.get_dfem_derivatives_at_integration_points(m_basis_deriv_at_quad);
}

void V_Matrix_Construction::construct_r_cv(Eigen::MatrixXd& r_cv)
{
  m_materials_ptr->get_cv(m_xs_evals);
  construct_reaction_matrix(r_cv,m_xs_evals);
  
 r_cv *= m_materials_ptr->get_cell_width()/2.;
  
  return;
}

void V_Matrix_Construction::construct_r_sigma_a(Eigen::MatrixXd& r_sig_a, const int grp)
{
  m_materials_ptr->get_sigma_a(grp, m_xs_evals);
  construct_reaction_matrix(r_sig_a,m_xs_evals);
  
  r_sig_a *= m_materials_ptr->get_cell_width()/2.;
  
  return;
}
  
void V_Matrix_Construction::construct_r_sigma_s(Eigen::MatrixXd& r_sig_s, const int grp, const int l_mom)
{
  m_materials_ptr->get_sigma_s(grp, l_mom, m_xs_evals);
  construct_reaction_matrix(r_sig_s,m_xs_evals);
  
  r_sig_s *= m_materials_ptr->get_cell_width()/2.;
  
  return;
}

void V_Matrix_Construction::construct_mass_matrix(Eigen::MatrixXd& mass_mat)
{
  construct_dimensionless_mass_matrix(mass_mat);
  mass_mat *= m_materials_ptr->get_cell_width()/2.;
  return;
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

void V_Matrix_Construction::construct_temperature_source_moments(Eigen::VectorXd& s_t, const double time)
{
  m_materials_ptr->get_temperature_source(time,m_source_evals);
  construct_source_moments(s_t,m_source_evals);
  return;
}
  
void V_Matrix_Construction::construct_radiation_source_moments(Eigen::VectorXd& s_i, const double time, const int dir, const int grp)
{
  m_materials_ptr->get_intensity_source(time, grp, dir, m_source_evals);
  construct_source_moments(s_i, m_source_evals);
  
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
    for(int q=0;q<m_n_source_quad_pts;q++)
      source_mom(j) += m_source_weights[q]*m_dfem_at_source_quad[q+j*m_n_source_quad_pts]*source_evals[q];
  }
  
  source_mom *= m_materials_ptr->get_cell_width()/2.;
  
  return;
}