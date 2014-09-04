/** @file   V_Sweep_Matrix_Creator.cc
  *   @author pmaginot
  *   @brief Implement the (partially) virtual Matrix_Construction class, 
  *     get all data from Fem_Quadrature object, calculate gradient matrix, and calculate upwinding vectors 
  *       does not change amongst the different matrix integration schemes
 */

#include "V_Sweep_Matrix_Creator.h"

V_Sweep_Matrix_Creator::V_Sweep_Matrix_Creator(const Fem_Quadrature& fem_quadrature, Materials* const materials,
  Cell_Data* const cell_data, const int n_stages)
:
  m_matrix_type{fem_quadrature.get_integration_type() },
  m_np{fem_quadrature.get_number_of_interpolation_points()},  
  
  m_r_sig_t{ Eigen::MatrixXd::Zero(m_np,m_np) },
  m_identity_matrix{ Eigen::MatrixXd::Identity(m_np,m_np) },  
  m_r_cv{ Eigen::MatrixXd::Zero(m_np,m_np) },
  m_r_sig_a{ Eigen::MatrixXd::Zero(m_np,m_np) },
  m_d_matrix{ Eigen::MatrixXd::Zero(m_np,m_np) },
  m_coefficent{ Eigen::MatrixXd::Zero(m_np,m_np) },
  m_mass{ Eigen::MatrixXd::Zero(m_np,m_np) },
  
  m_no_mu_pos_l_matrix{ Eigen::MatrixXd::Zero(m_np,m_np) },
  m_no_mu_neg_l_matrix{ Eigen::MatrixXd::Zero(m_np,m_np) },
  m_no_mu_pos_f_vector{ Eigen::VectorXd::Zero(m_np) },
  m_no_mu_neg_f_vector{ Eigen::VectorXd::Zero(m_np) },
  
  m_materials{ materials},
  
  m_cell_data{ cell_data } ,
  
  m_c{materials->get_c() },
  
  m_dt{-1.},
  m_time{-1.},
  m_stage{-1},
  
  m_t_old{nullptr},
  m_t_star{nullptr},  
  m_kt{nullptr},
  
  m_i_old{nullptr},
  m_ki{nullptr},
  
  m_dx{-1.}  
{  
  m_rk_a.resize(n_stages-1,0.);
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
  
  m_mtrx_builder->construct_mass_matrix(m_mass);
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

void V_Sweep_Matrix_Creator::set_thermal_iteration_data(const Temperature_Data* t_eval, const Temperature_Data* t_old, 
    const K_Temperature* kt )
{
  m_t_old = t_old;
  m_t_star = t_eval;
  m_kt = kt;
  
  return;
}

void V_Sweep_Matrix_Creator::set_intensity_iteration_data(const Intensity_Data* i_old, const K_Intensity* ki)
{
  m_i_old = i_old;
  m_ki = ki;
  
  return;
}

void V_Sweep_Matrix_Creator::set_stage_data(const int stage, const std::vector<double>& rk_a, const double time)
{
  m_stage = stage;
  m_time = time;
  
  for(int s=0; s<m_stage ; s++)
    m_rk_a[s] = rk_a[s];
 
  return;
}

void V_Sweep_Matrix_Creator::set_timestep_data(const double dt)
{
  m_dt = dt;
  return;
}

void V_Sweep_Matrix_Creator::get_cell_size(const int cell)
{
  m_dx = m_cell_data->get_cell_width(cell);
  return;
}

void V_Sweep_Matrix_Creator::get_r_sig_t(Eigen::MatrixXd& r_sig_t)
{
  r_sig_t = m_r_sig_t;
  return;
}
void V_Sweep_Matrix_Creator::get_r_sig_s(Eigen::MatrixXd& r_sig_s)
{
  return;
}

void V_Sweep_Matrix_Creator::get_r_sig_s_higher_moment(const int l_mom, Eigen::MatrixXd& r_sig_s_ho)
{
  return;
}

void V_Sweep_Matrix_Creator::contstruct_sweep_matrices(const int cell, const int grp)
{
  return;
}