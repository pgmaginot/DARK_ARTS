/** @file   V_Sweep_Matrix_Creator.cc
  *   @author pmaginot
  *   @brief Implement the (partially) virtual Matrix_Construction class, 
  *     get all data from Fem_Quadrature object, calculate gradient matrix, and calculate upwinding vectors 
  *       does not change amongst the different matrix integration schemes
 */

#include "V_Sweep_Matrix_Creator.h"

V_Sweep_Matrix_Creator::V_Sweep_Matrix_Creator(const Fem_Quadrature& fem_quadrature, Materials& materials,
  const int n_stages, 
  const double sn_w, 
  const int n_l_mom,
  const Temperature_Data& t_old, 
  const Intensity_Data& i_old,
  const K_Temperature& kt, 
  const K_Intensity& ki,
  const Temperature_Data& t_star)
:
  m_matrix_type{fem_quadrature.get_integration_type() },
  m_sn_w{ sn_w },
  m_n_l_mom{ n_l_mom } , 
  m_np{fem_quadrature.get_number_of_interpolation_points()}, 
  
  /// matrices/vectors for linearization definitions
  m_r_sig_t{ Eigen::MatrixXd::Zero(m_np,m_np) },
  m_r_sig_s(m_n_l_mom,Eigen::MatrixXd::Zero(m_np,m_np)),
  m_s_i( Eigen::VectorXd(m_np,m_np) ),
  
  /// additional matrices needed to get k_I
  m_k_i_r_sig_t(Eigen::MatrixXd::Zero(m_np,m_np)),
  /// just scattering.  No absorption or re-emission or linearization terms
  m_k_i_r_sig_s_zero(Eigen::MatrixXd::Zero(m_np,m_np)),
  /// don't need a vector for S_I, m_driving_source was calculated at last call
  
  m_r_sig_a(Eigen::MatrixXd::Zero(m_np,m_np) ),
  m_identity_matrix( Eigen::MatrixXd::Identity(m_np,m_np) ),  
  m_r_cv( Eigen::MatrixXd::Zero(m_np,m_np) ),  
  m_d_matrix( Eigen::MatrixXd::Zero(m_np,m_np) ),
  m_coefficient( Eigen::MatrixXd::Zero(m_np,m_np) ),
  m_mass( Eigen::MatrixXd::Zero(m_np,m_np) ),
  m_dx_div_2_mass( Eigen::MatrixXd::Zero(m_np,m_np) ), 
  m_mass_inv( Eigen::MatrixXd::Zero(m_np,m_np) ),
  m_xi_isotropic( Eigen::VectorXd(m_np) ),
  m_driving_source( Eigen::VectorXd(m_np) ),
  
  m_t_old_vec( Eigen::VectorXd::Zero(m_np) ),
  m_t_star_vec( Eigen::VectorXd::Zero(m_np) ),
  m_k_vec( Eigen::VectorXd::Zero(m_np) ),
  m_planck_vec( Eigen::VectorXd::Zero(m_np) ),
  m_temp_vec( Eigen::VectorXd::Zero(m_np) ),
  
  m_no_mu_pos_l_matrix( Eigen::MatrixXd::Zero(m_np,m_np) ),
  m_no_mu_neg_l_matrix(Eigen::MatrixXd::Zero(m_np,m_np) ),
  m_no_mu_pos_f_vector( Eigen::VectorXd::Zero(m_np) ),
  m_no_mu_neg_f_vector( Eigen::VectorXd::Zero(m_np) ),
  
  m_materials(materials),
  
  m_c{ materials.get_c() },
  
  m_dt{-1.},  
  m_stage{-1},
  m_time{-1.},
  
  m_t_old_ref(t_old),
  m_i_old_ref(i_old),  
  m_kt_ref(kt),  
  m_ki_ref(ki),
  

  
  m_dx{-1.}  ,
  m_cell_num{-1},
  m_group_num{ -1},
  
  m_t_star_ref(t_star),
    
  m_ard_phi_ptr{nullptr},
  
  m_rk_a(n_stages,0.)
{    
  /// initialize matrix constructor
  if(m_matrix_type ==  EXACT)
  {
    m_mtrx_builder = std::shared_ptr<V_Matrix_Construction> (new    
      Matrix_Construction_Exact(fem_quadrature,materials) );
  }
  else if(m_matrix_type == TRAD_LUMPING)
  {
    m_mtrx_builder = std::shared_ptr<V_Matrix_Construction> (new 
      Matrix_Construction_Trad_Lumping(fem_quadrature,materials) );
  }
  else if(m_matrix_type == SELF_LUMPING)
  {
    m_mtrx_builder = std::shared_ptr<V_Matrix_Construction> (new 
      Matrix_Construction_Self_Lumping(fem_quadrature,materials) );
  }
  
  /// construct gradient terms once
  m_mtrx_builder->construct_pos_gradient_matrix(m_no_mu_pos_l_matrix);
  m_mtrx_builder->construct_neg_gradient_matrix(m_no_mu_neg_l_matrix);
  
  m_mtrx_builder->construct_pos_upwind_vector(m_no_mu_pos_f_vector);
  m_mtrx_builder->construct_neg_upwind_vector(m_no_mu_neg_f_vector);  
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

void V_Sweep_Matrix_Creator::set_time_data( const double dt, const double time_stage, const std::vector<double>& rk_a_of_stage_i, const int stage )
{
  m_stage = stage;
  m_time = time_stage;
  
  for(int s=0; s<m_stage ; s++)
    m_rk_a[s] = rk_a_of_stage_i[s];
 
  m_dt = dt;
  return;
}

void V_Sweep_Matrix_Creator::get_r_sig_t(Eigen::MatrixXd& r_sig_t) const 
{
  r_sig_t = m_r_sig_t;
  return;
}

void V_Sweep_Matrix_Creator::get_r_sig_s(Eigen::MatrixXd& r_sig_s, const int l_mom) const 
{
  r_sig_s = m_r_sig_s[l_mom];
  
  return;
}

void V_Sweep_Matrix_Creator::get_s_i(Eigen::VectorXd& s_i) const
{
  s_i = m_s_i;
  return;
}

void V_Sweep_Matrix_Creator::get_mass_inverse(Eigen::MatrixXd& mass_inv) const
{
  mass_inv = m_mass_inv;
  return;
}

 void V_Sweep_Matrix_Creator::set_ard_phi_ptr(Intensity_Moment_Data* ard_phi_ptr)
 {
  m_ard_phi_ptr = ard_phi_ptr;
  return;
 }
 
void V_Sweep_Matrix_Creator::calculate_k_i_quantities(void)
{
  m_k_i_r_sig_t = m_r_sig_t - 1./(m_rk_a[m_stage]*m_c*m_dt)*m_dx_div_2_mass;
  m_k_i_r_sig_s_zero = m_k_i_r_sig_t - m_r_sig_a;
  return;
}

void V_Sweep_Matrix_Creator::k_i_get_r_sig_a(Eigen::MatrixXd& r_sig_a) const
{
  r_sig_a = m_r_sig_a;
  return;
}

void V_Sweep_Matrix_Creator::k_i_get_r_sig_s_zero(Eigen::MatrixXd& r_sig_s_zero) const
{
  r_sig_s_zero = m_k_i_r_sig_s_zero;
  return;
}

void V_Sweep_Matrix_Creator::k_i_get_r_sig_t(Eigen::MatrixXd& r_sig_t) const
{
  r_sig_t = m_k_i_r_sig_t;
  return;
}

void V_Sweep_Matrix_Creator::k_i_get_s_i(Eigen::VectorXd& s_i) const
{
  s_i = m_driving_source;
  return;
}

void V_Sweep_Matrix_Creator::k_i_get_planck_vec(Eigen::VectorXd& planck) const
{
  planck = m_planck_vec;
  return;
}
