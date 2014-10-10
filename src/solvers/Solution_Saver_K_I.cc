#include "Solution_Saver_K_I.h"

Solution_Saver_K_I::Solution_Saver_K_I(const Fem_Quadrature& fem_quadrature, 
  std::shared_ptr<V_Sweep_Matrix_Creator> matrix_creator_ptr,
  const Angular_Quadrature& angular_quadrature, 
  K_Intensity& k_i_ref, const double c )
:
  V_Solution_Saver(fem_quadrature,angular_quadrature),
  m_sweep_matrix_ptr{matrix_creator_ptr},
  m_k_i_ref(k_i_ref),
  m_local_ki{Eigen::VectorXd::Zero(m_np)},
  m_scratch_vec{Eigen::VectorXd::Zero(m_np)},
  m_scratch_mat{Eigen::MatrixXd::Zero(m_np,m_np)}  ,
  m_zero{0},
  m_c{c}
{
 
}


void Solution_Saver_K_I::save_local_solution(Intensity_Moment_Data& phi_new, 
  const Eigen::VectorXd& local_intensity, 
  Psi_In& psi_in,
  const int cell, 
  const int grp, 
  const int dir)
{
  /**
    \f[
      k_I = c \mathbf{M}^{-1}\left[ \frac{1}{\text{m_sn_w}}\mathbf{R}_{\sigma_s} \vec{\phi} + \mathbf{R}_{\sigma_a}\vec{\widehat{B}} 
      - \mathbf{L}\vec{I} + I_{in} \vec{f} + \vec{S}_I - \mathbf{R}_{\sigma_t} \vec{I}  \right]    
    \f]
    -- this means that using the local solution, phi_new, we need to prepare in 
  */
  /// get 0 moment scatter source
  phi_new.get_cell_angle_integrated_intensity(cell,grp, 0 , m_scratch_vec);
  m_sweep_matrix_ptr->k_i_get_r_sig_s_zero(m_scratch_mat);  
  m_local_ki = m_scratch_mat*m_scratch_vec;
  
  /// add in higher order scattering moments
  for(int l=1; l<m_n_l_mom;l++)
  {
    phi_new.get_cell_angle_integrated_intensity(cell,grp,l, m_scratch_vec);
    m_sweep_matrix_ptr->get_r_sig_s(m_scratch_mat,l);  
    m_local_ki += m_quad_ref.get_leg_poly(dir,l)*m_scratch_mat*m_scratch_vec;
  }
  
  /// add in planck term
  m_sweep_matrix_ptr->k_i_get_planck_vec(m_scratch_vec);
  m_sweep_matrix_ptr->k_i_get_r_sig_a(m_scratch_mat);
  m_local_ki += m_scratch_mat*m_scratch_vec;
  
  /// subtract gradient matrix
  m_sweep_matrix_ptr->construct_l_matrix(m_quad_ref.get_mu(dir),m_scratch_mat);
  m_local_ki -= m_scratch_mat*local_intensity;
  
  /// add in upwind contribution
  m_sweep_matrix_ptr->construct_f_vector(m_quad_ref.get_mu(dir),m_scratch_vec); 
  m_local_ki += psi_in(dir,grp)*m_scratch_vec;
  
  /// add in driving source
  m_sweep_matrix_ptr->k_i_get_s_i(m_scratch_vec);
  m_local_ki += m_scratch_vec;
  
  /// subtract removal term
  m_sweep_matrix_ptr->k_i_get_r_sig_t(m_scratch_mat);
  m_local_ki -= m_scratch_mat*local_intensity;
  
  /// apply \f$ c\mathbf{M}^{-1} \f$
  m_sweep_matrix_ptr->get_mass_inverse(m_scratch_mat);
  m_local_ki *= m_c*m_scratch_mat;
  
  /// save the local k_i in K_Intensity object
  m_k_i_ref.set_ki(cell,grp,dir,m_stage,m_local_ki);
  
  /// calculate cell outflow (next cell's inflow)
  psi_in(dir,grp) = calculate_outflow(dir,local_intensity);
  return;
}


