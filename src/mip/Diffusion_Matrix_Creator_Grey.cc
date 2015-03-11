/** @file   Diffusion_Matrix_Creator_Grey.cc
  *   @author pmaginot
  *   @brief Implement the Grey_Diffusion_Matrix_Cretator class, this accounts for the full SDRIK linearization
  *   this class will evaluate the following :
  *   \f$ \mathbf{R}_{\widetilde{\Sigma}_a} \f$, 
  *   \f$ \mathbf{R}_{\widetilde{\Sigma}_s} \f$, 
  *  and $\widetilde{D}$ at cell edges and quadrature integration points
 */

#include "Diffusion_Matrix_Creator_Grey.h"

Diffusion_Matrix_Creator_Grey::Diffusion_Matrix_Creator_Grey(const Fem_Quadrature& fem_quadrature, Materials& materials,
  const Angular_Quadrature& angular_quadrature, const Temperature_Data& t_star)
:
  V_Diffusion_Matrix_Creator(fem_quadrature,materials,angular_quadrature,t_star),
  m_rk_a_ii(-1.),
  m_dt(-1.),
  m_c_speed( materials.get_c() ),
  m_sum_sn_w(m_angular_quadrature.get_sum_w() ),
  m_mass(Eigen::MatrixXd::Zero(m_np,m_np)),
  m_dx_mass(Eigen::MatrixXd::Zero(m_np,m_np)),
  m_r_cv(Eigen::MatrixXd::Zero(m_np,m_np)),  
  m_d_matrix(Eigen::MatrixXd::Zero(m_np,m_np)),
  m_coefficient(Eigen::MatrixXd::Zero(m_np,m_np)),
  m_identity_matrix(Eigen::MatrixXd::Identity(m_np,m_np)),
  m_temporary_mat(Eigen::MatrixXd::Zero(m_np,m_np))
{    
  /// calculate dimensionless mass matrix one time
  m_mtrx_builder->construct_dimensionless_mass_matrix(m_mass);
}

void Diffusion_Matrix_Creator_Grey::set_time_data( const double dt, const double time_stage, const double sdirk_a_of_stage )
{
  m_dt = dt;
  m_rk_a_ii = sdirk_a_of_stage;
  return;
}

void Diffusion_Matrix_Creator_Grey::calculate_pseudo_r_sig_a_and_r_sig_s(Eigen::MatrixXd& r_sig_a, Eigen::MatrixXd& r_sig_s)
{
  m_dx_mass = m_dx/2.*m_mass;
  
  m_mtrx_builder->construct_r_sigma_a(m_r_sig_a,m_group_num);
  m_mtrx_builder->construct_r_sigma_s(m_r_sig_s,m_group_num,0);
  
  /// calculate R_{cv} then invert
  m_mtrx_builder->construct_r_cv(m_coefficient);  
  
  // std::cout << "Cell: " << m_cell_num << std::endl;
  // std::cout << "r_sig_a: \n" << m_r_sig_a << std::endl;
  // std::cout << "r_sig_s: \n" << m_r_sig_s << std::endl;
  // std::cout << "r_cv: \n" << m_coefficient << std::endl;
  
 /// then invert it and store in a temporary matrix
  m_r_cv = m_coefficient.fullPivLu().solve(m_identity_matrix);   
  
  /// calculate \f$ \mathbf{D} \f$
  m_materials.get_grey_planck_derivative(m_t_eval_vec,m_d_matrix);
  
  m_coefficient = m_r_sig_a*m_d_matrix;
  m_temporary_mat = m_identity_matrix + m_sum_sn_w*m_dt*m_rk_a_ii*m_r_cv*m_coefficient;
  /**
    \f[
        \text{m_coefficient} = \left[ 
              \mathbf{I} + \text{m_sn_w} \Delta t a_{ii} \mathbf{R}_{C_v}^{-1} \mathbf{R}_{\sigma_a} \mathbf{D}_*
              \right]^{-1} \mathbf{R}_{C_v}^{-1} \mathbf{R}_{\sigma_a}
    \f]
  */
  m_coefficient = m_temporary_mat.fullPivLu().solve( (m_r_cv*m_r_sig_a) );
  
  m_temporary_mat = m_d_matrix*m_coefficient;
  r_sig_s = m_sum_sn_w*m_dt*m_rk_a_ii*m_r_sig_a*m_temporary_mat + m_r_sig_s;
  
  r_sig_a = m_r_sig_s + m_r_sig_a + (1./(m_rk_a_ii*m_dt*m_c_speed))*m_dx_mass - r_sig_s;
}

void Diffusion_Matrix_Creator_Grey::evaluate_all_pseudo_d_coefficients(void)
{
  return;
}
    
