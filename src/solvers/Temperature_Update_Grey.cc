#include "Temperature_Update_Grey.h"

Temperature_Update_Grey::Temperature_Update_Grey(const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* material,
  const Angular_Quadrature& angular_quadrature)
  :
  V_Temperature_Update(fem_quadrature, cell_data, material,angular_quadrature)
{

}

void Temperature_Update_Grey::update_temperature(const Intensity_Data& intensity, 
  Temperature_Data& t_new, const Temperature_Data& t_star, const Temperature_Data& t_n,
  const K_Temperature& k_t, const int stage, const double time, const double dt)
{
  return;
}

void  Temperature_Update_Grey::calculate_local_matrices(const int cell_num, Eigen::VectorXd& t_eval,
  const double dt, const double a_ii)
{
  /// Calculate R_cv, R_sig_a, D, and the constant matrix
  
  /// First we must calculate material properties at DFEM integration points
  m_material->calculate_local_temp_and_position(cell_num , t_eval);
  
  /**
    m_material has temperature and material number already after call to calculate_local_temp_and_position()
  */
  /// get sigma_a now
  /// implicitly we know that since this is grey, we want group 0
  m_material->get_sigma_a(0,m_temp_mat_vec);  
  m_mtrx_builder->construct_reaction_matrix(m_r_sig_a,m_temp_mat_vec);
  
  /// get cv now
  m_material->get_cv(m_temp_mat_vec);  
  m_mtrx_builder->construct_reaction_matrix(m_r_cv,m_temp_mat_vec);
  
  /// get derivative of planck function WRT temperature
  m_material->get_planck_derivative(m_d_matrix);
  
  /// temporarily store \$f \mathbf{R}_{C_v}^{-1}  \$f
  m_coeff_matrix = m_r_cv.inverse();
  
  m_r_cv = m_coeff_matrix;
  
  m_coeff_matrix = m_i_matrix + m_sn_w*dt*a_ii*m_r_cv*m_r_sig_a*m_d_matrix;  
  
  return;
}