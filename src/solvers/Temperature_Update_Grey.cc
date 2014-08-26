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

void  Temperature_Update_Grey::calculate_local_matrices(const int cell_num, Eigen::VectorXd& t_eval)
{
  /// Calculate R_cv, R_sig_a, D, and the constant matrix
  
  /// First we must calculate material properties at DFEM integration points
  m_material->calculate_local_temp_and_position(cell_num , t_eval);
  
  /// get sigma_a now
  m_material->update_sigma_a();
  /// implicitly we know that since this is grey, we want group 0
  m_material->get_sigma_a(cell_num,0,m_temp_mat_vec);  
  m_mtrx_builder->construct_reaction_matrix(m_r_sig_a,m_temp_mat_vec);
  
  /// get cv now
  m_material->update_cv();
  m_material->get_cv(cell_num,m_temp_mat_vec);  
  m_mtrx_builder->construct_reaction_matrix(m_r_cv,m_temp_mat_vec);
  
  
  return;
}