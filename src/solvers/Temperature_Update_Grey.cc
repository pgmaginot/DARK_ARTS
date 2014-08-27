#include "Temperature_Update_Grey.h"

Temperature_Update_Grey::Temperature_Update_Grey(const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* material,
  const Angular_Quadrature& angular_quadrature, Time_Stepper* time_stepper)
  :
  V_Temperature_Update(fem_quadrature, cell_data, material,angular_quadrature,time_stepper)
{

}

void Temperature_Update_Grey::update_temperature(const Intensity_Data& intensity, 
  Temperature_Data& t_new, const Temperature_Data& t_star, const Temperature_Data& t_n,
  const K_Temperature& k_t, const int stage, const double time, const double dt)
{
  double a_ii = m_time_stepper->get_a(stage,stage);
  for(int c=0;c<m_n_cells ; c++)
  {
    t_star.get_cell_temperature(c,m_t_star);
    calculate_local_matrices(c , m_t_star ,dt, a_ii, time);
  }
  return;
}

void  Temperature_Update_Grey::calculate_local_matrices(const int cell_num, Eigen::VectorXd& t_eval,
  const double dt, const double a_ii, const double time)
{
  /// Calculate R_cv, R_sig_a, D, and the constant matrix
  
  /// First we must calculate material properties at DFEM integration points
  m_material->calculate_local_temp_and_position(cell_num , t_eval);
  
  /// cell width
  double dx = m_cell_data->get_cell_width(cell_num);
  
  /**
    m_material has temperature and material number already after call to calculate_local_temp_and_position()
  */
  /// get sigma_a now
  /// implicitly we know that since this is grey, we want group 0
  m_material->get_sigma_a(0,m_temp_mat_vec);  
  m_mtrx_builder->construct_reaction_matrix(m_r_sig_a,m_temp_mat_vec);
  m_r_sig_a *= dx;
  
  /// get cv now
  m_material->get_cv(m_temp_mat_vec);  
  m_mtrx_builder->construct_reaction_matrix(m_r_cv,m_temp_mat_vec);
  m_r_cv *= dx;
  
  /// get derivative of planck function WRT temperature
  get_planck_derivative_matrix(t_eval);
  
  /// get planck vector
  get_planck_vector(t_eval);
  
  /// temporarily store \$f \mathbf{R}_{C_v}^{-1}  \$f
  m_coeff_matrix = m_r_cv.inverse();
  
  m_r_cv = m_coeff_matrix;
  
  m_coeff_matrix = m_i_matrix + m_sn_w*dt*a_ii*m_r_cv*m_r_sig_a*m_d_matrix;  
  
  /// get temperature driving source
  // m_material->get_
  
  return;
}

void Temperature_Update_Grey::get_planck_vector(const Eigen::VectorXd& t_eval)
{
  for(int i=0;i<m_np ; i++)
    m_planck(i) = m_material->planck.integrate_B_grey(t_eval(i));

  return;
}

void Temperature_Update_Grey::get_planck_derivative_matrix(const Eigen::VectorXd& t_eval)
{
  for(int i=0;i <m_np ; i++)
    m_d_matrix(i,i) = m_material->planck.integrate_dBdT_grey(t_eval(i));
  
  return;
}

