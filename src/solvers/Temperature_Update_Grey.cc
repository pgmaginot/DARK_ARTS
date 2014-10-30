#include "Temperature_Update_Grey.h"

Temperature_Update_Grey::Temperature_Update_Grey(const Fem_Quadrature& fem_quadrature, const Cell_Data& cell_data, Materials& material, 
    const Angular_Quadrature& angular_quadrature, const int n_stages)
  :
  V_Temperature_Update(fem_quadrature, cell_data, material,angular_quadrature,n_stages)
{
  std::cout << "Made it to grey Temperature update constrcutor body" << std::endl;
}

void Temperature_Update_Grey::update_temperature(const Intensity_Moment_Data& phi, Temperature_Data& t_star, 
    const Temperature_Data& t_n, const K_Temperature& k_t, const double damping, Err_Temperature& err_temperature)
{  
  for(int cell=0; cell <m_n_cells ; cell++)
  {
    /// m_t_star_vec is an Eigen::VectorXd
    t_star.get_cell_temperature(cell,m_t_star_vec);
    
    t_n.get_cell_temperature(cell,m_t_old_vec);
    
    /// Since we are going to evaluate material properties, first we need to populate local data in the Materials object
    m_material.calculate_local_temp_and_position(cell, m_t_star_vec);
    
    /** this routine will calculate 
     1. \f$ \mathbf{R}_{\sigma_a}  \f$ (m_r_sig_a)
     2. \f$ \mathbf{R}_{C_v}^{-1} \f$ (m_r_cv)
     3. \f$ \left[ \mathbf{I} + \text{m_sn_w}*\Delta t a_{ii} \mathbf{R}_{C_v}^{-1} \mathbf{R}_{\sigma_a} \mathbf{D} \right] \f$ (m_coeff_matrix)    
    */
    calculate_local_matrices();
    
    ///calculate \f$ \vec{T}_{n} - \vec{T}^* + \Delta t \sum_{j=1}^{stage-1} a_{ij k_{T,j} \f$
    m_t_old_vec -= m_t_star_vec;
    for(int i=0; i< m_stage; i++)
    {
      k_t.get_kt(cell, i, m_k_vec);
      m_t_old_vec += m_dt*m_rk_a[i]*m_k_vec;
    }
    
    /// calculate \f$ \mathbf{R}_{\sigma_a} \left(\vec{\phi}_i - \text{m_sn_w} \mathbf{B}^*   \right) + \vec{S}_T \f$
    /// store quantity in m_phi
    phi.get_cell_angle_integrated_intensity(cell,0,0,m_phi_vec);
    m_material.get_grey_planck(m_t_star_vec,m_planck_vec);
    
    m_phi_vec -= m_sn_w * m_planck_vec;
    m_phi_vec *= m_r_sig_a*m_phi_vec;
    m_phi_vec += m_driving_source_vec;
        
    /// use all these quantites and calculate \f$ \vec{T}_i \f$
    m_delta = m_coeff_matrix*m_t_old_vec + m_dt*m_rk_a[m_stage]*m_coeff_matrix*m_r_cv*m_phi_vec;
    
    /// check convergence of temperature
    
    m_t_star_vec += damping*m_delta;
    
    /// save the local Eigen::VectorXd T_i in the t_new Temperature_Data object
    t_star.set_cell_temperature(cell,m_t_star_vec);
  }
  return;
}

void  Temperature_Update_Grey::calculate_local_matrices(void)
{
  /// Calculate R_cv, R_sig_a, D, and the constant matrix
  
  /**
    m_material has temperature and material number already after call to calculate_local_temp_and_position()
  */
  /// get sigma_a now
  /// implicitly we know that since this is grey, we want group 0
  m_mtrx_builder->construct_r_sigma_a(m_r_sig_a,0);
  
  /// get cv now
  m_mtrx_builder->construct_r_cv(m_r_cv);
  
  /// get derivative of planck function WRT temperature
  m_material.get_grey_planck_derivative(m_t_star_vec,m_d_matrix);
   
  /// temporarily store \$f \mathbf{R}_{C_v}^{-1}  \$f
  m_coeff_matrix = m_r_cv.inverse();
  
  m_r_cv = m_coeff_matrix;
  
  m_coeff_matrix = m_i_matrix + m_sn_w*m_dt*m_rk_a[m_stage]*m_r_cv*m_r_sig_a*m_d_matrix;  
  
  /// get temperature driving source
  m_mtrx_builder->construct_temperature_source_moments(m_driving_source_vec,m_time);
  
  return;
}

void Temperature_Update_Grey::calculate_k_t(const Temperature_Data& t_star, K_Temperature& k_t, const Intensity_Moment_Data& ard_phi)
{
  /**
    For the grey case:
    \f[
      k_T = \mathbf{R}_{C_v}^{-1}\left[ \mathbf{R}_{\sigma_a}\left( \vec{\phi} - \text{m_sn_w} \vec{\widehat{B}} \right) + \vec{S}_T \right]
    \f]
  */
  for(int cell = 0; cell < m_n_cells ; cell ++)
  {
    t_star.get_cell_temperature(cell,m_t_star_vec);
    
    m_material.calculate_local_temp_and_position(cell, m_t_star_vec);
    
    m_mtrx_builder->construct_r_sigma_a(m_r_sig_a,0);
    
    /// get \f$ \mathbf{R}_{C_v}^{-1} \f$ 
    m_mtrx_builder->construct_r_cv(m_coeff_matrix);   
    m_r_cv = m_coeff_matrix.inverse();
  
    /// get temperature driving source
    m_mtrx_builder->construct_temperature_source_moments(m_driving_source_vec,m_time);
    
    /// get planck and angle integrated intensity
    ard_phi.get_cell_angle_integrated_intensity(cell,0,0,m_phi_vec);
    m_material.get_grey_planck(m_t_star_vec,m_planck_vec);
    
    m_k_vec = m_r_cv*( m_r_sig_a*(m_phi_vec - m_sn_w*m_planck_vec) + m_driving_source_vec );
    
    k_t.set_kt(cell, m_stage, m_k_vec);
  }
  
  return;
}
