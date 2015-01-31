#include "Temperature_Update_Grey.h"

Temperature_Update_Grey::Temperature_Update_Grey(const Fem_Quadrature& fem_quadrature, const Cell_Data& cell_data, Materials& material, 
    const Angular_Quadrature& angular_quadrature, const int n_stages)
  :
  V_Temperature_Update(fem_quadrature, cell_data, material,angular_quadrature,n_stages)
{

}

void Temperature_Update_Grey::update_temperature(const Intensity_Moment_Data& phi, Temperature_Data& t_star, 
    const Temperature_Data& t_n, const K_Temperature& k_t, const double damping, Err_Temperature& err_temperature)
{  
  double worst_err = -1.;
  const double small_num = err_temperature.get_small_number() ;
  
  if(!m_time_data_set)
    throw Dark_Arts_Exception(TIME_MARCHER , "Must set stage number and time_stage for each call to update_temperature");
  
  
  
  for(int cell=0; cell <m_n_cells ; cell++)
  {
    t_star.get_cell_temperature(cell,m_t_star_vec);    
    t_n.get_cell_temperature(cell,m_t_old_vec);
    
    phi.get_cell_angle_integrated_intensity(cell,0,0,m_phi_vec);
    m_material.get_grey_planck(m_t_star_vec,m_planck_vec);    
    
    
    m_material.get_grey_planck_derivative(m_t_star_vec,m_d_matrix);
    
    /// Since we are going to evaluate material properties, first we need to populate local data in the Materials object
    m_material.calculate_local_temp_and_position(cell, m_t_star_vec);

    
    m_mtrx_builder->construct_r_sigma_a(m_r_sig_a,0);  
    m_mtrx_builder->construct_r_cv(m_r_cv);    
    m_mtrx_builder->construct_temperature_source_moments(m_driving_source_vec,m_time);
            
    ///calculate \f$ \vec{T}_{n} - \vec{T}^* + \Delta t \sum_{j=1}^{stage-1} a_{ij k_{T,j} \f$
    Eigen::VectorXd scratch1 = Eigen::VectorXd::Zero(m_np);
    Eigen::VectorXd scratch2 = Eigen::VectorXd::Zero(m_np);
    
    scratch1 = m_r_sig_a*(m_phi_vec - m_planck_vec) + m_driving_source_vec;
    scratch2 = m_r_cv.fullPivLu().solve(scratch1);
    scratch1 = m_rk_a[m_stage]*m_dt*scratch2;
    
    
    m_t_old_vec -= m_t_star_vec;
    scratch2 = Eigen::VectorXd::Zero(m_np);
    for(int i=0; i< m_stage; i++)
    {
      k_t.get_kt(cell, i, m_phi_vec);
      scratch2 += m_dt*m_rk_a[i]*m_phi_vec;
    }
    scratch2 += m_t_old_vec;
        
    Eigen::MatrixXd scratch_mat = Eigen::MatrixXd(m_np,m_np);
    scratch_mat = m_r_cv.fullPivLu().solve(m_r_sig_a*m_d_matrix);
    
    Eigen::MatrixXd coefficient = Eigen::MatrixXd::Zero(m_np,m_np);
    coefficient = m_i_matrix + m_sn_w*m_dt*m_rk_a[m_stage]*scratch_mat;
    
    Eigen::VectorXd delta_vec = Eigen::VectorXd::Zero(m_np);
    scratch2 += scratch1;
    delta_vec = coefficient.fullPivLu().solve( scratch2);
    
    /**
      Use a relative temperature difference of m_t_star_vec and m_delta
    */
    for(int el = 0 ; el < m_np ; el++)
    {
      double err;
      if( fabs( m_t_star_vec(el) > small_num) )
      {
        err = fabs(delta_vec(el)/m_t_star_vec(el));
      }
      else
      {
        err = fabs(delta_vec(el));
      }
      
      if(err > worst_err)
      {
        
        err_temperature.set_error(cell,el,err,delta_vec);
        worst_err = err;
      }
    }
    
    m_t_star_vec += damping*delta_vec;
    
    /// save the local Eigen::VectorXd T_i in the t_new Temperature_Data object
    t_star.set_cell_temperature(cell,m_t_star_vec);
  }
  

  return;
}

void Temperature_Update_Grey::calculate_k_t(const Temperature_Data& t_star, K_Temperature& k_t, const Intensity_Moment_Data& ard_phi)
{
  m_time_data_set = false;
  
  /**
    For the grey case:
    \f[
      k_T = \mathbf{R}_{C_v}^{-1}\left[ \mathbf{R}_{\sigma_a}\left( \vec{\phi} - \text{m_sn_w} \vec{\widehat{B}} \right) + \vec{S}_T \right]
    \f]
  */
  for(int cell = 0; cell < m_n_cells ; cell++)
  {    
    t_star.get_cell_temperature(cell,m_t_star_vec);
    
    m_material.calculate_local_temp_and_position(cell, m_t_star_vec);
    
    m_mtrx_builder->construct_r_sigma_a(m_r_sig_a,0);
    
    /// get \f$ \mathbf{R}_{C_v}^{-1} \f$ 
    m_mtrx_builder->construct_r_cv(m_r_cv);   
  
    /// get temperature driving source
    m_mtrx_builder->construct_temperature_source_moments(m_driving_source_vec,m_time);
    
    /// get planck and angle integrated intensity
    ard_phi.get_cell_angle_integrated_intensity(cell,0,0,m_phi_vec);
    m_material.get_grey_planck(m_t_star_vec,m_planck_vec);
    
    Eigen::VectorXd scratch;
    scratch = m_r_sig_a*( m_phi_vec - m_sn_w*m_planck_vec ) + m_driving_source_vec;
    
    m_k_vec = m_r_cv.fullPivLu().solve(scratch);
    
    k_t.set_kt(cell, m_stage, m_k_vec);    
  }
  check_eigen_variables_finite();
  
  return;
}
