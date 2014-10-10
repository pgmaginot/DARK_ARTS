#include "Temperature_Update_MF.h"

Temperature_Update_MF::Temperature_Update_MF(const Fem_Quadrature& fem_quadrature, const Cell_Data& cell_data, Materials& material, 
    const Angular_Quadrature& angular_quadrature, const int n_stages)
  :
  V_Temperature_Update(fem_quadrature, cell_data, material, angular_quadrature, n_stages),  
  m_n_groups{angular_quadrature.get_number_of_groups() }  ,
  m_spectrum{ Eigen::VectorXd::Zero(m_np) }
{

}

void Temperature_Update_MF::update_temperature(const Intensity_Moment_Data& phi, Temperature_Data& t_star, 
    const Temperature_Data& t_n, const K_Temperature& k_t, const double damping, Err_Temperature& err_temperature)
{ 
  for(int cell=0;cell<m_n_cells ; cell++)
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
    calculate_local_matrices(cell,phi);
    
    ///calculate \f$ \vec{T}_{n} - \vec{T}^* + \Delta t \sum_{j=1}^{stage-1} a_{ij k_{T,j} \f$
    m_t_old_vec -= m_t_star_vec;
    for(int i=0; i< m_stage; i++)
    {
      k_t.get_kt(cell, i, m_k_vec);
      m_t_old_vec += m_dt*m_rk_a[i]*m_k_vec;
    }
           
    /** use calculate local \f$ \vec{T}_i \f$ (of stage i)
      -remember that m_driving_source holds the following term:
      \f[
        \sum_{g=1}^G{ \mathbf{R}_{\sigma_{a,g}} \left( \vec{\phi}_{g,i} - \text{m_sn_w}\vec{\widehat{B}}_g \right) } + \vec{S}_T
      \f]
    */
    m_delta = m_coeff_matrix*m_t_old_vec + m_dt*m_rk_a[m_stage]*m_coeff_matrix*m_r_cv*m_driving_source_vec;
    m_t_star_vec += damping*m_delta; 
    
    /// save the local Eigen::VectorXd T_i in the t_new Temperature_Data object
    t_star.set_cell_temperature(cell,m_t_star_vec);
  }
  return;
}

void Temperature_Update_MF::calculate_local_matrices(const int cell, const Intensity_Moment_Data& phi )
{  
  /// calculate the m_spectrum matrix (planck deriviative and r_sig_a)
  /// first zero out m_spectrum
  m_spectrum = Eigen::MatrixXd::Zero(m_np,m_np);
  
  /// get temperature driving source, store all phi - planck values in this vector, to avoid calculating extra copies of m_r_sig_a
  m_mtrx_builder->construct_temperature_source_moments(m_driving_source_vec,m_time);
  
  for(int g=0;g<m_n_groups;g++)
  {
    m_mtrx_builder->construct_r_sigma_a(m_r_sig_a,g);
    
    m_material.get_mf_planck_derivative(m_t_star_vec,g,m_d_matrix);
    
    /// \f$ \sum_{g=1}^G{ \mathbf{R}_{\sigma_{a,g}} \mathbf{D}_g} \f$
    m_spectrum += m_r_sig_a*m_d_matrix;
    
    /// group g anagle integrated intensity
    phi.get_cell_angle_integrated_intensity(cell, g, 0, m_phi_vec);
    
    /// group g Planck integration
    m_material.get_mf_planck(m_t_star_vec,g,m_planck_vec);
    
    /// \f$ \sum_{g=1}^G{ \mathbf{R}_{\sigma_{a,g}} \left( \vec{\phi}_g - \text{m_sn_w} \vec{\widehat{B}}_g \right) } \f$
    m_driving_source_vec += m_r_sig_a*( m_phi_vec - m_sn_w*m_planck_vec );
  }
  
  /// calculate \f$ \mathbf{R}_{C_v} \f$
  m_mtrx_builder->construct_r_cv(m_r_cv);
  
  /// temporarily store the inverse in m_coefficient_matrix;  
  m_coeff_matrix = m_r_cv.inverse();
  
  /// store \f$ \mathbf{R}_{C_v}^{-1} \f$ in m_r_cv
  m_r_cv = m_coeff_matrix;
  
  m_coeff_matrix = m_i_matrix + m_sn_w*m_rk_a[m_stage]*m_dt*m_r_cv*m_spectrum;  
  
  return;
}
