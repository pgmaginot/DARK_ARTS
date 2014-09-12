#include "Temperature_Update_MF.h"

Temperature_Update_MF::Temperature_Update_MF(const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* const material, 
    const Angular_Quadrature& angular_quadrature, const Time_Stepper& time_stepper)
  :
  V_Temperature_Update(fem_quadrature, cell_data,material, angular_quadrature, time_stepper),  
  m_n_groups{angular_quadrature.get_number_of_groups() }  ,
  m_spectrum{ Eigen::VectorXd::Zero(m_np) }
{

}

void Temperature_Update_MF::update_temperature(const Intensity_Data& intensity, 
  Temperature_Data& t_new, const Temperature_Data& t_star, const Temperature_Data& t_n,
  const K_Temperature& k_t, const int stage, const std::vector<double>& outside_rk_a, const double time, const double dt)
{
  /// load outside values into local memory
  load_rk_a(stage,outside_rk_a);
  
  for(int c=0;c<m_n_cells ; c++)
  {
    /// m_t_star is an Eigen::VectorXd
    t_star.get_cell_temperature(c,m_t_star);
    
    t_n.get_cell_temperature(c,m_t_old);
    
    
    /// Since we are going to evaluate material properties, first we need to populate local data in the Materials object
    m_material->calculate_local_temp_and_position(c, m_t_star);
    
    /** this routine will calculate 
     1. \f$ \mathbf{R}_{\sigma_a}  \f$ (m_r_sig_a)
     2. \f$ \mathbf{R}_{C_v}^{-1} \f$ (m_r_cv)
     3. \f$ \left[ \mathbf{I} + \text{m_sn_w}*\Delta t a_{ii} \mathbf{R}_{C_v}^{-1} \mathbf{R}_{\sigma_a} \mathbf{D} \right] \f$ (m_coeff_matrix)    
    */
    calculate_local_matrices(c , m_t_star ,dt, m_rk_a[stage] , time, intensity);
    
    ///calculate \f$ \vec{T}_{n} - \vec{T}^* + \Delta t \sum_{j=1}^{stage-1} a_{ij k_{T,j} \f$
    m_t_old -= m_t_star;
    for(int i=0; i< stage; i++)
    {
      k_t.get_kt(c, i, m_k_vec);
      m_t_old += dt*m_rk_a[i]*m_k_vec;
    }
           
    /** use calculate local \f$ \vec{T}_i \f$ (of stage i)
      -remember that m_driving_source holds the following term:
      \f[
        \sum_{g=1}^G{ \mathbf{R}_{\sigma_{a,g}} \left( \vec{\phi}_{g,i} - \text{m_sn_w}\vec{\widehat{B}}_g \right) } + \vec{S}_T
      \f]
    */
    m_t_new = m_t_star + m_coeff_matrix*m_t_old + dt*m_rk_a[stage]*m_coeff_matrix*m_r_cv*m_driving_source;
    
    /// save the local Eigen::VectorXd T_i in the t_new Temperature_Data object
    t_new.set_cell_temperature(c,m_t_new);
  }
  return;
}

void Temperature_Update_MF::calculate_local_matrices(const int cell , const Eigen::VectorXd& m_t_star ,
  const double dt, const double a_ii , const double time, const Intensity_Data& intensity)
{  
  /// calculate the m_spectrum matrix (planck deriviative and r_sig_a)
  /// first zero out m_spectrum
  m_spectrum = Eigen::MatrixXd::Zero(m_np,m_np);
  
  /// get temperature driving source, store all phi - planck values in this vector, to avoid calculating extra copies of m_r_sig_a
  m_mtrx_builder->construct_temperature_source_moments(m_driving_source,time);
  
  for(int g=0;g<m_n_groups;g++)
  {
    m_mtrx_builder->construct_r_sigma_a(m_r_sig_a,g);
    
    m_material->get_mf_planck_derivative(m_t_star,g,m_d_matrix);
    
    /// \f$ \sum_{g=1}^G{ \mathbf{R}_{\sigma_{a,g}} \mathbf{D}_g} \f$
    m_spectrum += m_r_sig_a*m_d_matrix;
    
    /// group g anagle integrated intensity
    intensity.get_cell_angle_integrated_intensity(cell, g, 0, m_phi);
    
    /// group g Planck integration
    m_material->get_mf_planck(m_t_star,g,m_planck);
    
    /// \f$ \sum_{g=1}^G{ \mathbf{R}_{\sigma_{a,g}} \left( \vec{\phi}_g - \text{m_sn_w} \vec{\widehat{B}}_g \right) } \f$
    m_driving_source += m_r_sig_a*( m_phi - m_sn_w*m_planck );
  }
  
  /// calculate \f$ \mathbf{R}_{C_v} \f$
  m_mtrx_builder->construct_r_cv(m_r_cv);
  
  /// temporarily store the inverse in m_coefficient_matrix;  
  m_coeff_matrix = m_r_cv.inverse();
  
  /// store \f$ \mathbf{R}_{C_v}^{-1} \f$ in m_r_cv
  m_r_cv = m_coeff_matrix;
  
  m_coeff_matrix = m_i_matrix + m_sn_w*a_ii*dt*m_r_cv*m_spectrum;  
  
  return;
}
