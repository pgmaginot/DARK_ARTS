#include "Temperature_Update_Grey.h"

Temperature_Update_Grey::Temperature_Update_Grey(const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* const material,
  const Angular_Quadrature& angular_quadrature, const Time_Stepper& time_stepper)
  :
  V_Temperature_Update(fem_quadrature, cell_data, material,angular_quadrature,time_stepper)
{

}

void Temperature_Update_Grey::update_temperature(const Intensity_Moment_Data& phi, 
  Temperature_Data& t_new, const Temperature_Data& t_star, const Temperature_Data& t_n,
  const K_Temperature& k_t, const int stage, const std::vector<double>& outside_rk_a,
  const double time, const double dt)
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
    calculate_local_matrices(c , m_t_star ,dt, m_rk_a[stage] , time);
    
    ///calculate \f$ \vec{T}_{n} - \vec{T}^* + \Delta t \sum_{j=1}^{stage-1} a_{ij k_{T,j} \f$
    m_t_old -= m_t_star;
    for(int i=0; i< stage; i++)
    {
      k_t.get_kt(c, i, m_k_vec);
      m_t_old += dt*m_rk_a[i]*m_k_vec;
    }
    
    /// calculate \f$ \mathbf{R}_{\sigma_a} \left(\vec{\phi}_i - \text{m_sn_w} \mathbf{B}^*   \right) + \vec{S}_T \f$
    /// store quantity in m_phi
    phi.get_cell_angle_integrated_intensity(c,0,0,m_phi);
    m_material->get_grey_planck(m_t_star,m_planck);
    
    m_phi -= m_sn_w * m_planck;
    m_phi *= m_r_sig_a*m_phi;
    m_phi += m_driving_source;
        
    /// use all these quantites and calculate \f$ \vec{T}_i \f$
    m_t_new = m_t_star + m_coeff_matrix*m_t_old + dt*m_rk_a[stage]*m_coeff_matrix*m_r_cv*m_phi;
    
    /// save the local Eigen::VectorXd T_i in the t_new Temperature_Data object
    t_new.set_cell_temperature(c,m_t_new);
  }
  return;
}

void  Temperature_Update_Grey::calculate_local_matrices(const int cell_num, Eigen::VectorXd& t_eval,
  const double dt, const double a_ii, const double time)
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
  m_material->get_grey_planck_derivative(t_eval,m_d_matrix);
   
  /// temporarily store \$f \mathbf{R}_{C_v}^{-1}  \$f
  m_coeff_matrix = m_r_cv.inverse();
  
  m_r_cv = m_coeff_matrix;
  
  m_coeff_matrix = m_i_matrix + m_sn_w*dt*a_ii*m_r_cv*m_r_sig_a*m_d_matrix;  
  
  /// get temperature driving source
  m_mtrx_builder->construct_temperature_source_moments(m_driving_source,time);
  
  return;
}
