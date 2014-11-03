/** @file   Sweep_Matrix_Creator_MF.cc
  *   @author pmaginot
  *   @brief Implement the Sweep_Matrix_Creator_MF class (make pseudo matrices and sources for transport sweep)
  *    Creates the "funny" matrices that arise in the Planck/temperature linearization when accounting for spatially varying material properties
 */

#include "Sweep_Matrix_Creator_MF.h"

Sweep_Matrix_Creator_MF::Sweep_Matrix_Creator_MF(const Fem_Quadrature& fem_quadrature, 
  Materials& materials,
  const int n_stages, const double sn_w, 
  const int n_l_mom,
  const int n_groups,
  const Temperature_Data& t_old,
  const Intensity_Data& i_old,
  const K_Temperature& kt, 
  const K_Intensity& ki,
  const Temperature_Data& t_star)
:
  V_Sweep_Matrix_Creator( fem_quadrature, materials, n_stages , sn_w, n_l_mom, t_old, i_old, kt, ki,t_star),
  m_n_groups{ n_groups },
  m_group_num{0},
  m_spectrum( Eigen::MatrixXd::Zero(m_np,m_np)),
  m_phi_vec( Eigen::VectorXd(m_np) ),
  m_sigma_planck( Eigen::VectorXd(m_np) ),
  m_local_ard( Eigen::VectorXd(m_np) ),
  m_group_independent_xi( Eigen::VectorXd(m_np) )
{  

}

    /// calculate \f$ \mathbf{R}_{C_v}^{-1} \f$, \f$ \mathbf{M} \f$, get \f$ \vec{T}^*,~\vec{T}_n \f$
void Sweep_Matrix_Creator_MF::update_cell_dependencies(const int cell)
{
  /// set cell number
  m_cell_num = cell;
  
  /// get \f$ \vec{T}^n \f$
  m_t_old_ref.get_cell_temperature(m_cell_num,m_t_old_vec) ;
  
  /// get \f$ \vec{T}^* \f$
  m_t_star_ref.get_cell_temperature(m_cell_num,m_t_star_vec) ;
  
  /// populate Materials object with local temperature and position to evalaute material properties
  m_materials.calculate_local_temp_and_position(cell,m_t_star_vec);
   
  /// set cell width
  m_dx = m_materials.get_cell_width();
  
  /// calculate \f$ \mathbf{M} \f$ by getting generic mass matrix and multiplying by cell width
  m_dx_div_2_mass = m_dx/2.*m_mass;
  
  m_mass_inv = m_dx_div_2_mass.inverse();
    
  /// load \f$ \mathbf{R}_{C_v} \f$ here
  m_mtrx_builder->construct_r_cv(m_r_cv);
 /// then invert it and store in a temporary matrix
  m_coefficient = m_r_cv.inverse(); 
  m_r_cv = m_coefficient;  
  

  
  /** build m_spectrum matrix ( \f$ \matbf{R}_{\sigma_{a,g}} \mathbf{D} \f$ ) by looping over groups.  
    -- also calculate \f$ \bar{\bar{\Sigma \Phi}} \f$  (m_local_ard)
    -- and \f$ \sum_{g}{\mathbf{R}_{\sigma_{a,g}} \vec{\widehat{B}}_g} \f$
  */
  m_spectrum = Eigen::MatrixXd::Zero(m_np,m_np);
  m_local_ard = Eigen::VectorXd::Zero(m_np);
  m_sigma_planck = Eigen::VectorXd::Zero(m_np);
  for(int g = 0; g < m_n_groups ; g++)
  {
    m_mtrx_builder->construct_r_sigma_a(m_r_sig_a,g);
    m_materials.get_mf_planck(m_t_star_vec, g, m_planck_vec);
    m_materials.get_mf_planck_derivative(m_t_star_vec, g, m_d_matrix);
    m_ard_phi_ptr->get_cell_angle_integrated_intensity(m_cell_num, g, 0, m_phi_vec);    
    
    m_spectrum += m_r_sig_a * m_d_matrix;
    m_local_ard += m_r_sig_a*m_phi_vec;
    m_sigma_planck += m_r_sig_a*m_planck_vec;
  }
  
  /// calculate \f$ \left[ \mathbf{I} + \text{m_sn_w}a_{ii} \Delta \mathbf{R}_{C_v}^{-1} \text{m_spectrum}  \right]
  m_coefficient = m_identity_matrix + m_sn_w*m_rk_a[m_stage]*m_dt*m_r_cv*m_spectrum;
  
  /** calculate the isotropic, group independent portions of xi:
    \f[
      \vec{T}_n - \vec{T}^* + \Delta t \sum_{i=1}^{ \text{m_stage}-1}{a_{stage,i} k_{T,j} } - 
        \text{m_sn_w} a_{ii} \Delta t \mathbf{R}_{C_v}^{-1} \left[ \sum_{ \text{m_sigma_planck} -  } \right]  
    \f]
  */  
  m_group_independent_xi = m_t_old_vec - m_t_star_vec;
  for(int s = 0; s< m_stage; s++)
  {
    m_kt_ref.get_kt(m_cell_num, s, m_temp_vec);
    m_group_independent_xi += m_dt*m_rk_a[s]*m_temp_vec;
  }
  
  m_mtrx_builder->construct_temperature_source_moments(m_driving_source,m_time);
  
  m_group_independent_xi -= m_sn_w*m_rk_a[m_stage]*m_dt*m_r_cv*( m_sigma_planck + m_driving_source );
  
  return;
}
  
  /**
    get \f$ \vec{\widehat{B}}_g, ~\mathbf{D}_G^* \f$
    calculate \f$ \bar{\bar{\mathbf{R}}}_{\sigma_{t,g}} \f$ , the isotropic components of \f$ \vec{S}_I \f$, and 
    if grey \f$ \mathbf{R}_{\sigma_{s,0}} + \bar{\bar{\nu}}\mathbf{R}_{\sigma_a} \f$ 
  */  

  void Sweep_Matrix_Creator_MF::update_group_dependencies(const int grp)
{  
  m_group_num = grp;
  m_materials.get_mf_planck_derivative(m_t_star_vec, m_group_num, m_d_matrix);  
  m_materials.get_mf_planck(m_t_star_vec, m_group_num , m_planck_vec);
  m_mtrx_builder->construct_r_sigma_a(m_r_sig_a,m_group_num);
  for(int l=0; l< m_n_l_mom ; l++)
    m_mtrx_builder->construct_r_sigma_s(m_r_sig_s,m_group_num,l);
  
  m_r_sig_t = m_r_sig_a + m_r_sig_s[0];
  /// calculate \f$ \bar{\bar{\mathbf{R}}}_{\sigma_t} \f$
  m_r_sig_t += 1./(m_c*m_dt*m_rk_a[m_stage])*m_dx_div_2_mass;  
  
  /// initalize/reset
  m_xi_isotropic = m_r_sig_a*m_d_matrix*m_coefficient*m_group_independent_xi;
  
    /// m_xi_isotropic can contain the ARD terms, \f$ \bar{\bar{\chi}}_g \bar{\bar{\nu}} \bar{\bar{\Sigma \Pih}} \f$ 
  m_xi_isotropic += m_sn_w*m_rk_a[m_stage]*m_dt*m_r_sig_a*m_d_matrix*m_coefficient*m_r_cv*m_local_ard;
  
  /// add in planck term
  m_xi_isotropic += m_r_sig_a*m_planck_vec;    
  
  return;
}

void Sweep_Matrix_Creator_MF::update_direction_dependencies(const int dir)
{

  /// complete m_s_i
  m_s_i = m_xi_isotropic;
  
  /// save \f$ \vec{S}_I \f$ in m_driving_source
  m_mtrx_builder->construct_radiation_source_moments(m_driving_source,m_time,dir,m_group_num);
  
  /// add this contribution to m_s_i as well
  m_s_i += m_driving_source;
  
  /// add in \f$ \vec{I}_n \f$ contribution
  m_i_old_ref.get_cell_intensity(m_cell_num, m_group_num, dir, m_temp_vec);  
  m_s_i += 1./(m_c*m_dt*m_rk_a[m_stage])*m_dx_div_2_mass*m_temp_vec;
  
  /// add in \f$ \vec{k}_{I,d,g,s} \f$ terms
  m_temp_vec = Eigen::VectorXd::Zero(m_np);
  for(int s=0 ; s< m_stage ; s++)
  {
    m_ki_ref.get_ki(m_cell_num,m_group_num,dir,s,m_k_vec);
    m_temp_vec += m_rk_a[s]*m_k_vec;
  }
  
  return;
}

