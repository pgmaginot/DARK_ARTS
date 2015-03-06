/** @file   Sweep_Matrix_Creator_Grey.cc
  *   @author pmaginot
  *   @brief Implement the Sweep_Matrix_Creator_Grey class (make pseudo matrices and sources for transport sweep)
  *    Creates the "funny" matrices that arise in the Planck/temperature linearization when accounting for spatially varying material properties
 */

#include "Sweep_Matrix_Creator_Grey.h"



Sweep_Matrix_Creator_Grey::Sweep_Matrix_Creator_Grey(const Fem_Quadrature& fem_quadrature, 
  Materials& materials,
  const int n_stages, 
  const double sn_w, 
  const int n_l_mom,
  const Temperature_Data& t_old, 
  const Intensity_Data& i_old,
  const K_Temperature& kt,
  const K_Intensity& ki,
  const Temperature_Data& t_star)
:
 V_Sweep_Matrix_Creator( fem_quadrature, materials, n_stages , sn_w, n_l_mom, t_old,  i_old, kt, ki, t_star ),
  m_group_num(0),
  m_hold_matrix( Eigen::MatrixXd::Zero(m_np,m_np) )
{  
}


  /// calculate \f$ \mathbf{R}_{C_v}^{-1} \f$, \f$ \mathbf{M} \f$, get \f$ \vec{T}^*,~\vec{T}_n \f$
void Sweep_Matrix_Creator_Grey::update_cell_dependencies(const int cell)
{  
  /// set cell number
  m_cell_num = cell;
  
  /// get \f$ \vec{T}^n \f$
  m_t_old_ref.get_cell_temperature(m_cell_num,m_t_old_vec) ;
  
  /// get \f$ \vec{T}^* \f$
  m_t_star_ref.get_cell_temperature(m_cell_num,m_t_star_vec) ;
    
  m_materials.calculate_local_temp_and_position(cell,m_t_star_vec);
    
  /// set cell width
  m_dx = m_materials.get_cell_width();
  
  /// calculate \f$ \mathbf{M} \f$ by getting generic mass matrix and multiplying by cell width
  m_dx_div_2_mass = (m_dx/2.)*m_mass;
  
  m_mass_inv = m_dx_div_2_mass.fullPivLu().solve(m_identity_matrix);
    
  /// load \f$ \mathbf{R}_{C_v} \f$ here
  m_mtrx_builder->construct_r_cv(m_coefficient);
  
 /// then invert it and store in a temporary matrix
  m_r_cv = m_coefficient.fullPivLu().solve(m_identity_matrix); 
  /// put this back into m_r_cv
  
  return;
}

  /**
    calculate \f$ \vec{\widehat{B}}_g, ~\mathbf{D}_G^* \f$
    calculate \f$ \bar{\bar{\mathbf{R}}}_{\sigma_{t,g}} \f$ , the isotropic components of \f$ \vec{S}_I \f$, and 
    if grey \f$ \mathbf{R}_{\sigma_{s,0}} + \bar{\bar{\nu}}\mathbf{R}_{\sigma_a} \f$ 
  */  
void Sweep_Matrix_Creator_Grey::update_group_dependencies(const int grp)
{  
  /** calculate the "coefficient matrix":
    \f[ \mathbf{I} + \text{m_sn_w} \Delta t a_{ii} \mathbf{R}_{C_v}^{-1} \mathbf{R}_{\sigma_a} \mathbf{D}
    \f]  ^{-1}  
  */
  m_mtrx_builder->construct_r_sigma_a(m_r_sig_a,m_group_num);
  
  /// get Planck vector since it won't change
  /// Planck grey_planck is \f$ acT^4 \f$  not \f$ \frac{1}{\text{m_sn_w}} acT^4 \f$
  m_materials.get_grey_planck(m_t_star_vec, m_planck_vec);
  
  /// calculate \f$ \mathbf{D} \f$
  m_materials.get_grey_planck_derivative(m_t_star_vec,m_d_matrix);
  
  m_coefficient = m_identity_matrix;
  m_hold_matrix = m_r_sig_a*m_d_matrix;
  
  Eigen::MatrixXd hold2 = Eigen::MatrixXd::Zero(m_np,m_np);
  hold2 = m_r_cv*m_hold_matrix;
  m_coefficient += m_sn_w*m_dt*m_rk_a[m_stage]*hold2;
  
  m_hold_matrix = m_coefficient.fullPivLu().solve(m_identity_matrix); //inverse();
  m_coefficient = m_hold_matrix;
  
  for(int l=0; l< m_n_l_mom ; l++)
    m_mtrx_builder->construct_r_sigma_s(m_r_sig_s,m_group_num,l);
  
  m_r_sig_t = m_r_sig_a + m_r_sig_s[0];
  
  /// calculate \f$ \bar{\bar{\mathbf{R}}}_{\sigma_t} \f$
  m_r_sig_t += 1./(m_c*m_dt*m_rk_a[m_stage])*m_dx_div_2_mass; 

  /// add \f$ \bar{\bar{\mathbf \nu}} \mathbf{R}_{\sigma_a} \f$ contribution to m_sig_s
  ///   m_r_sig_s[0] += (m_sn_w*m_dt*m_rk_a[m_stage])*m_r_sig_a*m_d_matrix*m_coefficient*m_r_cv*m_r_sig_a
  m_hold_matrix = m_r_cv*m_r_sig_a;
  
  hold2 = m_d_matrix*m_coefficient;
  Eigen::MatrixXd hold3 = Eigen::MatrixXd::Zero(m_np,m_np);
  hold3 = hold2*m_hold_matrix;
  hold2 = m_r_sig_a*hold3;
  m_r_sig_s[0] += (m_sn_w*m_dt*m_rk_a[m_stage])*hold2;
    
  /// start building the isotropic portion of \f$ \bar{\bar{\xi}} \f$
  
  /// \f$ \text{m_xi_isotropic} = \vec{T}^n - - \vec{T}^* \f$
  m_xi_isotropic = m_t_old_vec - m_t_star_vec;
  
  /// \f$ \text{m_xi_isotropic} += \Delta t \sum_{j=1}^{i-1}{a_{ij} \vec{k}_{T,j} } \f$
  for(int s=0; s< m_stage ; s++)
  {
    m_kt_ref.get_kt(m_cell_num, s, m_temp_vec);
    m_xi_isotropic += m_dt*m_rk_a[s]*m_temp_vec;
  }
  
  /// get \f$ \vec{S}_T \f$
  m_mtrx_builder->construct_temperature_source_moments(m_driving_source,m_time);
  
  /** \f$ \text{m_xi_isotropic} += 
    \Delta t a_{ii} \mathbf{R}_{C_v}^{-1} \left[ \vec{S}_T - \text{m_sn_w} \mathbf{R}_{\sigma_a} \vec{\widehat{B}} \right] \f$ 
  */
  m_temp_vec = m_driving_source - m_sn_w*m_r_sig_a*m_planck_vec;
  m_xi_isotropic += m_dt*m_rk_a[m_stage]*m_r_cv*m_temp_vec;
  
  m_hold_matrix = m_r_sig_a*m_d_matrix*m_coefficient;
  
  m_temp_vec = m_hold_matrix*m_xi_isotropic;
  
  m_xi_isotropic = m_temp_vec + m_r_sig_a*m_planck_vec;
  
  return;
}

void Sweep_Matrix_Creator_Grey::update_direction_dependencies(const int dir)
{
  m_s_i = m_xi_isotropic;
  
  m_i_old_ref.get_cell_intensity(m_cell_num, m_group_num, dir, m_temp_vec);
  
  m_s_i += (1./(m_c*m_dt*m_rk_a[m_stage]))*m_dx_div_2_mass*m_temp_vec;
  
  m_temp_vec = Eigen::VectorXd::Zero(m_np);
  for(int s=0 ; s< m_stage ; s++)
  {
    m_ki_ref.get_ki(m_cell_num,m_group_num,dir,s,m_k_vec);
    m_temp_vec += m_rk_a[s]*m_k_vec;
  }
  
  m_s_i += (1./(m_c*m_rk_a[m_stage]))*m_dx_div_2_mass*m_temp_vec;
  
  m_mtrx_builder->construct_radiation_source_moments(m_driving_source,m_time,dir,m_group_num);
  
  m_s_i += m_driving_source;
  
  return;
}
