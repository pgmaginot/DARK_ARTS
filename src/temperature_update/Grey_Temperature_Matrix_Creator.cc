#include "Grey_Temperature_Matrix_Creator.h"

Grey_Temperature_Matrix_Creator::Grey_Temperature_Matrix_Creator(const Fem_Quadrature& fem_quadrature, const Cell_Data& cell_data, 
  Materials& materials, const Angular_Quadrature& angular_quadrature, const int n_stages, const Temperature_Data& t_old,
  const Intensity_Moment_Data& phi)
  :
  V_Temperature_Matrix_Creator(fem_quadrature, cell_data, materials, angular_quadrature, n_stages, t_old, phi),
  scratch1(Eigen::VectorXd::Zero(m_np) ),
  scratch2(Eigen::VectorXd::Zero(m_np) ),
  scratch_mat(Eigen::MatrixXd::Zero(m_np,m_np) )
{
  
}

void Grey_Temperature_Matrix_Creator::calculate_update_quantities(const int cell, const Eigen::VectorXd& t_star, 
  const K_Temperature& k_t,
  Eigen::MatrixXd& coefficient , Eigen::VectorXd& rhs)
{  
  rhs = Eigen::VectorXd::Zero(m_np);
  coefficient = Eigen::MatrixXd::Zero(m_np,m_np);
  
  m_t_old.get_cell_temperature(cell,m_t_old_vec);
    
  m_phi.get_cell_angle_integrated_intensity(cell,0,0,m_phi_vec);
    
  m_materials.get_grey_planck(t_star,m_planck_vec);        
    
  m_materials.get_grey_planck_derivative(t_star,m_d_matrix);
    
    /// Since we are going to evaluate material properties, first we need to populate local data in the Materials object
  m_materials.calculate_local_temp_and_position(cell, t_star);

  m_mtrx_builder->construct_r_sigma_a(m_r_sig_a,0);  
  m_mtrx_builder->construct_r_cv(m_r_cv);    
  m_mtrx_builder->construct_temperature_source_moments(m_driving_source_vec,m_time);
    
  // std::cout << "time: " << m_time << std::endl;
  // std::cout << "s_t\n" << m_driving_source_vec << std::endl;
  // std::cout << "r_cv\n" << m_r_cv << std::endl;
  // std::cout << "r_sig_a\n" << m_r_sig_a << std::endl;
  // std::cout << "d_matrix\n" << m_d_matrix << std::endl;
    
  scratch1 = Eigen::VectorXd::Zero(m_np);
  scratch2 = Eigen::VectorXd::Zero(m_np);  
  /** calculate 
   \f[
    \mathbf{R}_{C_v}^{-1} \left[ \mathbf{R}_{\sigma_a}\left(\vec{\phi} - \text{m_sn_w} \text{m_planck_vec}  \right) + \vec{S}_T  \right]
   \f]
  */
  scratch1 = m_r_sig_a*(m_phi_vec - m_sn_w*m_planck_vec) + m_driving_source_vec;
  scratch2 = m_r_cv.fullPivLu().solve(scratch1);
  scratch1 = Eigen::VectorXd::Zero(m_np);
  scratch1 = m_rk_a[m_stage]*m_dt*scratch2;
    
  ///calculate \f$ \vec{T}_{n} - \vec{T}^* + \Delta t \sum_{j=1}^{stage-1} a_{ij k_{T,j} \f$
  m_t_old_vec -= t_star;
  
  scratch2 = Eigen::VectorXd::Zero(m_np);  
  for(int i=0; i< m_stage; i++)
  {
    k_t.get_kt(cell, i, m_phi_vec);
    scratch2 += m_dt*m_rk_a[i]*m_phi_vec;
  }
  scratch2 += m_t_old_vec;

  rhs = scratch2 + scratch1;
     
  /** calculate
   \f[
    \left[ \mathbf{I} + \text{m_sn_w}\Delta t a_{ii} \mathbf{R}_{C_v}^{-1}\mathbf{R}_{\sigma_a} \mathbf{D} \right]
   \f]
  */
  scratch_mat = Eigen::MatrixXd::Zero(m_np,m_np);
  scratch_mat = m_r_cv.fullPivLu().solve(m_r_sig_a*m_d_matrix);
  
  coefficient = m_identity_matrix + m_sn_w*m_dt*m_rk_a[m_stage]*scratch_mat;

  return;
}

void Grey_Temperature_Matrix_Creator::calculate_k_t(const int cell, const Eigen::VectorXd& t_star, Eigen::VectorXd& k_t)
{
  k_t = Eigen::VectorXd::Zero(m_np);
  m_phi.get_cell_angle_integrated_intensity(cell,0,0,m_phi_vec);
    
  m_materials.get_grey_planck(t_star,m_planck_vec);        
    
    /// Since we are going to evaluate material properties, first we need to populate local data in the Materials object
  m_materials.calculate_local_temp_and_position(cell, t_star);

  m_mtrx_builder->construct_r_sigma_a(m_r_sig_a,0);  
  m_mtrx_builder->construct_r_cv(m_r_cv);    
  m_mtrx_builder->construct_temperature_source_moments(m_driving_source_vec,m_time);
  
  // std::cout << "Calculating k_t\n" ;
  // std::cout << "time: " << m_time << std::endl;
  // std::cout << "s_t\n" << m_driving_source_vec << std::endl;
  // std::cout << "r_cv\n" << m_r_cv << std::endl;
  // std::cout << "r_sig_a\n" << m_r_sig_a << std::endl;
  
  
  
  scratch1 = Eigen::VectorXd::Zero(m_np);
  scratch1 = m_r_sig_a*(m_phi_vec - m_sn_w*m_planck_vec) + m_driving_source_vec;
  // std::cout << "Terms before hitting with r_cv _inv:\n" << scratch1 << std::endl;
  
  k_t = m_r_cv.fullPivLu().solve( scratch1 );
  // std::cout << "Inner k_t: \n" << k_t << std::endl;
  // std::cout << "Leaving kt calculation\n";
  
  return;
}


