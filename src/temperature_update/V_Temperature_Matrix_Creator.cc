#include "V_Temperature_Matrix_Creator.h"

V_Temperature_Matrix_Creator::V_Temperature_Matrix_Creator(const Fem_Quadrature& fem_quadrature, const Cell_Data& cell_data, 
  Materials& materials, const Angular_Quadrature& angular_quadrature, const int n_stages, const Temperature_Data& t_old,
  const Intensity_Moment_Data& phi)
  :
  m_np(fem_quadrature.get_number_of_interpolation_points() ),
  m_dt(-1.),
  m_time(-1.),
  m_stage(-1),
  m_materials(materials),
  m_sn_w(angular_quadrature.get_sum_w() ), 
  m_t_old(t_old) , 
  m_phi(phi),
  m_identity_matrix(Eigen::MatrixXd::Identity(m_np,m_np)),
  m_d_matrix(Eigen::MatrixXd::Zero(m_np,m_np)),
  m_r_sig_a(Eigen::MatrixXd::Zero(m_np,m_np)),
  m_r_cv(Eigen::MatrixXd::Zero(m_np,m_np)),  
  m_t_old_vec(Eigen::VectorXd::Zero(m_np)),
  m_driving_source_vec(Eigen::VectorXd::Zero(m_np)),
  m_phi_vec(Eigen::VectorXd::Zero(m_np)),
  m_planck_vec(Eigen::VectorXd::Zero(m_np)),
  m_k_vec(Eigen::VectorXd::Zero(m_np)),
  m_time_data_set(false)
{
  MATRIX_INTEGRATION matrix_type = fem_quadrature.get_integration_type();   
  /// initialize matrix constructor
  if(matrix_type ==  EXACT)
  {
    m_mtrx_builder = std::make_shared<Matrix_Construction_Exact>(fem_quadrature,materials );
  }
  else if(matrix_type == TRAD_LUMPING)
  {
    m_mtrx_builder = std::make_shared<Matrix_Construction_Trad_Lumping>(fem_quadrature,materials) ;
  }
  else if(matrix_type == SELF_LUMPING)
  {
    m_mtrx_builder = std::make_shared<Matrix_Construction_Self_Lumping>(fem_quadrature,materials );
  }

  m_rk_a.resize(n_stages,0.);
}

void V_Temperature_Matrix_Creator::set_time_data(const std::vector<double>& rk_a, const double dt, const double time_stage, const int stage_num)
{
  m_dt = dt;
  m_time = time_stage;
  m_stage = stage_num;
  
  for(int i = 0 ; i <= stage_num ; i++)
    m_rk_a[i] = rk_a[i];
    
  m_time_data_set=true;
  return;
}
