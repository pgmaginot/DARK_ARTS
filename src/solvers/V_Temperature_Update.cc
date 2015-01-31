#include "V_Temperature_Update.h"

V_Temperature_Update::V_Temperature_Update(const Fem_Quadrature& fem_quadrature, const Cell_Data& cell_data, Materials& material, 
    const Angular_Quadrature& angular_quadrature, const int n_stages)
  :
  m_material(material),
  m_sn_w(angular_quadrature.get_sum_w() ),  
  m_np(fem_quadrature.get_number_of_interpolation_points() ) , 
  m_n_cells(cell_data.get_total_number_of_cells() ),
  m_n_source_quad_pts(fem_quadrature.get_number_of_integration_points() ),
  m_i_matrix(Eigen::MatrixXd::Identity(m_np,m_np)),
  m_time_data_set(false),
  m_t_star_vec(Eigen::VectorXd::Zero(m_np)),
  m_t_old_vec(Eigen::VectorXd::Zero(m_np)),
  m_driving_source_vec(Eigen::VectorXd::Zero(m_np)),
  m_phi_vec(Eigen::VectorXd::Zero(m_np)),
  m_planck_vec(Eigen::VectorXd::Zero(m_np)),
  m_k_vec(Eigen::VectorXd::Zero(m_np)),
  m_r_sig_a(Eigen::MatrixXd::Zero(m_np,m_np)),    
  m_r_cv(Eigen::MatrixXd::Zero(m_np,m_np)),
  m_d_matrix(Eigen::MatrixXd::Zero(m_np,m_np) )
{
  MATRIX_INTEGRATION matrix_type = fem_quadrature.get_integration_type();   
  /// initialize matrix constructor
  if(matrix_type ==  EXACT)
  {
    m_mtrx_builder = std::make_shared<Matrix_Construction_Exact>(fem_quadrature,material );
  }
  else if(matrix_type == TRAD_LUMPING)
  {
    m_mtrx_builder = std::make_shared<Matrix_Construction_Trad_Lumping>(fem_quadrature,material) ;
  }
  else if(matrix_type == SELF_LUMPING)
  {
    m_mtrx_builder = std::make_shared<Matrix_Construction_Self_Lumping>(fem_quadrature,material );
  }

  m_rk_a.resize(n_stages,0.);
}

bool V_Temperature_Update::check_eigen_variables_finite(void) const
{
  bool good_vals = true;
  std::stringstream err;
  
  if(!m_i_matrix.allFinite() )
  {
    err << "m_i_matrix is not finite\n";
    good_vals = false;
  }
  
  std::cout << err.str() << std::endl;
  if(!good_vals)
    throw Dark_Arts_Exception(TIME_MARCHER , "Eigen variables of V_Temperature are non-finite");
  
  return good_vals;
}

void V_Temperature_Update::set_time_data( const double dt, const double time_stage, const std::vector<double>& rk_a_of_stage_i, const int stage )
{
  m_dt = dt;
  m_time = time_stage;
  m_stage = stage;
  
  for(int i = 0 ; i <= m_stage; i++)
    m_rk_a[i] = rk_a_of_stage_i[i];
    
  m_time_data_set = true;
    
  return;
}

