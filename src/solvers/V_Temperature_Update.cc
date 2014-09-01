#include "V_Temperature_Update.h"

V_Temperature_Update::V_Temperature_Update(const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* material, 
    const Angular_Quadrature& angular_quadrature, const Time_Stepper& time_stepper)
  :
  m_material{material},
  m_sn_w{angular_quadrature.get_sum_w() },  
  m_np{fem_quadrature.get_number_of_interpolation_points() } , 
  m_cell_data{cell_data},
  m_n_cells{cell_data->get_total_number_of_cells() },
  m_n_source_quad_pts{fem_quadrature.get_number_of_integration_points() },
  m_r_sig_a{ Eigen::MatrixXd::Zero(m_np,m_np) },
  m_r_cv{ Eigen::MatrixXd::Zero(m_np,m_np) },
  m_i_matrix{Eigen::MatrixXd::Identity(m_np,m_np)},
  m_d_matrix{ Eigen::MatrixXd::Zero(m_np,m_np)},
  m_coeff_matrix{ Eigen::MatrixXd::Zero(m_np,m_np)},
  m_phi{Eigen::VectorXd::Zero(m_np)},
  m_planck{Eigen::VectorXd::Zero(m_np)},
  m_t_old{Eigen::VectorXd::Zero(m_np)},
  m_t_star{Eigen::VectorXd::Zero(m_np)},
  m_driving_source{ Eigen::VectorXd::Zero(m_np)},
  m_t_new{ Eigen::VectorXd::Zero(m_np)},
  m_k_vec{ Eigen::VectorXd::Zero(m_np)},
    
  m_matrix_type{fem_quadrature.get_integration_type() }    
{
  /// initialize matrix constructor
  if(m_matrix_type ==  EXACT)
  {
    m_mtrx_builder = std::shared_ptr<V_Matrix_Construction> (new    
      Matrix_Construction_Exact(fem_quadrature) );
  }
  else if(m_matrix_type == TRAD_LUMPING)
  {
    m_mtrx_builder = std::shared_ptr<V_Matrix_Construction> (new 
      Matrix_Construction_Trad_Lumping(fem_quadrature) );
  }
  else if(m_matrix_type == SELF_LUMPING)
  {
    m_mtrx_builder = std::shared_ptr<V_Matrix_Construction> (new 
      Matrix_Construction_Self_Lumping(fem_quadrature) );
  }
  
  /// resize STL vectors that we use for data storage
  m_temp_mat_vec.resize(m_n_source_quad_pts,0.);

  m_rk_a.resize(time_stepper.get_number_of_stages() );
}

void V_Temperature_Update::load_rk_a(const int stage, const std::vector<double>& outside_rk_a)
{
  for(int i =0 ; i < stage; i++)
    m_rk_a[i] = outside_rk_a[i];
    
  return;
}

