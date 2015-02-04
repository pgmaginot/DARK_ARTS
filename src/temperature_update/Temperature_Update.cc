#include "Temperature_Update.h"

Temperature_Update::Temperature_Update(const Fem_Quadrature& fem_quadrature, const Cell_Data& cell_data, 
  Materials& materials, const Angular_Quadrature& angular_quadrature, const int n_stages, const Temperature_Data& t_old,
  const Intensity_Moment_Data& phi)
  :  
  m_np(fem_quadrature.get_number_of_interpolation_points() ) , 
  m_n_cells(cell_data.get_total_number_of_cells() ),
  m_time_data_set(false),
  m_stage(-1),
  m_materials(materials),
  m_t_star_vec(Eigen::VectorXd::Zero(m_np)),  
  m_delta_vec(Eigen::VectorXd::Zero(m_np)),  
  m_rhs_vec(Eigen::VectorXd::Zero(m_np)),  
  m_k_vec(Eigen::VectorXd::Zero(m_np)),
  m_coefficient_matrix(Eigen::MatrixXd::Zero(m_np,m_np))
{
  if( angular_quadrature.get_number_of_groups() > 1 )
  {
    m_matrix_creator = std::make_shared<MF_Temperature_Matrix_Creator>(fem_quadrature,cell_data, materials, 
      angular_quadrature, n_stages, t_old, phi);
  }
  else{
    m_matrix_creator = std::make_shared<Grey_Temperature_Matrix_Creator>(fem_quadrature,cell_data, materials, 
      angular_quadrature, n_stages, t_old, phi);
  }
}

bool Temperature_Update::check_eigen_variables_finite(void) const
{
  bool good_vals = true;
  std::stringstream err;
    
  if(!m_t_star_vec.allFinite() )
  {
    err << "m_t_star_vec is not finite\n";
    good_vals = false;
  }
  
  if(!m_delta_vec.allFinite() )
  {
    err << "m_delta_vec is not finite\n";
    good_vals = false;
  }
  
  if(!m_rhs_vec.allFinite() )
  {
    err << "m_rhs_vec is not finite\n";
    good_vals = false;
  }
  
  if(!m_coefficient_matrix.allFinite() )
  {
    err << "m_coefficient is not finite\n";
    good_vals = false;
  }
  
  std::cout << err.str() << std::endl;
  if(!good_vals)
    throw Dark_Arts_Exception(TIME_MARCHER , "Eigen variables of Temperature_Update are non-finite");
  
  return good_vals;
}

void Temperature_Update::set_time_data( const double dt, const double time_stage, const std::vector<double>& rk_a_of_stage_i, const int stage )
{
  m_matrix_creator->set_time_data(rk_a_of_stage_i, dt, time_stage, stage);
  m_stage = stage;
  m_time_data_set=true;
  return;
}

void Temperature_Update::update_temperature(Temperature_Data& t_star, const K_Temperature& k_t, 
  const double damping, Err_Temperature& err_temperature)
{
  double worst_err = -1.;
  const double small_num = err_temperature.get_small_number() ;
  
  if(!m_time_data_set)
    throw Dark_Arts_Exception(TIME_MARCHER , "Must set stage number and time_stage for each call to update_temperature");
    
  for(int cell=0; cell <m_n_cells ; cell++)
  {
    t_star.get_cell_temperature(cell,m_t_star_vec);    
    
    m_matrix_creator->calculate_update_quantities(cell,m_t_star_vec,k_t,m_coefficient_matrix,m_rhs_vec);
    
    m_delta_vec = m_coefficient_matrix.fullPivLu().solve(m_rhs_vec);
    
    /**
      Use a relative temperature difference of m_t_star_vec and m_delta
    */
    for(int el = 0 ; el < m_np ; el++)
    {
      double err;
      if( fabs( m_t_star_vec(el) > small_num) )
      {
        err = fabs(m_delta_vec(el)/m_t_star_vec(el));
      }
      else
      {
        err = fabs(m_delta_vec(el));
      }
      
      if(err > worst_err)
      {        
        err_temperature.set_error(cell,el,err,m_delta_vec);
        worst_err = err;
      }
    }
    
    m_t_star_vec += damping*m_delta_vec;
    
    /// save the local Eigen::VectorXd T_i in the t_new Temperature_Data object
    t_star.set_cell_temperature(cell,m_t_star_vec);
  }
  
  return;
}
    

void Temperature_Update::calculate_k_t(const Temperature_Data& t_star, K_Temperature& k_t)
{  
  /**
    For the grey case:
    \f[
      k_T = \mathbf{R}_{C_v}^{-1}\left[ \mathbf{R}_{\sigma_a}\left( \vec{\phi} - \text{m_sn_w} \vec{\widehat{B}} \right) + \vec{S}_T \right]
    \f]
  */
  for(int cell = 0; cell < m_n_cells ; cell++)
  {    
    t_star.get_cell_temperature(cell,m_t_star_vec);
    
    m_materials.calculate_local_temp_and_position(cell, m_t_star_vec);
    
    m_matrix_creator->calculate_k_t(cell, m_t_star_vec, m_k_vec);
    
    k_t.set_kt(cell, m_stage, m_k_vec);    
  }
  
  m_time_data_set = false;
  
  return;
}
