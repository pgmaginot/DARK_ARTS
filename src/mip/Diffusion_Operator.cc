/** @file   Diffusion_Operator.cc
  *   @author pmaginot
  *   @brief Implement a MIP DSA operator
 */

#include "Diffusion_Operator.h"

Diffusion_Operator::Diffusion_Operator(const Input_Reader& input_reader, const Fem_Quadrature& fem_quadrature, 
  const Cell_Data& cell_data, Materials& materials, const Angular_Quadrature& angular_quadrature,
  const int n_groups, const Temperature_Data& t_eval, const bool is_wg_solve)
:
  m_sn_w( angular_quadrature.get_sum_w() ),
  m_np(fem_quadrature.get_number_of_interpolation_points()),
  
  m_r_sig_s(Eigen::MatrixXd::Zero(m_np,m_np) ),
  m_r_sig_a(Eigen::MatrixXd::Zero(m_np,m_np) ), 
  
  m_no_mu_pos_l_matrix(Eigen::MatrixXd::Zero(m_np,m_np) ),
  m_no_mu_neg_l_matrix(Eigen::MatrixXd::Zero(m_np,m_np) ),
  
  m_sdirk_a_stage(-1.),
  m_dt(-1.),  
  m_time_stage(-1.)
{    
  if(is_wg_solve)
  {
    if(n_groups ==1)
    {
      m_diffusion_matrix_creator = std::make_shared<Diffusion_Matrix_Creator_Grey>
        (fem_quadrature,materials,t_eval);
    }
    else
    {
      throw Dark_Arts_Exception(TIME_MARCHER , "MF within group scattering MIP not coded");
    }  
  }
  else
  {
    LMFGA_STRUCTURE lmfga_type = input_reader.get_lmfga_structure();
    if( lmfga_type == GROUP_COLLAPSE)
    {
      throw Dark_Arts_Exception(TIME_MARCHER , "MF LMFGA with Group Collapse Not Coded");
    }
    else if( lmfga_type == NO_COLLAPSE)
    {
      throw Dark_Arts_Exception(TIME_MARCHER , "MF LMFGA without Group Collapse Not Coded");
    }  
  }
  
  if( !m_diffusion_matrix_creator )
    throw Dark_Arts_Exception(TRANSPORT, "Unable to allocate a Diffusion_Matrix_Creator in Diffusion_Operator");
}

bool Diffusion_Operator::check_all_eigen_variables_for_finite(void)
{
  bool bad_eigen_variables = false;
  std::stringstream err;
  err << std::endl;
  
  // if(!m_r_sig_t.allFinite())
  // {
    // bad_eigen_variables = true;
    // err << "m_r_sig_t has non finite values!\n";
  // }
  
  std::cout << err.str();
    
  return bad_eigen_variables;
}
