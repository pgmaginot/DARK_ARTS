/** @file   V_Diffusion_Matrix_Creator.cc
  *   @author pmaginot
  *   @brief Implement the (partially) virtual Diffusion_Matrix_Cretator class, 
  *     The goal of this clas is to construct MIP diffusion operator matrices in a generic format
  * concrete instantiations of this class will evaluate the following :
  *   \f$ \mathbf{R}_{\widetilde{\Sigma}_a} \f$, \f$ \mathbf{R}_{\widetilde{\Sigma}_s} \f$, 
  *  and $\widetilde{D}$ at cell edges and quadratureintegration points
 */

#include "V_Diffusion_Matrix_Creator.h"

V_Diffusion_Matrix_Creator::V_Diffusion_Matrix_Creator(const Fem_Quadrature& fem_quadrature, Materials& materials,
  const Angular_Quadrature& angular_quadrature,  const Temperature_Data& t_star, const int n_cells)
:
  m_np(fem_quadrature.get_number_of_interpolation_points()), 
  
  m_n_integration_pts(fem_quadrature.get_number_of_integration_points() ),
  m_n_cells(n_cells),
  m_r_sig_s(Eigen::MatrixXd::Zero(m_np,m_np) ),
  
  m_r_sig_a(Eigen::MatrixXd::Zero(m_np,m_np) ), 
  
  m_d_r_cm1(-1) ,
  m_d_l_c(-1) ,
  m_d_r_c(-1) ,
  m_d_l_cp1(-1) ,
  
  m_dfem_derivatives(m_n_integration_pts*m_np,0.),
  m_integration_wts(m_n_integration_pts,0.),
  m_d_at_integration_pts(m_n_integration_pts,0.),
  
  m_t_eval(t_star),
  m_materials(materials),
  m_angular_quadrature(angular_quadrature),
  m_t_eval_vec(Eigen::VectorXd::Zero(m_np)),
  m_dx(-1.),
  m_cell_num(-1)  ,
  m_group_num(-1),
  
  m_mip_assembler(fem_quadrature)
{    
  MATRIX_INTEGRATION matrix_type = fem_quadrature.get_integration_type() ;
  /// initialize matrix constructor
  if(matrix_type ==  EXACT)
  {
    m_mtrx_builder = std::make_shared<Matrix_Construction_Exact>(fem_quadrature,materials) ;
  }
  else if(matrix_type == TRAD_LUMPING)
  {
    m_mtrx_builder = std::make_shared<Matrix_Construction_Trad_Lumping>(fem_quadrature,materials) ;
  }
  else if(matrix_type == SELF_LUMPING)
  {
    m_mtrx_builder = std::make_shared<Matrix_Construction_Self_Lumping>(fem_quadrature,materials) ;
  }
  
  fem_quadrature.get_dfem_derivatives_at_integration_points(m_dfem_derivatives);
  fem_quadrature.get_integration_weights(m_integration_wts);
}

void V_Diffusion_Matrix_Creator::set_cell_group_information( const int cell, const int group)
{
  m_cell_num = cell;
  m_group_num = group;
  
  /// get \f$ \vec{T} \f$ for material evaluation
  m_t_eval.get_cell_temperature(m_cell_num,m_t_eval_vec) ;
  
  m_materials.calculate_local_temp_and_position(cell,m_t_eval_vec);
    
  /// set cell width
  m_dx = m_materials.get_cell_width();
  return;
}

void V_Diffusion_Matrix_Creator::calculate_d_dependent_quantities(
  double& d_r_cm1, double& d_l_c , double& d_r_c , double& d_l_cp1, Eigen::MatrixXd& s_matrix)
{
  evaluate_diffusion_coefficents(d_r_cm1, d_l_c , d_r_c , d_l_cp1);
  calculate_s_matrix(s_matrix);
  
  return;
}

void V_Diffusion_Matrix_Creator::calculate_s_matrix(Eigen::MatrixXd& s_matrix)
{
  s_matrix = Eigen::MatrixXd::Zero(m_np,m_np);
  for(int i = 0 ; i < m_np ; i++)
  {
    for(int j = 0 ; j < m_np ; j++)
    {
      for(int q = 0 ; q < m_n_integration_pts ; q++)
      {
        // std::cout << "Inside the loop\n";
        s_matrix(i,j) += m_integration_wts[q]*m_d_at_integration_pts[q]*
          m_dfem_derivatives[q+i*m_n_integration_pts]*m_dfem_derivatives[q+j*m_n_integration_pts];
      }      
    }
  }
  s_matrix *= 2./m_dx;

}
