/** @file   Local_MIP_Assembler.cc
  *   @author pmaginot
  *   @brief Implement the (partially) virtual Diffusion_Matrix_Cretator class, 
  *     The goal of this clas is to construct MIP diffusion operator matrices in a generic format
  * concrete instantiations of this class will evaluate the following :
  *   \f$ \mathbf{R}_{\widetilde{\Sigma}_a} \f$, \f$ \mathbf{R}_{\widetilde{\Sigma}_s} \f$, 
  *  and $\widetilde{D}$ at cell edges and quadratureintegration points
 */

#include "Local_MIP_Assembler.h"

Local_MIP_Assembler::Local_MIP_Assembler(const Fem_Quadrature& fem_quadrature, const Input_Reader& input_reader)
:
  m_L( fem_quadrature.get_dfem_at_left_edge() ),
  m_Ls( fem_quadrature.get_dfem_deriv_at_left_edge() ),
  m_R( fem_quadrature.get_dfem_at_right_edge() ),
  m_Rs( fem_quadrature.get_dfem_deriv_at_right_edge() ),
  
  m_Lt_R( m_L.transpose() * m_R ),
  m_Lst_R(m_Ls.transpose() * m_R ),
  m_Lt_Rs(m_L.transpose() * m_Rs),
  
  m_Lt_L(m_L.transpose() * m_L),
  m_Rt_R(m_R.transpose() * m_R),
  m_Lst_L(m_Ls.transpose() * m_L),
  m_Rst_R(m_Rs.transpose() * m_R),
  m_Lt_Ls(m_L.transpose() * m_Ls),
  m_Rt_Rs(m_R.transpose() * m_Rs),
  
  m_Rt_L(m_R.transpose() * m_L),
  m_Rst_L(m_Rs.transpose() * m_L),
  m_Rt_Ls(m_R.transpose() * m_Ls)
  {
    std::cout << "m_L: " << m_L << std::endl;
    std::cout << "m_R: " << m_R << std::endl;
    std::cout << "m_Ls: " << m_Ls << std::endl;
    std::cout << "m_Rs: " << m_Rs << std::endl;
    
    /// set MIP boundary condition functions
    if( input_reader.get_radiation_bc_type_left() == REFLECTIVE_BC)
    {
      ///adds C_{MIP} \vec{L}^T to the RHS . C_{MIP} is going to need to come from somewhere else (Transport_Sweep)
      m_left_boundary = std::make_shared<MIP_Left_Boundary_Reflective>(fem_quadrature);
    }
    else
    {
      m_left_boundary = std::make_shared<MIP_Left_Boundary_Incident_Flux>(fem_quadrature);
    }
    
     m_right_boundary = std::make_shared<MIP_Right_Boundary_Incident_Flux>(fem_quadrature);
  }
  
void Local_MIP_Assembler::calculate_left_boundary_matrices(const double kappa_12, const double kappa_32 ,   
  const double dx_c, const double dx_cp1, 
  const double d_c_l , const double d_c_r , const double d_cp1_l,
  const Eigen::MatrixXd& r_sig_a, const Eigen::MatrixXd& s_matrix,
  Eigen::MatrixXd& cell_c, Eigen::MatrixXd& cell_cp1)
{  
  cell_c = r_sig_a + s_matrix;
  
  m_left_boundary->add_left_boundary_contributions(kappa_12,kappa_32, d_c_l, d_c_r , d_cp1_l, dx_c, dx_cp1 , cell_c , cell_cp1);  
  return;
}

void Local_MIP_Assembler::additive_update_left_boundary_rhs(Eigen::VectorXd& rhs)
{
  m_left_boundary->add_left_boundary_rhs_contributions(rhs);
}
    
void Local_MIP_Assembler::calculate_interior_matrices(const double kappa_cm12, const double kappa_cp12 , 
  const double dx_cm1, const double dx_c, const double dx_cp1, 
  const double d_cm1_r, const double d_c_l , const double d_c_r , const double d_cp1_l,
  const Eigen::MatrixXd& r_sig_a, const Eigen::MatrixXd& s_matrix,
  Eigen::MatrixXd& cell_cm1, Eigen::MatrixXd& cell_c, Eigen::MatrixXd& cell_cp1)
{
  /**
    \f[
      -\kappa_{c-1/2}\vec{L}^T\vec{R} - \frac{D(x_{c-1/2}^+)}{\Delta x_c} \vec{L}_s^T \vec{R} + \frac{ D(x_{c-1/2}^-) }{\Delta x_{c-1}} \vec{L}^T \vec{R}_s
    \f]
  */
  cell_cm1 = -1.*kappa_cm12*m_Lt_R - d_c_l/dx_c*m_Lst_R + d_cm1_r/dx_cm1*m_Lt_Rs;

  cell_c = r_sig_a + s_matrix;
  cell_c += kappa_cm12*m_Lt_L + kappa_cp12*m_Rt_R + d_c_l/dx_c*m_Lst_L - d_c_r/dx_c*m_Rst_R
    + d_c_l/dx_c*m_Lt_Ls - d_c_r/dx_c*m_Rt_Rs;
    
  cell_cp1 = -1.*kappa_cp12*m_Rt_L + d_c_r/dx_c*m_Rst_L - d_cp1_l/dx_cp1*m_Rt_Ls;
  return;
}
    
void Local_MIP_Assembler::calculate_right_boundary_matrices(const double kappa_nm12, const double kappa_np12 ,   
  const double dx_cm1, const double dx_c, 
  const double d_cm1_r, const double d_c_l , const double d_c_r , 
  const Eigen::MatrixXd& r_sig_a, const Eigen::MatrixXd& s_matrix,
  Eigen::MatrixXd& cell_cm1, Eigen::MatrixXd& cell_c)
{
  cell_c = r_sig_a + s_matrix;
  
  m_right_boundary->add_right_boundary_contributions(kappa_nm12, kappa_np12 ,d_cm1_r, d_c_l ,d_c_r , 
    dx_cm1, dx_c, cell_cm1, cell_c);
    
  return;
}    

void Local_MIP_Assembler::additive_update_right_boundary_rhs(Eigen::VectorXd& rhs)
{
  m_right_boundary->add_right_boundary_rhs_contributions(rhs);
}

