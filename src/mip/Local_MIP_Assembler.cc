/** @file   Local_MIP_Assembler.cc
  *   @author pmaginot
  *   @brief Implement the (partially) virtual Diffusion_Matrix_Cretator class, 
  *     The goal of this clas is to construct MIP diffusion operator matrices in a generic format
  * concrete instantiations of this class will evaluate the following :
  *   \f$ \mathbf{R}_{\widetilde{\Sigma}_a} \f$, \f$ \mathbf{R}_{\widetilde{\Sigma}_s} \f$, 
  *  and $\widetilde{D}$ at cell edges and quadratureintegration points
 */

#include "Local_MIP_Assembler.h"

Local_MIP_Assembler::Local_MIP_Assembler(const Fem_Quadrature& fem_quadrature)
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
    std::cout << "This is m_L:\n" << m_L << std::endl;
  }
  
void Local_MIP_Assembler::calculate_left_boundary_matrices(const double kappa_cm12, const double kappa_cp12 , 
  const double dx_c, const double dx_cp1,
  const double d_c_l , const double d_c_r , const double d_cp1_l,
  const Eigen::MatrixXd& r_sig_a, const Eigen::MatrixXd& s_matrix,
  Eigen::MatrixXd& cell_c, Eigen::MatrixXd& cell_cp1, Eigen::VectorXd& rhs)
{
  
  cell_c = r_sig_a + s_matrix;
  
  return;
}
    
void Local_MIP_Assembler::calculate_interior_matrices(const double kappa_cm12, const double kappa_cp12 , 
  const double dx_cm1, const double dx_c, const double dx_cp1, 
  const double d_cm12_r, const double d_c_l , const double d_c_r , const double d_cp1_l,
  const Eigen::MatrixXd& r_sig_a, const Eigen::MatrixXd& s_matrix,
  Eigen::MatrixXd& cell_cm1, Eigen::MatrixXd& cell_c, Eigen::MatrixXd& cell_cp1)
{
  cell_cm1 = -kappa_cm12*m_Lt*m_R - d_c_l/dx_c*m_Lst*m_R + d_cm1_r/dx_cm1*m_Lt*m_Rs;

  cell_c = r_sig_a + s_matrix;
  cell_c += kappa_cm12*m_Lt*m_L + kappa_cp12*m_Rt*m_R + d_c_l/dx_c*m_Lst*m_L - d_c_r/dx_c*m_Rst*m_R
    + d_c_l/dx_c*m_Lt*m_Ls - d_cp1_l/dx_c*m_Rt*m_Rs;
    
  cell_cp1 = -kappa_cp12*m_Rt*m_L + d_c_r/dx_c*m_Rst*m_L - d_cp1_l/dx_cp1*m_Rt*m_Ls;
  return;
}
    
void Local_MIP_Assembler::calculate_right_boundary_matrices(const double kappa_cm12, const double kappa_cp12 , 
  const double dx_cm1, const double dx_c,
  const double d_cm12_r, const double d_c_l , const double d_c_r ,
  const Eigen::MatrixXd& r_sig_a, const Eigen::MatrixXd& s_matrix,
  Eigen::MatrixXd& cell_cm1, Eigen::MatrixXd& cell_c)
{
  cell_c = r_sig_a + s_matrix;
  return;
}    
