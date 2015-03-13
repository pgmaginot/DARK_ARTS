/** @file   MIP_Left_Boundary_Reflective.cc
  *   @author pmaginot
  *   @brief boundary conditions base class
 */

#include "MIP_Left_Boundary_Reflective.h"

MIP_Left_Boundary_Reflective::MIP_Left_Boundary_Reflective(const Fem_Quadrature& fem_quadrature)
:
  V_MIP_Left_Boundary(fem_quadrature),
    
  m_Rt_L(m_R.transpose() * m_L),
  m_Rt_R(m_R.transpose() * m_R),
  m_Rst_L(m_Rs.transpose() * m_L),
  m_Rst_R(m_Rs.transpose() * m_R),
  m_Rt_Rs(m_R.transpose() * m_Rs),
  m_Rt_Ls(m_R.transpose() * m_Ls)
{

}

void MIP_Left_Boundary_Reflective::add_left_boundary_contributions(const double kappa_12, const double kappa_32 ,   
    const double d_1_l , const double d_1_r , const double d_2_l,
    const double dx_1, const double dx_2, 
    Eigen::MatrixXd& cell_c, Eigen::MatrixXd& cell_cp1, Eigen::VectorXd& rhs) 
{  
  cell_c += kappa_32*m_Rt_R - d_1_r/dx_1*m_Rst_R - d_1_r/dx_1*m_Rt_Rs;
  
  cell_cp1 +=  -kappa_32*m_Rt_L + d_1_r/dx_1*m_Rst_L - d_2_l/dx_2*m_Rt_Ls;

  rhs += m_incoming_current*m_L.transpose();
  return;  
}


