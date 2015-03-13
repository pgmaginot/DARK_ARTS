/** @file   MIP_Right_Boundary_Incident_Flux.cc
  *   @author pmaginot
  *   @brief boundary conditions base class
 */

#include "MIP_Right_Boundary_Incident_Flux.h"

MIP_Right_Boundary_Incident_Flux::MIP_Right_Boundary_Incident_Flux(const Fem_Quadrature& fem_quadrature)
:
V_MIP_Right_Boundary(fem_quadrature),
m_Lt_R(m_L.transpose() * m_R),
m_Lt_L(m_L.transpose() * m_L),

m_Lst_R(m_Ls.transpose() * m_R),
m_Lst_L(m_Ls.transpose() * m_L),

m_Lt_Rs(m_L.transpose() * m_Rs),
m_Lt_Ls(m_L.transpose() * m_Ls),

m_Rt_R(m_R.transpose() * m_R),
m_Rst_R(m_Rs.transpose() * m_R),
m_Rt_Rs(m_R.transpose() * m_Rs)
{

}

void MIP_Right_Boundary_Incident_Flux::add_right_boundary_contributions(const double kappa_nm12, const double kappa_np12 ,   
    const double d_cm1_r, const double d_c_l , const double d_c_r , 
    const double dx_cm1, const double dx_c, 
    Eigen::MatrixXd& cell_cm1, Eigen::MatrixXd& cell_c, Eigen::VectorXd& rhs)
{
  cell_cm1 += -kappa_nm12*m_Lt_R - d_c_l/dx_c*m_Lst_R + d_cm1_r/dx_cm1*m_Lt_Rs;
  
  cell_c += kappa_nm12*m_Lt_L + d_c_l/dx_c*m_Lst_L + d_c_l/dx_c*m_Lt_Ls + kappa_np12*m_Rt_R - d_c_r/dx_c*m_Rst_R - d_c_r/dx_c*m_Rt_Rs;
  return;
}
