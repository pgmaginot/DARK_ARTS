/** @file   MIP_Left_Boundary_Incident_Flux.cc
  *   @author pmaginot
  *   @brief boundary conditions base class
 */

#include "MIP_Left_Boundary_Incident_Flux.h"

MIP_Left_Boundary_Incident_Flux::MIP_Left_Boundary_Incident_Flux(const Fem_Quadrature& fem_quadrature)
:
  V_MIP_Left_Boundary(fem_quadrature),
   
  m_Lt_L(m_L.transpose() * m_L),
  m_Lst_L(m_Ls.transpose() * m_L),
  m_Lt_Ls(m_L.transpose() * m_Ls),
    
  m_Rt_L(m_R.transpose() * m_L),
  m_Rt_R(m_R.transpose() * m_R),
  m_Rst_L(m_Rs.transpose() * m_L),
  m_Rst_R(m_Rs.transpose() * m_R),
  m_Rt_Rs(m_R.transpose() * m_Rs),
  m_Rt_Ls(m_R.transpose() * m_Ls)
{

}

void MIP_Left_Boundary_Incident_Flux::add_left_boundary_contributions(const double kappa_12, const double kappa_32 ,   
    const double d_1_l , const double d_1_r , const double d_2_l,
    const double dx_1, const double dx_2, 
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& cell_c, 
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& cell_cp1) 
{
  cell_c += kappa_32*m_Rt_R - d_1_r/dx_1*m_Rst_R - d_1_r/dx_1*m_Rt_Rs + kappa_12*m_Lt_L + d_1_l/dx_1*m_Lst_L + d_1_l/dx_1*m_Lt_Ls;
  
  cell_cp1 += -kappa_32*m_Rt_L + d_1_r/dx_1*m_Rst_L - d_2_l/dx_2*m_Rt_Ls;
  
  return; 
}

void MIP_Left_Boundary_Incident_Flux::add_left_boundary_rhs_contributions(Eigen::VectorXd& rhs)
{
  return;
}
