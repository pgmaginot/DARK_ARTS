/** @file   Sweep_Matrix_Creator_Grey.cc
  *   @author pmaginot
  *   @brief Implement the Sweep_Matrix_Creator_Grey class (make pseudo matrices and sources for transport sweep)
  *    Creates the "funny" matrices that arise in the Planck/temperature linearization when accounting for spatially varying material properties
 */

#include "Sweep_Matrix_Creator_Grey.h"

Sweep_Matrix_Creator_Grey::Sweep_Matrix_Creator_Grey(const Fem_Quadrature& fem_quadrature)
:
  V_Sweep_Matrix_Creator( fem_quadrature )
{  
  
}
void Sweep_Matrix_Creator_Grey::construct_r_sig_t(const int cell, const int grp, Eigen::MatrixXd& r_sig_t)
{
  return;
}

void Sweep_Matrix_Creator_Grey::construct_r_sig_s(const int cell, const int grp, const int l_mom, Eigen::MatrixXd& r_sig_s)
{
  return;
}

void Sweep_Matrix_Creator_Grey::construct_s_i(const int cell,const int grp, const int l_mom, Eigen::VectorXd& s_i)
{
  return;
}
