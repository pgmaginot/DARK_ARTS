/** @file   MIP_Left_Boundary_Reflective.cc
  *   @author pmaginot
  *   @brief boundary conditions base class
 */

#include "MIP_Left_Boundary_Reflective.h"

MIP_Left_Boundary_Reflective::MIP_Left_Boundary_Reflective(const Angular_Quadrature& angular_quadrature)
:
V_MIP_Left_Boundary(angular_quadrature)
{

}

void MIP_Left_Boundary_Reflective::add_left_boundary_contributions(
  Eigen::MatrixXd& cell_c, Eigen::MatrixXd& cell_cp1, Eigen::VectorXd& rhs) 
{
  return;
  
}


