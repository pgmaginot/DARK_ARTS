/** @file   MIP_Right_Boundary_Incident_Flux.cc
  *   @author pmaginot
  *   @brief boundary conditions base class
 */

#include "MIP_Right_Boundary_Incident_Flux.h"

MIP_Right_Boundary_Incident_Flux::MIP_Right_Boundary_Incident_Flux(const Angular_Quadrature& angular_quadrature)
:
V_MIP_Right_Boundary(angular_quadrature)
{

}

void MIP_Right_Boundary_Incident_Flux::add_right_boundary_contributions(
  Eigen::MatrixXd& cell_c, Eigen::MatrixXd& cell_cp1, Eigen::VectorXd& rhs)
{
  return;
}
