/** @file   V_MIP_Right_Boundary.cc
  *   @author pmaginot
  *   @brief boundary conditions base class
 */

#include "V_MIP_Right_Boundary.h"

V_MIP_Right_Boundary::V_MIP_Right_Boundary(const Fem_Quadrature& fem_quadrature)
:
m_time(-1.),
m_incoming_current(-1.),
m_L(fem_quadrature.get_dfem_at_left_edge() ),
m_R(fem_quadrature.get_dfem_at_right_edge() ),
m_Ls(fem_quadrature.get_dfem_deriv_at_left_edge() ),
m_Rs(fem_quadrature.get_dfem_deriv_at_right_edge() )
{

}

void V_MIP_Right_Boundary::set_time(const double time)
{
  m_time = time;
  return;
}

void V_MIP_Right_Boundary::set_incoming_current(const double incoming_current)
{
  m_incoming_current = incoming_current;
  return;
}

