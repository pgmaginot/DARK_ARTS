#include "Transport_BC_MMS.h"

Transport_BC_MMS::Transport_BC_MMS(const Angular_Quadrature& angular_quadrature,
  const Input_Reader& input_reader, const double x_b)
:
  V_Transport_BC(),  
  m_x_boundary(x_b),
  m_intensity(input_reader,angular_quadrature)
{
  
}

double Transport_BC_MMS::get_boundary_condition(const double mu, const int grp , const double time) 
{
  return m_intensity.get_mms_intensity(m_x_boundary , time , mu);
}

