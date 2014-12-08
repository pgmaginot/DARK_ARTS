#include "MMS_Space_Cos.h"

MMS_Space_Cos::MMS_Space_Cos(const std::vector<double>& coeff) 
:
  V_MMS_Space(),
  m_cos_coeff(coeff)
{
}

double MMS_Space_Cos::get_position_component(const double position)
{
  m_val = m_cos_coeff[0]*cos(position*m_pi/m_cos_coeff[1] + m_cos_coeff[2] ) + m_cos_coeff[3];
  return m_val;
}

double MMS_Space_Cos::get_position_derivative(const double position)
{
  m_val = -m_cos_coeff[0]*sin(position*m_pi/m_cos_coeff[1] + m_cos_coeff[2] )*m_pi/m_cos_coeff[1];
  return m_val;
}
