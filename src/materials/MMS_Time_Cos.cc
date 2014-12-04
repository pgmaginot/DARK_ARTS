#include "MMS_Time_Cos.h"

MMS_Time_Cos::MMS_Time_Cos(const std::vector<double>& coeff) 
:
  V_MMS_Time(),
  m_cos_coeff(coeff)
{
  
}

double MMS_Time_Cos::get_time_component(const double time)
{
  m_val = m_cos_coeff[0]*cos(time*m_pi/m_cos_coeff[1] + m_cos_coeff[2] ) + m_cos_coeff[3];
  return m_val;
}

double MMS_Time_Cos::get_time_derivative(const double time)
{
  m_val = -m_cos_coeff[0]*sin(time*m_pi/m_cos_coeff[1] + m_cos_coeff[2] )*m_pi/m_cos_coeff[1];
  return m_val;
}