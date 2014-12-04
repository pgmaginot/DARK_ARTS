#include "MMS_Time_Poly.h"

MMS_Time_Poly::MMS_Time_Poly(const std::vector<double>& coeff) 
:
  V_MMS_Time(),
  m_poly_coeff(coeff),
  m_max_poly_degree_p1( m_poly_coeff.size() )
{
  
}

double MMS_Time_Poly::get_time_component(const double time)
{
  m_val = 0.;
  m_val += m_poly_coeff[0];
  m_pow = time;
  for(int i = 1; i < m_max_poly_degree_p1 ; i++)
  {
    m_val += m_poly_coeff[i]*m_pow;
    m_pow *= time;
  }  
  return m_val;
}

double MMS_Time_Poly::get_time_derivative(const double time)
{
  m_val = 0.;
  m_pow = 1.;
  for(int i = 1; i < m_max_poly_degree_p1 ; i++)
  {
    m_val += double(i)*m_poly_coeff[i]*m_pow;
    m_pow *= time;
  }  
  return m_val;
}