#include "MMS_Space_Poly.h"

MMS_Space_Poly::MMS_Space_Poly(const std::vector<double>& coeff) 
:
  V_MMS_Space(),
  m_poly_coeff(coeff),
  m_max_poly_degree_p1( m_poly_coeff.size() )
{
  
}

double MMS_Space_Poly::get_position_component(const double position)
{
  m_val = m_poly_coeff[0];
  m_pow = position;
  for(int i = 1; i < m_max_poly_degree_p1 ; i++)
  {
    m_val += m_poly_coeff[i]*m_pow;
    m_pow *= position;
  }  
  return m_val;
}

double MMS_Space_Poly::get_position_derivative(const double position)
{
  m_val = 0.;
  m_pow = 1.;
  for(int i = 1; i < m_max_poly_degree_p1 ; i++)
  {
    m_val += double(i)*m_poly_coeff[i]*m_pow;
    m_pow *= position;
  }  
  return m_val;
}