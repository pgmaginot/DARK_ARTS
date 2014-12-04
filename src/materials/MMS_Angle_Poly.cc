#include "MMS_Angle_Poly.h"

MMS_Angle_Poly::MMS_Angle_Poly(const std::vector<double>& poly_coeff, const Angular_Quadrature& angular_quadrature) 
:
  V_MMS_Angle(),
  m_poly_coeff(poly_coeff),
  m_max_poly_degree_p1(m_poly_coeff.size() )
{
  angular_quadrature.get_all_mu(m_angles);
}

double MMS_Angle_Poly::get_angle_component(const int dir)
{ 
  m_val = m_poly_coeff[0];
  m_mu = m_angles[dir];
  m_pow = m_mu;
  for(int i = 1; i < m_max_poly_degree_p1 ; i++)
  {
    m_val += m_poly_coeff[i]*m_pow;
    m_pow *= m_mu;
  }  
  return m_val;
}

double MMS_Angle_Poly::get_angle_component(const double mu)
{ 
  m_val = m_poly_coeff[0];
  m_mu = mu;
  m_pow = mu;
  for(int i = 1; i < m_max_poly_degree_p1 ; i++)
  {
    m_val += m_poly_coeff[i]*m_pow;
    m_pow *= m_mu;
  }  
  return m_val;
}