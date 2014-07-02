/** @file   Angular_Quadrature.cc
  *   @author pmaginot
  *   @brief Implement the Angular_Quadrature class, stores discrete ordinates, weights, Legendre polynomials
 */

#include "Angular_Quadrature.h"

Angular_Quadrature::Angular_Quadrature(const Input_Reader& input_reader, const Quadrule_New& quad_fun)
{
  /// get number of groups, number of angles, quadrature type from Input_Reader
  m_n_dir = input_reader.get_number_of_angles();
  m_n_groups = input_reader.get_number_of_groups();
  
  ANGULAR_QUADRATURE_TYPE quad_type = input_reader.get_angular_quadrature_type();
  
  m_n_legendre_moments = input_reader.get_number_of_legendre_moments();
  
  m_mu.resize(m_n_dir,0.);
  m_w.resize(m_n_dir,0.);
  
  const int n_leg_evals = m_n_legendre_moments * m_n_dir;
  m_legendre_poly.resize(n_leg_evals,0.);
  
  if(quad_type == GAUSS_ANGLE)
  {
  
  }
  else if(quad_type == LOBATTO_ANGLE)
  {
  
  }
}
    
int Angular_Quadrature::get_number_of_dir(void) const
{
  return m_n_dir;
}

int Angular_Quadrature::get_number_of_groups(void) const{
  return m_n_groups;
}

int Angular_Quadrature::get_number_of_leg_moments(void) const
{
  return m_n_legendre_moments;
}

