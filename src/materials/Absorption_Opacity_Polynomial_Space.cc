/** @file   Scattering_Opacity_Rational.cc
  *   @author pmaginot
  *   @brief Implement the Scattering_Opacity_Rational class
  *   \f$ \sigma_s = \frac{c_1}{c_2 + T^P} \f$
*/
#include "Absorption_Opacity_Polynomial_Space.h"

Absorption_Opacity_Polynomial_Space::Absorption_Opacity_Polynomial_Space(
  const Input_Reader& input_reader, const int mat_num) 
  :
  m_high_poly_degree{ input_reader.get_abs_int_constant(mat_num) }
{  
  input_reader.get_absorption_poly_coeff(mat_num , m_poly_coeff);
}

double  Absorption_Opacity_Polynomial_Space::get_absorption_opacity(const int group, 
  const double temperature, const double position)
{
  double val = 0.;
  
  for(int p=0 ; p <= m_high_poly_degree ; p++)
    val += m_poly_coeff[p]*pow(position , p);
  
  return val;
}
