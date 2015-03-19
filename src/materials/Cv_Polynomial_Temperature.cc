/** @file   Scattering_Opacity_Rational.cc
  *   @author pmaginot
  *   @brief Implement the Scattering_Opacity_Rational class
  *   \f$ \sigma_s = \frac{c_1}{c_2 + T^P} \f$
*/
#include "Cv_Polynomial_Temperature.h"

Cv_Polynomial_Temperature::Cv_Polynomial_Temperature(
  const Input_Reader& input_reader, const int mat_num) 
  :
  m_high_poly_degree( input_reader.get_cv_poly_power(mat_num) )
{  
  input_reader.get_cv_poly_coefficients(mat_num , m_poly_coeff);
}

double  Cv_Polynomial_Temperature::get_cv(
  const double position, const double temperature) 
{
  double val = 0.;
  
  double t=1.;
  for(int p=0 ; p <= m_high_poly_degree ; p++)
  {  
    val += m_poly_coeff[p]*t;
    t *= temperature;
  }
  
  return val;
}
