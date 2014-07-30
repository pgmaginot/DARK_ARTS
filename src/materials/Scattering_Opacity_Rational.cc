/** @file   Scattering_Opacity_Rational.cc
  *   @author pmaginot
  *   @brief Implement the Scattering_Opacity_Rational class
  *   \f$ \sigma_s = \frac{c_1}{c_2 + T^P} \f$
*/
#include "Scattering_Opacity_Rational.h"

Scattering_Opacity_Rational::Scattering_Opacity_Rational(){}

Scattering_Opacity_Rational::~Scattering_Opacity_Rational(){}

double  Scattering_Opacity_Rational::get_scattering_opacity(const int l_mom, const int group, 
  const double temperature, const double position)
{
  return 0.;
}
