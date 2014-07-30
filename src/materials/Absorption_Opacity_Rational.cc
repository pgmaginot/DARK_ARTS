/** @file   Absorption_Opacity_Rational.cc
  *   @author pmaginot
  *   @brief Implement the concrete Absorption_Opacity_Rational class
  *   \f$ \sigma_a = \frac{c_1}{c_2 + T^P} \f$
*/
#include "Absorption_Opacity_Rational.h"

Absorption_Opacity_Rational::Absorption_Opacity_Rational()
{

}

Absorption_Opacity_Rational::~Absorption_Opacity_Rational()
{

}

double Absorption_Opacity_Rational::get_absorption_opacity(const int group, const double temperature, const double position)
{
  return 1./( pow(temperature,3));
}

