/** @file   Absorption_Opacity_Ideal.cc
  *   @author pmaginot
  *   @brief Implement the concrete Absorption_Opacity_Ideal class
  *   \f$ \sigma_a = \frac{1}{T^3} \f$ of the abstract class VAbsorption_Opacity
*/
#include "Absorption_Opacity_Ideal.h"

Absorption_Opacity_Ideal::Absorption_Opacity_Ideal()
{

}

Absorption_Opacity_Ideal::~Absorption_Opacity_Ideal()
{

}

double Absorption_Opacity_Ideal::get_absorption_opacity(const int group, const double temperature, const double position)
{
  return 1./( pow(temperature,3));
}

