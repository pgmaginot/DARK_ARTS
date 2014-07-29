/** @file   Scattering_Opacity_Ideal.cc
  *   @author pmaginot
  *   @brief Implement the Scattering_Opacity_Ideal class
  *   No scattering (Marshak wave problem)
*/
#include "Scattering_Opacity_Ideal.h"

Scattering_Opacity_Ideal::Scattering_Opacity_Ideal(){}

Scattering_Opacity_Ideal::~Scattering_Opacity_Ideal(){}

double  Scattering_Opacity_Ideal::get_scattering_opacity(const int l_mom, const int group, 
  const double temperature, const double position)
{
  return 0.;
}
