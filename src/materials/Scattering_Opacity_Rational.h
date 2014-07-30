#ifndef Scattering_Opacity_Rational_h
#define Scattering_Opacity_Rational_h

#include "VScattering_Opacity.h"
#include <vector>
#include <stdlib.h>
#include <math.h>

class Scattering_Opacity_Rational: public VScattering_Opacity
{

public:
  Scattering_Opacity_Rational();
  virtual ~Scattering_Opacity_Rational();

  double get_scattering_opacity(const int l_mom, const int group, const double temperature, const double position) override;
};

#endif