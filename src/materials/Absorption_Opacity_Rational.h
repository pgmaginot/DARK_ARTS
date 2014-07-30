#ifndef Absorption_Opacity_Rational_h
#define Absorption_Opacity_Rational_h

#include "VAbsorption_Opacity.h"
#include <vector>
#include <stdlib.h>
#include <math.h>

class Absorption_Opacity_Rational: public VAbsorption_Opacity
{

public:
  Absorption_Opacity_Rational();
  virtual ~Absorption_Opacity_Rational();

  double get_absorption_opacity(const int group, const double temperature, const double position) override;
};

#endif