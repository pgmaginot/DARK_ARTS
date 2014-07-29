#ifndef Absorption_Opacity_Ideal_h
#define Absorption_Opacity_Ideal_h

#include "VAbsorption_Opacity.h"
#include <vector>
#include <stdlib.h>
#include <math.h>

class Absorption_Opacity_Ideal: public VAbsorption_Opacity
{

public:
  Absorption_Opacity_Ideal();
  virtual ~Absorption_Opacity_Ideal();

  double get_absorption_opacity(const int group, const double temperature, const double position) override;
};

#endif