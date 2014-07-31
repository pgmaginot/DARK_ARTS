#ifndef VAbsorption_Opacity_h
#define VAbsorption_Opacity_h

#include <vector>
#include <stdlib.h>

class VAbsorption_Opacity
{

public:
  VAbsorption_Opacity();
  virtual ~VAbsorption_Opacity();

  virtual double get_absorption_opacity(const int group, 
    const double temperature, const double position) = 0;
};

#endif