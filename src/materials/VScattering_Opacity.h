#ifndef VScattering_Opacity_h
#define VScattering_Opacity_h

#include <vector>
#include <stdlib.h>

class VScattering_Opacity
{

public:
  VScattering_Opacity();
  virtual ~VScattering_Opacity();

  virtual double get_scattering_opacity(const int l_mom, const int group, const double temperature, const double position) = 0;
};

#endif