#ifndef VSource_I_h
#define VSource_I_h

#include <vector>
#include <stdlib.h>

class VSource_I
{

public:
  VSource_I();
  virtual ~VSource_I();

  virtual double get_intensity_source(const double position, const int group, const int dir) = 0;
};

#endif