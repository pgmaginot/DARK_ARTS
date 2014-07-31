#ifndef VSource_T_h
#define VSource_T_h

#include <vector>
#include <stdlib.h>

class VSource_T
{

public:
  VSource_T();
  virtual ~VSource_T();

  virtual double get_temperature_source(const double position, const double time) = 0;
};

#endif