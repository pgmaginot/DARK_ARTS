#ifndef VCv_h
#define VCv_h

#include <vector>
#include <stdlib.h>

class VCv
{

public:
  VCv();
  virtual ~VCv();

  virtual double get_cv(const double position, const double temperature) = 0;
};

#endif