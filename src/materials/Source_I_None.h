#ifndef Source_I_None_h
#define Source_I_None_h

#include <vector>
#include <stdlib.h>
#include "VSource_I.h"
#include "Input_Reader.h"

class Source_I_None : public VSource_I
{
public:
  Source_I_None(const Input_Reader& input_reader, const int mat_num);
  virtual ~Source_I_None(){}

  double get_intensity_source(const double position, 
    const int group, const int dir, const double time) override;
private:
  
};

#endif