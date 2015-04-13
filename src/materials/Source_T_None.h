#ifndef Source_T_None_h
#define Source_T_None_h

#include <vector>
#include <stdlib.h>
#include "VSource_T.h"
#include "Input_Reader.h"

class Source_T_None : public VSource_T
{

public:
  Source_T_None(const Input_Reader& input_reader, const int mat_num);
  virtual ~Source_T_None(){}

  double get_temperature_source(const double position, const double time) override;
private:

};

#endif