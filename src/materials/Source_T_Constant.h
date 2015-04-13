#ifndef Source_T_Constant_h
#define Source_T_Constant_h

#include <vector>
#include <stdlib.h>
#include "VSource_T.h"
#include "Input_Reader.h"

class Source_T_Constant : public VSource_T
{

public:
  Source_T_Constant(const Input_Reader& input_reader, const int mat_num);
  virtual ~Source_T_Constant(){}

  double get_temperature_source(const double position, const double time) override;
private:
  const double m_t_start;
  const double m_t_end;
  const double m_output;
};

#endif