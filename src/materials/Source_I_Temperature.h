#ifndef Source_I_Temperature_h
#define Source_I_Temperature_h

#include <vector>
#include <stdlib.h>
#include "VSource_I.h"
#include "Input_Reader.h"

class Source_I_Temperature : public VSource_I
{
public:
  Source_I_Temperature(const Input_Reader& input_reader, const int mat_num, const double a, const double c, const double sn_sum);
  virtual ~Source_I_Temperature(){}

  double get_intensity_source(const double position, 
    const int group, const int dir, const double time) override;
private:
  const double m_t_start;
  const double m_t_end;
  const double m_output;
};

#endif