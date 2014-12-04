#ifndef MMS_Space_Cos_h
#define MMS_Space_Cos_h

#include <vector>
#include <stdlib.h>
#include <math.h>   
#include "V_MMS_Space.h"

class MMS_Space_Cos : public V_MMS_Space
{
public:
  MMS_Space_Cos(const std::vector<double>& coeff);
  virtual ~MMS_Space_Cos(){}
  
  double get_position_component(const double position) override;
  double get_position_derivative(const double position) override;
private:
  const std::vector<double> m_cos_coeff;
  double m_val;
  const double m_pi = 3.14159265358979323846264;
};

#endif