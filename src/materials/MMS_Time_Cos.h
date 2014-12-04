#ifndef MMS_Time_Cos_h
#define MMS_Time_Cos_h

#include <vector>
#include <stdlib.h>
#include <math.h>   
#include "V_MMS_Time.h"

class MMS_Time_Cos : public V_MMS_Time
{
public:
  MMS_Time_Cos(const std::vector<double>& coeff);
  virtual ~MMS_Time_Cos(){}
  
  double get_time_component(const double time) override;
  double get_time_derivative(const double time) override;
private:
  const std::vector<double> m_cos_coeff;
  double m_val;
  const double m_pi = 3.14159265358979323846264;
};

#endif