#ifndef V_MMS_Time_h
#define V_MMS_Time_h

#include <memory>

class V_MMS_Time
{

public:
  V_MMS_Time(){}
  virtual ~V_MMS_Time(){}

  virtual double get_time_component(const double time) = 0;
  virtual double get_time_derivative(const double time) = 0;
};

#endif