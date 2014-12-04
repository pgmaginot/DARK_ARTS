#ifndef V_MMS_Space_h
#define V_MMS_Space_h

#include <memory>

class V_MMS_Space 
{

public:
  V_MMS_Space(){}
  virtual ~V_MMS_Space(){}

  virtual double get_position_component(const double position) = 0;
  virtual double get_position_derivative(const double position) = 0;
};

#endif