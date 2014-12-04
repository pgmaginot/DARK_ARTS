#ifndef V_MMS_Angle_h
#define V_MMS_Angle_h

#include <memory>


class V_MMS_Angle 
{

public:
  V_MMS_Angle(){}
  virtual ~V_MMS_Angle(){}

  virtual double get_angle_component(const int dir) = 0;
  virtual double get_angle_component(const double mu) = 0;
};

#endif