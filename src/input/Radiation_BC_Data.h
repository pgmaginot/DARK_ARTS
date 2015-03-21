#ifndef Radiation_BC_Data_h
#define Radiation_BC_Data_h

#include "Inputs_Allowed.h"

struct Radiation_BC_Data{
  RADIATION_BC_TYPE type;
  double value;
  INCIDENT_BC_VALUE_TYPE value_type; 
  BC_ANGLE_DEPENDENCE angle_dependence;
  BC_ENERGY_DEPENDENCE energy_dependence;
  BC_TIME_DEPENDENCE time_dependence;
  double start_time;
  double end_time;
  
  /// Default to every value be invalid
  Radiation_BC_Data() :
    type(INVALID_RADIATION_BC_TYPE),
    value(-1.) ,
    value_type(INVALID_INCIDENT_BC_VALUE_TYPE), 
    angle_dependence(INVALID_BC_ANGLE_DEPENDENCE),
    energy_dependence(INVALID_BC_ENERGY_DEPENDENCE),
    time_dependence(INVALID_BC_TIME_DEPENDENCE),
    start_time(-1.),
    end_time(-2.)
  {}

};

#endif