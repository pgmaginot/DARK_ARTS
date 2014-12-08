/** @file   MMS_Temperature.cc
  *   @author pmaginot
  *   @brief A convenience class for accesing a manufactured solution temperature profile (primarily for IC)
*/
#include "MMS_Temperature.h"

MMS_Temperature::MMS_Temperature(const Input_Reader& input_reader)
{
  TEMPERATURE_SPACE_MMS t_space_dependence = input_reader.get_mms_temperature_space_dependence() ; 
  TIME_MMS_TYPE time_dependence = input_reader.get_mms_time_dependence();
  
  std::vector<double> time_coeff;
  input_reader.get_mms_time_coeff(time_coeff);
  
  std::vector<double> t_space_coeff;
  input_reader.get_mms_temperature_space_coeff(t_space_coeff);  
  switch(time_dependence)
  {
    case POLY_TIME:
    {      
      m_time_dep = std::shared_ptr<V_MMS_Time> (new MMS_Time_Poly(time_coeff) );
      break;
    }
    case COS_TIME:
    {
      m_time_dep = std::shared_ptr<V_MMS_Time> (new MMS_Time_Cos(time_coeff) );
      break;
    }
    case INVALID_TIME_MMS_TYPE:
    {
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Trying to create an invalid MMS temperature time dependence");
      break;
    }
  }
  switch(t_space_dependence)
  {
    case TEMP_POLY_SPACE:
    {
      m_temp_space_dep = std::shared_ptr<V_MMS_Space> (new MMS_Space_Poly(t_space_coeff) );
      break;
    }
    case TEMP_COS_SPACE:
    {
      m_temp_space_dep = std::shared_ptr<V_MMS_Space> (new MMS_Space_Cos(t_space_coeff) );
      break;
    }
    case INVALID_TEMPERATURE_SPACE_MMS:
    {
      throw Dark_Arts_Exception(SUPPORT_OBJECT , "Trying to create a temperature MMS with an invalid spatial dependence");
      break;
    }
  }    
  

}

double MMS_Temperature::get_mms_temperature(const double position, const double time)
{
  double time_component = m_time_dep->get_time_component(time);
  double position_component = m_temp_space_dep->get_position_component(position);

  double val = time_component*position_component;
  
  return val;
}

double MMS_Temperature::get_mms_temperature_time_derivative(const double position , const double time)
{
  double time_component = m_time_dep->get_time_derivative(time);
  double position_component = m_temp_space_dep->get_position_component(position);

  double val = time_component*position_component;
  
  return val;
}

