#include "Transport_BC_MMS.h"

Transport_BC_MMS::Transport_BC_MMS(const Angular_Quadrature& angular_quadrature,
  const Input_Reader& input_reader, const double x_b)
:
  V_Transport_BC(),  
  m_x_boundary(x_b)
{
  RADIATION_SPACE_MMS r_space_dependence = input_reader.get_mms_radiation_space_dependence() ; 
  TIME_MMS_TYPE time_dependence = input_reader.get_mms_time_dependence();
  RADIATION_ANGLE_MMS angle_dependence = input_reader.get_mms_radiation_angle_dependence();
  
  std::vector<double> time_coeff;
  input_reader.get_mms_time_coeff(time_coeff);
  
  std::vector<double> r_space_coeff;
  input_reader.get_mms_radiation_space_coeff(r_space_coeff);  
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
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Trying to create an invalid MMS radiation time dependence");
      break;
    }
  }
  switch(r_space_dependence)
  {
    case RAD_POLY_SPACE:
    {
      m_rad_space_dep = std::shared_ptr<V_MMS_Space> (new MMS_Space_Poly(r_space_coeff) );
      break;
    }
    case RAD_COS_SPACE:
    {
      m_rad_space_dep = std::shared_ptr<V_MMS_Space> (new MMS_Space_Cos(r_space_coeff) );
      break;
    }
    case INVALID_RADIATION_SPACE_MMS:
    {
      throw Dark_Arts_Exception(SUPPORT_OBJECT , "Trying to create a radiation MMS with an invalid spatial dependence");
      break;
    }
  }    
  switch(angle_dependence)
  {
    case MMS_ISOTROPIC:
    {
      m_angle_dep = std::shared_ptr<V_MMS_Angle> (new MMS_Angle_Isotropic(angular_quadrature.get_sum_w() ) );
      break;
    }
    case MMS_ANGLE_POLY:
    {
      std::vector<double> angle_coeff;
      input_reader.get_mms_angle_coeff(angle_coeff);
      m_angle_dep = std::shared_ptr<V_MMS_Angle> (new MMS_Angle_Poly(angle_coeff,angular_quadrature) );
      break;
    }
    case INVALID_RADIATION_ANGLE_MMS:
    {
      throw Dark_Arts_Exception(SUPPORT_OBJECT, "Invalid radiation angle mms dependence");
      break;
    }
  }
}

double Transport_BC_MMS::get_boundary_condition(const double mu, const int grp , const double time) 
{
  return m_angle_dep->get_angle_component(mu)*m_rad_space_dep->get_position_component(m_x_boundary)*m_time_dep->get_time_component(time);
}

