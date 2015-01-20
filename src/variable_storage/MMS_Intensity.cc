/** @file   MMS_Intensity.cc
  *   @author pmaginot
  *   @brief A convenience class for accesing a manufactured solution angular intensity
*/
#include "MMS_Intensity.h"

MMS_Intensity::MMS_Intensity(const Input_Reader& input_reader, const Angular_Quadrature& angular_quadrature)
:
m_ang_quad(angular_quadrature)
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
      m_angle_dep = std::shared_ptr<V_MMS_Angle> (new MMS_Angle_Poly(angle_coeff, angular_quadrature) );
      break;
    }
    case INVALID_RADIATION_ANGLE_MMS:
    {
      throw Dark_Arts_Exception(SUPPORT_OBJECT, "Invalid radiation angle mms dependence");
      break;
    }
  }

  m_angle_integration = 0.;
  for(int d=0; d<angular_quadrature.get_number_of_dir() ; d++)
    m_angle_integration += angular_quadrature.get_w(d) * m_angle_dep->get_angle_component(d);
    
}

double MMS_Intensity::get_mms_intensity(const double position , const double time , const int dir) const
{
  double time_component = m_time_dep->get_time_component(time);
  double angle_comp = m_angle_dep->get_angle_component(dir);
  double i_position_comp = m_rad_space_dep->get_position_component(position);
  
  double val = time_component*i_position_comp*angle_comp;
  
  return val;
}

double MMS_Intensity::get_mms_intensity(const double position , const double time , const double mu) const
{
  double time_component = m_time_dep->get_time_component(time);
  double angle_comp = m_angle_dep->get_angle_component(mu);
  double i_position_comp = m_rad_space_dep->get_position_component(position);
  
  double val = time_component*i_position_comp*angle_comp;
  
  return val;
}

double MMS_Intensity::get_mms_intensity_space_derivative(const double position , const double time , const int dir) const
{
  double time_component = m_time_dep->get_time_component(time);
  double angle_comp = m_angle_dep->get_angle_component(dir);
  double i_position_comp = m_rad_space_dep->get_position_derivative(position);
  
  double val = time_component*i_position_comp*angle_comp;
  
  return val;
}

double MMS_Intensity::get_mms_intensity_time_derivative(const double position , const double time , const int dir) const
{
  double time_component = m_time_dep->get_time_derivative(time);
  double angle_comp = m_angle_dep->get_angle_component(dir);
  double i_position_comp = m_rad_space_dep->get_position_component(position);
  
  double val = time_component*i_position_comp*angle_comp;
  
  return val;
}

double MMS_Intensity::get_mms_phi(const double position , const double time ) const
{
  double time_component = m_time_dep->get_time_component(time);
  double i_position_comp = m_rad_space_dep->get_position_component(position);

  double val = time_component*i_position_comp*m_angle_integration;
  
  return val;
}



