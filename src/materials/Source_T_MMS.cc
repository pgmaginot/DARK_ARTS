/** @file   Source_T_Constant.cc
  *   @author pmaginot
  *   @brief Implement the Source_T_Constant class
  *   \f$ Source_T_Constant = constant\f$
*/
#include "Source_T_MMS.h"

Source_T_MMS::Source_T_MMS(const Input_Reader& input_reader, 
  const Angular_Quadrature& angular_quadrature, 
  std::vector<std::shared_ptr<VAbsorption_Opacity> >& abs_op, 
  std::vector<std::shared_ptr<VCv> >& cv, 
  const int mat_num,
  const Planck& planck):
  m_planck(planck),
  m_sn_w(angular_quadrature.get_sum_w() )
{
  m_abs_op = abs_op[mat_num];
  m_cv = cv[mat_num];
  
  m_angle_integration = 0.;
  int n_dir = angular_quadrature.get_number_of_dir();
  std::shared_ptr<V_MMS_Angle> angle_dependence;
  RADIATION_ANGLE_MMS angle_type = input_reader.get_mms_radiation_angle_dependence();
  switch(angle_type)
  {
    case MMS_ISOTROPIC:
    {
      angle_dependence = std::shared_ptr<V_MMS_Angle> (new MMS_Angle_Isotropic(angular_quadrature.get_sum_w() ) );
      break;
    }
    case MMS_ANGLE_POLY:
    {
      std::vector<double> angle_coeff;
      input_reader.get_mms_angle_coeff(angle_coeff);
      angle_dependence = std::shared_ptr<V_MMS_Angle> (new MMS_Angle_Poly(angle_coeff, angular_quadrature) );
      break;
    }
    case INVALID_RADIATION_ANGLE_MMS:
    {
      throw Dark_Arts_Exception(SUPPORT_OBJECT, "Invalid radiation angle mms dependence");
      break;
    }
  }
  for(int d=0; d<n_dir; d++)
    m_angle_integration += angular_quadrature.get_sum_w() * angle_dependence->get_angle_component(d);
  
  TIME_MMS_TYPE time_dependence = input_reader.get_mms_time_dependence();
  TEMPERATURE_SPACE_MMS temp_space_dependence = input_reader.get_mms_temperature_space_dependence();
  RADIATION_SPACE_MMS rad_space_dependence = input_reader.get_mms_radiation_space_dependence();
  
  std::vector<double> time_coeff;  
  input_reader.get_mms_time_coeff(time_coeff);
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
  
  std::vector<double> t_space_coeff; 
  input_reader.get_mms_temperature_space_coeff(t_space_coeff);
  switch(temp_space_dependence)
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
  
  std::vector<double> r_space_coeff; 
  input_reader.get_mms_radiation_space_coeff(r_space_coeff);
  switch(rad_space_dependence)
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
}

Source_T_MMS::~Source_T_MMS(){}

double  Source_T_MMS::get_temperature_source(const double position, const double time)
{  
  double time_component = m_time_dep->get_time_component(time) ;
  double dt_time = m_time_dep->get_time_derivative(time);
  double temperature_space = m_temp_space_dep->get_position_component(position);
  double temperature = temperature_space * time_component;
  double phi = m_rad_space_dep->get_position_component(position)*m_angle_integration* time_component;
  double val = 0.;
  /** \f[
    S_T = C_v \frac{\partial T}{\partial t} - \sigma_a \left( \phi - m_sn_w B \right)
    
  */
  double cv = m_cv->get_cv(position,temperature);
  double sig_a = m_abs_op->get_absorption_opacity(0,temperature,position);
  val = cv*dt_time*temperature_space - sig_a*( phi - m_planck.integrate_B_grey(temperature) );
  return val;
}
