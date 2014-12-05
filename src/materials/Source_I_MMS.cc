/** @file   Source_I_MMS.cc
  *   @author pmaginot
  *   @brief Implement the Source_I_Constant class
  *   \f$ Source_I_Constant = constant\f$
*/
#include "Source_I_MMS.h"

Source_I_MMS::Source_I_MMS( const Input_Reader& input_reader, const Angular_Quadrature& angular_quadrature,
  std::vector<std::shared_ptr<VAbsorption_Opacity>>& abs_op , 
  std::vector<std::shared_ptr<VScattering_Opacity>>& scat_op, 
  const int mat_num , const Planck& planck)
  :
  VSource_I(),
  m_planck(planck),
  m_ang_quad(angular_quadrature),
  m_c(planck.get_c() )
{
  m_scat_op = scat_op[mat_num];
  m_abs_op = abs_op[mat_num];
  
  TEMPERATURE_SPACE_MMS t_space_dependence = input_reader.get_mms_temperature_space_dependence();
  RADIATION_SPACE_MMS r_space_dependence = input_reader.get_mms_radiation_space_dependence() ; 
  TIME_MMS_TYPE time_dependence = input_reader.get_mms_time_dependence();
  RADIATION_ANGLE_MMS angle_dependence = input_reader.get_mms_radiation_angle_dependence();
  
  std::vector<double> time_coeff;
  input_reader.get_mms_time_coeff(time_coeff);
  
  std::vector<double> r_space_coeff;
  input_reader.get_mms_radiation_space_coeff(r_space_coeff);  
  
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
  
  for(int d = 0; d < m_ang_quad.get_number_of_dir() ; d++)
    m_angle_integral += m_ang_quad.get_w(d)*m_angle_dep->get_angle_component(d);
}

double  Source_I_MMS::get_intensity_source(const double position, 
  const int group, const int dir, const double time)
{
  /**
    \mu \frac{\partial I}{\partial x} + \sigma_t I - \sigma_s \phi - \sigma_a planck = S_I
  */
  double val = 0.;
  double time_component = m_time_dep->get_time_component(time);
  double angle_comp = m_angle_dep->get_angle_component(dir);
  double i_position_comp = m_rad_space_dep->get_position_component(position);
  double temperature = time_component*m_temp_space_dep->get_position_component(position);
  
  val = m_time_dep->get_time_derivative(time)*angle_comp*i_position_comp/m_c;
  
  val+= m_ang_quad.get_mu(dir)*m_rad_space_dep->get_position_derivative(position)*time_component*angle_comp;
  
  double sig_s = m_scat_op->get_scattering_opacity( 0, 0, temperature, position);
  double sig_a = m_abs_op->get_absorption_opacity(0,temperature,position); 
  val += (sig_s + sig_a)*angle_comp*time_component*i_position_comp;
  val -= sig_s/2.*m_angle_integral*time_component*i_position_comp;
  val -= sig_a*m_planck.integrate_B_grey(temperature);  
  
  return val;
}
