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
  const Planck& planck)
  :
  VSource_T(),
  m_planck(planck),
  m_sn_w(angular_quadrature.get_sum_w() ) ,
  m_intensity(input_reader, angular_quadrature),
  m_temperature(input_reader) ,
  m_abs_op(abs_op[mat_num]),
  m_cv(cv[mat_num])
{
  
}

double  Source_T_MMS::get_temperature_source(const double position, const double time)
{  
  double val = 0.;
  /** \f[
    S_T = C_v \frac{\partial T}{\partial t} - \sigma_a \left( \phi - m_sn_w B \right)
    
  */
  double temperature = m_temperature.get_mms_temperature(position,time);
  double cv = m_cv->get_cv(position,temperature);
  double sig_a = m_abs_op->get_absorption_opacity(0,temperature,position);
  
  // std::cout << "cv: " << cv << "sig_a: " << sig_a <<std::endl;
  // std::cout << "d f(t)/dt: " << dt_time <<std::endl;
  // std::cout << "angle integration: " << m_angle_integration << std::endl;
 
  
  val = cv*m_temperature.get_mms_temperature_time_derivative(position,time)
    - sig_a*( m_intensity.get_mms_phi(position,time) 
               - m_sn_w*m_planck.integrate_B_grey(temperature) );
               
  return val;
}
