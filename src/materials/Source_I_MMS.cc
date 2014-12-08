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
  m_c(planck.get_c() ),
  m_ang_quad(angular_quadrature),
  m_intensity(input_reader, angular_quadrature),
  m_temperature(input_reader) ,
  m_abs_op(abs_op[mat_num]),
  m_scat_op(scat_op[mat_num])
{
  
}

double  Source_I_MMS::get_intensity_source(const double position, 
  const int group, const int dir, const double time)
{
  /**
    \mu \frac{\partial I}{\partial x} + \sigma_t I - \sigma_s \phi - \sigma_a planck = S_I
  */
  // std::cout << "Group: " << group << " Dir: " << dir << std::endl; 
  std::cout << "Position: " << position << std::endl;
  
  double val = 0.;
  double temperature = m_temperature.get_mms_temperature(position,time);
  
  val = m_intensity.get_mms_intensity_time_derivative(position,time,dir)/m_c;
  
  val+= m_ang_quad.get_mu(dir)*m_intensity.get_mms_intensity_space_derivative(position,time,dir);
  
  double sig_s = m_scat_op->get_scattering_opacity( 0, 0, temperature, position);
  double sig_a = m_abs_op->get_absorption_opacity(0,temperature,position);
  
  std::cout << "sig_s: " << sig_s <<std::endl;
  std::cout << "time: " << time <<std::endl;
  
  val += (sig_s + sig_a)*m_intensity.get_mms_intensity(position,time,dir) ;
  val -= sig_s/m_ang_quad.get_sum_w()*m_intensity.get_mms_phi(position,time);
  val -= sig_a*m_planck.integrate_B_grey(temperature);  
  
  return val;
}
