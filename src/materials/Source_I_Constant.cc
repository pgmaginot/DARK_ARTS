/** @file   Source_I_Constant.cc
  *   @author pmaginot
  *   @brief Implement the Source_I_Constant class
  *   \f$ Source_I_Constant = constant\f$
*/
#include "Source_I_Constant.h"

Source_I_Constant::Source_I_Constant(
  const Input_Reader& input_reader, const int mat_num) :
    VSource_I(),
    m_t_start( input_reader.get_rad_source_start(mat_num) ),
    m_t_end( input_reader.get_rad_source_end(mat_num) ),
    m_isotropic_output(input_reader.get_rad_source_output(mat_num) )
{

}

double  Source_I_Constant::get_intensity_source(const double position, 
  const int group, const int dir, const double time)
{
  double val=0.;
  if( ( (time > m_t_start) && (time < m_t_end)) || ( fabs( (time-m_t_end)/(m_t_end - m_t_start)) < 0.000)  )
  {
    val = m_isotropic_output;
  }
  else
  {
    val  =0.;
  }
  
  return val;
}
