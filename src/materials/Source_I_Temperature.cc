/** @file   Source_I_Temperature.cc
  *   @author pmaginot
  *   @brief Temperature specified radiation source
*/
#include "Source_I_Temperature.h"

Source_I_Temperature::Source_I_Temperature(const Input_Reader& input_reader, const int mat_num, const double a, const double c, const double sn_sum) :
    VSource_I(),
    m_t_start( input_reader.get_rad_source_start(mat_num) ),
    m_t_end( input_reader.get_rad_source_end(mat_num) ),
    m_output(a*c/sn_sum*pow( input_reader.get_rad_source_output(mat_num) , 4) )
{
  
}

double  Source_I_Temperature::get_intensity_source(const double position, 
  const int group, const int dir, const double time)
{
  double val=0.;
  if( ( (time > m_t_start) && (time < m_t_end)) || ( fabs( (time-m_t_end)/(m_t_end - m_t_start)) < 0.000)  )
  {
    val = m_output;
  }
  else
  {
    val  =0.;
  }
  
  return val;
}
