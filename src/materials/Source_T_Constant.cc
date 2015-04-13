/** @file   Source_T_Constant.cc
  *   @author pmaginot
  *   @brief Implement the Source_T_Constant class
  *   \f$ Source_T_Constant = constant\f$
*/
#include "Source_T_Constant.h"

Source_T_Constant::Source_T_Constant(
  const Input_Reader& input_reader, const int mat_num) :
    VSource_T(),    
    m_t_start( input_reader.get_temp_source_start(mat_num) ),
    m_t_end( input_reader.get_temp_source_end(mat_num) ),
    m_output(input_reader.get_temp_source_output(mat_num) )
{

}


double  Source_T_Constant::get_temperature_source(const double position, const double time)
{
  double val=0.;
  
  if( (time >= m_t_start) && (time <= m_t_end))
    val = m_output;
  else
    val  = 0.;
    
  return val;
}
