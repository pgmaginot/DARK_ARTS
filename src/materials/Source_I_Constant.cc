/** @file   Source_I_Constant.cc
  *   @author pmaginot
  *   @brief Implement the Source_I_Constant class
  *   \f$ Source_I_Constant = constant\f$
*/
#include "Source_I_Constant.h"

Source_I_Constant::Source_I_Constant(
  const Input_Reader& input_reader, const int mat_num) :
    VSource_I(),
    m_const{ 0. }
{
  if(m_const < 0. )
  {
    std::stringstream err;
    err    << "Invalid radiation source in material " << mat_num ;
    throw Dark_Arts_Exception( SUPPORT_OBJECT , err.str() );
  }
}

double  Source_I_Constant::get_intensity_source(const double position, 
  const int group, const int dir, const double time)
{
  return m_const;
}
