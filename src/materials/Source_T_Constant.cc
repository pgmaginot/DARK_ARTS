/** @file   Source_T_Constant.cc
  *   @author pmaginot
  *   @brief Implement the Source_T_Constant class
  *   \f$ Source_T_Constant = constant\f$
*/
#include "Source_T_Constant.h"

Source_T_Constant::Source_T_Constant(
  const Input_Reader& input_reader, const int mat_num) :
    m_const{ 0. }
{
  if(m_const < 0. )
  {
    std::stringstream err;
    err    << "Invalid temperature source in material " << mat_num ;
    throw Dark_Arts_Exception( SUPPORT_OBJECT , err.str() );
  }
}

Source_T_Constant::~Source_T_Constant(){}

double  Source_T_Constant::get_temperature_source(const double position, const double time)
{
  return m_const;
}
