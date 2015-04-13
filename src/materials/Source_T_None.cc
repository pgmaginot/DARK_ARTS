/** @file   Source_T_None.cc
  *   @author pmaginot
  *   @brief Implement the Source_T_Constant class
  *   \f$ Source_T_None = constant\f$
*/
#include "Source_T_None.h"

Source_T_None::Source_T_None(
  const Input_Reader& input_reader, const int mat_num) :
    VSource_T()
{
  
}

double  Source_T_None::get_temperature_source(const double position, const double time)
{
  return 0.;
}
