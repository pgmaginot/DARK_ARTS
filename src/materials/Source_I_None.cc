/** @file   Source_I_None.cc
  *   @author pmaginot
  *   @brief Implement the Source_I_None class
  *   \f$ Source_I_None =0\f$
*/
#include "Source_I_None.h"

Source_I_None::Source_I_None(
  const Input_Reader& input_reader, const int mat_num) :
    VSource_I()
{

}

double  Source_I_None::get_intensity_source(const double position, 
  const int group, const int dir, const double time)
{
  return 0.;
}
