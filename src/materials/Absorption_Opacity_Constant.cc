/** @file   Absorption_Opacity_Constant.cc
  *   @author pmaginot
  *   @brief Implement the Absorption_Opacity_Constant class
  *   \f$ \sigma_a = constant\f$
*/
#include "Absorption_Opacity_Constant.h"

Absorption_Opacity_Constant::Absorption_Opacity_Constant(
  const Input_Reader& input_reader, const int mat_num) :
    m_const{ input_reader.get_abs_double_constant_1(mat_num) }
{
  if(m_const < 0. )
  {
    std::stringstream err;
    err << "Invalid absorption opacity constant in material " << mat_num;
    throw Dark_Arts_Exception( SUPPORT_OBJECT, err.str() ) ;
  }
}

Absorption_Opacity_Constant::~Absorption_Opacity_Constant(){}

double  Absorption_Opacity_Constant::get_absorption_opacity(const int group, 
  const double temperature, const double position)
{
  return m_const;
}
