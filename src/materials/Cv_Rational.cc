/** @file   Cv_Constant.cc
  *   @author pmaginot
  *   @brief Implement the Cv_Constant class
  *   \f$ C_v = constant\f$
*/
#include "Cv_Rational.h"

Cv_Rational::Cv_Rational(
  const Input_Reader& input_reader, const int mat_num) :
    m_const{ input_reader.get_cv_constant(mat_num) },
    m_power{ input_reader.get_rational_cv_power(mat_num) },
    m_offset{ input_reader.get_rational_cv_offset(mat_num) }
{
}


double  Cv_Rational::get_cv( const double position, const double temperature)
{
  return   m_const/(pow(temperature, m_power) + m_offset);
}
