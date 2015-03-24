/** @file   Scattering_Opacity_Rational.cc
  *   @author pmaginot
  *   @brief Implement the Scattering_Opacity_Rational class
  *   \f$ \sigma_s = \frac{c_1}{c_2 + T^P} \f$
*/
#include "Scattering_Opacity_Rational.h"

Scattering_Opacity_Rational::Scattering_Opacity_Rational(
  const Input_Reader& input_reader, const int mat_num) :
    m_const{ input_reader.get_scat_double_constant_1(mat_num) },
    m_offset{ input_reader.get_scat_double_constant_2(mat_num) },
    m_p{ input_reader.get_scat_int_constant(mat_num) }
{
  if(m_const < 0. )
  {
    std::stringstream err;
    err <<  "Invalid scattering opacity constant in material " << mat_num ;
    throw Dark_Arts_Exception( SUPPORT_OBJECT , err.str() );
  }
  if(m_offset < 0. )
  {
    std::stringstream err;
    err <<  "Invalid scattering opacity denominator offset in material " << mat_num ;
    throw Dark_Arts_Exception( SUPPORT_OBJECT , err.str() );
  }
  if(m_p < 1 )
  {
    std::stringstream err;
    err <<  "Invalid scattering opacity temperature power in material " << mat_num ;
    throw Dark_Arts_Exception( SUPPORT_OBJECT , err.str() );
  }
}

Scattering_Opacity_Rational::~Scattering_Opacity_Rational(){}

double  Scattering_Opacity_Rational::get_scattering_opacity(const int l_mom, const int group, 
  const double temperature, const double position)
{
  return m_const/(m_offset + pow( fabs(temperature),m_p));
}
