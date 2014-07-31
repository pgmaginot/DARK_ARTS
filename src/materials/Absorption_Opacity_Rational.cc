/** @file   Absorption_Opacity_Rational.cc
  *   @author pmaginot
  *   @brief Implement the concrete Absorption_Opacity_Rational class
  *   \f$ \sigma_a = \frac{c_1}{c_2 + T^P} \f$
*/
#include "Absorption_Opacity_Rational.h"

Absorption_Opacity_Rational::Absorption_Opacity_Rational( const Input_Reader& input_reader, 
  const int mat_num) :
    m_const{ input_reader.get_abs_double_constant_1(mat_num) },
    m_offset{ input_reader.get_abs_double_constant_2(mat_num) },
    m_p{ input_reader.get_abs_int_constant(mat_num) }
{
  if(m_const < 0. )
  {
    std::cerr << "Invalid absorption opacity constant in material " << mat_num << std::endl;
    exit(EXIT_FAILURE);
  }
  if(m_offset < 0. )
  {
    std::cerr << "Invalid absorption opacity denominator offset in material " << mat_num << std::endl;
    exit(EXIT_FAILURE);
  }
  if(m_p < 1 )
  {
    std::cerr << "Invalid absorption opacity temperature power in material " << mat_num << std::endl;
    exit(EXIT_FAILURE);
  }
}

Absorption_Opacity_Rational::~Absorption_Opacity_Rational()
{

}

double Absorption_Opacity_Rational::get_absorption_opacity(const int group, const double temperature, const double position)
{
  return m_const/(m_offset + pow(temperature,m_p));
}

