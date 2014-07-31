/** @file   Scattering_Opacity_Constant.cc
  *   @author pmaginot
  *   @brief Implement the Scattering_Opacity_Constant class
  *   \f$ \sigma_s = constant\f$
*/
#include "Scattering_Opacity_Constant.h"

Scattering_Opacity_Constant::Scattering_Opacity_Constant(
  const Input_Reader& input_reader, const int mat_num) :
    m_const{ input_reader.get_scat_double_constant_1(mat_num) }
{
  if(m_const < 0. )
  {
    std::cerr << "Invalid scattering opacity constant in material " << mat_num << std::endl;
    exit(EXIT_FAILURE);
  }
}

Scattering_Opacity_Constant::~Scattering_Opacity_Constant(){}

double  Scattering_Opacity_Constant::get_scattering_opacity(const int l_mom, const int group, 
  const double temperature, const double position)
{
  return m_const;
}
