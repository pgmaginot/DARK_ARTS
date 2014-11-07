/** @file   Cv_Constant.cc
  *   @author pmaginot
  *   @brief Implement the Cv_Constant class
  *   \f$ C_v = constant\f$
*/
#include "Cv_Constant.h"

Cv_Constant::Cv_Constant(
  const Input_Reader& input_reader, const int mat_num) :
    m_const{ input_reader.get_cv_constant(mat_num) }
{
  if(m_const < 0. )
  {
    std::stringstream err;
    err << "Invalid absorption opacity constant in material " << mat_num ;
    throw Dark_Arts_Exception( SUPPORT_OBJECT , err.str() );
  }
}

Cv_Constant::~Cv_Constant(){}

double  Cv_Constant::get_cv( const double position, const double temperature)
{
  return m_const;
}
