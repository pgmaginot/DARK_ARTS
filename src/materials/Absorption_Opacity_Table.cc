/** @file   Absorption_Opacity_Table.cc
  *   @author pmaginot
  *   @brief Implement the Absorption_Opacity_Table class
  *   opacity is a table look-up value
*/
#include "Absorption_Opacity_Table.h"

Absorption_Opacity_Table::Absorption_Opacity_Table(
  const Input_Reader& input_reader, const int mat_num) 
{
  input_reader.get_abs_file_str(mat_num,m_filename);
}

Absorption_Opacity_Table::~Absorption_Opacity_Table(){}

double  Absorption_Opacity_Table::get_absorption_opacity(const int group, 
  const double temperature, const double position)
{
  throw Dark_Arts_Exception( SUPPORT_OBJECT , "Tabular absorption opacity objects not yet coded");
  
  return 0.;
}
