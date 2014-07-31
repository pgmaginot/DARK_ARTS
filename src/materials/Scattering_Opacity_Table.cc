/** @file   Scattering_Opacity_Table.cc
  *   @author pmaginot
  *   @brief Implement the Scattering_Opacity_Table class
  *   opacity is a table look-up value
*/
#include "Scattering_Opacity_Table.h"

Scattering_Opacity_Table::Scattering_Opacity_Table(
  const Input_Reader& input_reader, const int mat_num) 
{
  input_reader.get_scat_file_str(mat_num,m_filename);
}

Scattering_Opacity_Table::~Scattering_Opacity_Table(){}

double  Scattering_Opacity_Table::get_scattering_opacity(const int l_mom,
  const int group, const double temperature, const double position)
{
  std::cerr << "Table look-up has not been coded yet.  Error.\n";
  exit(EXIT_FAILURE);
  return 0.;
}
