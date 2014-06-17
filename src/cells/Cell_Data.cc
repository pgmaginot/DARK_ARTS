/** @file   Cell_Data.cc
  *   @author pmaginot
  *   @brief Implement the Cell_Data class, holds all cell information (xL,xR,material_number)
*/
#include "Cell_Data.h"

Cell_Data::Cell_Data(Input_Reader&  input_reader)
{
  std::vector<int> cells_per_region;
  input_reader.get_cells_per_region_vector(cells_per_region);
  
  std::vector<double> region_left_bound;
  std::vector<double> region_right_bound;
  input_reader.get_region_boundaries(region_left_bound,region_right_bound);
  
  std::vector<GRID_SPACING> region_spacing;
  input_reader.get_region_spacing(region_spacing);
  
  std::vector<int> region_material_num;
  input_reader.get_region_materials(region_material_num);
  
  int n_region = input_reader.get_n_regions();
  for(int i=0;i<n_region;i++)
  {
    m_total_cells += cells_per_region[i];
  }
  
  /// initialize vectors of data
  m_x_l.resize(m_total_cells,0.);
  m_x_r.resize(m_total_cells,0.);
  m_material_num.resize(m_total_cells,0);
}

