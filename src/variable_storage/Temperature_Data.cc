/** @file   Temperature_Data.cc
  *   @author pmaginot
  *   @brief Implement the Temperature_Data class
  *   Store temperature unknowns
*/
#include "Temperature_Data.h"

Temperature_Data::Temperature_Data(const Cell_Data& cell_data, const Fem_Quadrature& fem_quad)
  /// initilaize range members
  : m_cells{cell_data.get_total_number_of_cells() } ,     
    m_el_per_cell{fem_quad.get_number_of_interpolation_points() }  
  {
    m_t_length = m_cells*m_el_per_cell;
    
    m_t.resize(m_t_length,0.);
  }
  
bool Temperature_Data::temperature_range_check(const int el, const int cell) const
{
  bool is_bad = false;
  
  if( (el >= m_el_per_cell ) || (el < 0) )
    is_bad = true;
       
  if( (cell < 0) || (cell >= m_cells) )
    is_bad = true;
  
  return is_bad;
}
  
int Temperature_Data::temperature_data_locator(const int el, const int cell) const
{
  int loc_val = -1;
  
  /// layout temperature unknowns from left to right
  loc_val = cell*m_el_per_cell + el;
  
  return loc_val;
}
   
bool Temperature_Data::temperature_bounds_check(const int loc) const
{
  bool is_bad_loc = false;
  if( (loc < 0) || (loc >= m_t_length) )
    is_bad_loc = true;  
  
  return is_bad_loc;
}
