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
  
/// Public accessor functions
double Temperature_Data::get_temperature(const int el, const int cell) const
{  
  return m_t[temperature_data_locator(el,cell)];
}

/// Public functions to save values
void Temperature_Data::set_temperature(const int el, const int cell, const double val)
{
  int loc = temperature_data_locator(el,cell);
  m_t[loc] = val;
  return ;
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
  
  if( temperature_range_check(el,cell) )
  {
    std::cerr << "Attempting to access illogical temperature location\n";
    exit(EXIT_FAILURE);
  }
  
  /// layout temperature unknowns from left to right
  loc_val = cell*m_el_per_cell + el;
  
  if( temperature_bounds_check(loc_val) )
  {
    std::cerr << "Location out of bounds in temperature data\n";
    exit(EXIT_FAILURE);
  }
  
  return loc_val;
}
   
bool Temperature_Data::temperature_bounds_check(const int loc) const
{
  bool is_bad_loc = false;
  if( (loc < 0) || (loc >= m_t_length) )
    is_bad_loc = true;  
  
  return is_bad_loc;
}
