/** @file   K_Temperature.cc
  *   @author pmaginot
  *   @brief Implement the K_Temperature class
  *   Store k quantities for DIRK time integration, for temperature data
*/
#include "K_Temperature.h"

K_Temperature::K_Temperature(const Cell_Data& cell_data, const int n_stages, 
    const Fem_Quadrature& fem_quadrature)
  /// initilaize range members
  : m_cells(cell_data.get_total_number_of_cells() ) ,
    m_el_per_cell(fem_quadrature.get_number_of_interpolation_points() ) ,
    m_n_stages(n_stages ),
    m_k_length( m_cells*m_el_per_cell*m_n_stages )
    
  {   
    m_kt.resize(m_k_length,0.);
  }
  
/// Public accessor functions
void K_Temperature::get_kt(const int cell, const int stage, Eigen::VectorXd& kt) const
{  
  /// find first element location
  int base_loc = kt_data_locator(cell, stage);
  for(int i=0;i<m_el_per_cell;i++)
    kt(i) = m_kt[base_loc + i];
    
  return;
}

/// Public functions to save values
void K_Temperature::set_kt(const int cell, const int stage, Eigen::VectorXd& kt)
{
  /// find first element location
  int base_loc = kt_data_locator(cell, stage);
  for(int i=0;i<m_el_per_cell;i++)
    m_kt[base_loc + i]=kt(i);
   
   
  return ;
}
bool K_Temperature::kt_range_check(const int cell, const int stage) const
{
  bool is_bad = false;
  
  if( (stage >= m_n_stages ) || (stage < 0) )
    is_bad = true;
       
  if( (cell < 0) || (cell >= m_cells) )
    is_bad = true;
  
  return is_bad;
}
  
int K_Temperature::kt_data_locator(const int cell, const int stage) const
{
  int loc_val = -1;
  
  if( kt_range_check(cell,stage) )
  {
    std::cerr << "Attempting to access illogical temperature location\n";
    exit(EXIT_FAILURE);
  }
  
  /// layout temperature unknowns from left to right
  loc_val = cell*m_el_per_cell*m_n_stages + stage*m_el_per_cell;
  
  if( kt_bounds_check(loc_val) )
  {
    std::cerr << "Location out of bounds in temperature data\n";
    exit(EXIT_FAILURE);
  }
  
  return loc_val;
}
   
bool K_Temperature::kt_bounds_check(const int loc) const
{
  bool is_bad_loc = false;
  if( (loc < 0) || (loc >= m_k_length) )
    is_bad_loc = true;  
  
  return is_bad_loc;
}
