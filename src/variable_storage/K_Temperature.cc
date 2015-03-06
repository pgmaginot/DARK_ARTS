/** @file   K_Temperature.cc
  *   @author pmaginot
  *   @brief Implement the K_Temperature class
  *   Store k quantities for DIRK time integration, for temperature data
*/
#include "K_Temperature.h"

K_Temperature::K_Temperature(const int n_cells, const int n_stages, 
    const Fem_Quadrature& fem_quadrature)
  /// initilaize range members
  : m_cells(n_cells) ,
    m_el_per_cell(fem_quadrature.get_number_of_interpolation_points() ) ,
    m_n_stages(n_stages ),
    m_k_length( m_cells*m_el_per_cell*m_n_stages )   ,
    m_vec_sum( Eigen::VectorXd::Zero(m_el_per_cell) ),
    m_vec_retrieve( Eigen::VectorXd::Zero(m_el_per_cell) ),
    m_kt(m_k_length,0.),
    m_rk_b(m_n_stages,0.)
  {   
    
  }
  
void K_Temperature::clear_kt(void)
{
  for(int i = 0 ; i < m_k_length ; i++)
    m_kt[i] = 0.;
    
  return;
}
  
/// Public accessor functions
void K_Temperature::get_kt(const int cell, const int stage, Eigen::VectorXd& kt) const
{  
  /// find first element location
  int base_loc = kt_data_locator(cell, stage);
  for(int i=0; i<m_el_per_cell ; i++)
    kt(i) = m_kt[base_loc + i];
    
  return;
}

/// Public functions to save values
void K_Temperature::set_kt(const int cell, const int stage, Eigen::VectorXd& kt)
{
  /// find first element location
  int base_loc = kt_data_locator(cell, stage);
  for(int i=0; i<m_el_per_cell; i++)
    m_kt[base_loc + i]= kt(i);   
   
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
    throw Dark_Arts_Exception( VARIABLE_STORAGE , "Attempting to access illogical K_T location");
  
  /// layout temperature unknowns from left to right
  loc_val = cell*m_el_per_cell*m_n_stages + stage*m_el_per_cell;
  
  if( kt_bounds_check(loc_val) )
   throw Dark_Arts_Exception( VARIABLE_STORAGE , "Calculated K_T memory index exceeds the size of K_T data");
  
  return loc_val;
}
   
bool K_Temperature::kt_bounds_check(const int loc) const
{
  bool is_bad_loc = false;
  if( (loc < 0) || (loc >= m_k_length) )
    is_bad_loc = true;  
  
  return is_bad_loc;
}

void K_Temperature::advance_temperature(Temperature_Data& t_old, const double dt, const Time_Data& time_data)
{
  /// loop over all stages in a given cell then advance (based on layout of k_t  
  for(int cell = 0; cell < m_cells ; cell++)
  {    
    m_vec_sum = Eigen::VectorXd::Zero(m_el_per_cell);
    for(int stage = 0; stage < m_n_stages; stage++)
    {         
      m_vec_retrieve = Eigen::VectorXd::Zero(m_el_per_cell);
      get_kt(cell,stage,m_vec_retrieve);
      // std::cout << "K_T of cell " << cell << " stage: " << stage << std::endl << m_vec_retrieve << std::endl;
      m_vec_sum += dt*time_data.get_b(stage)*m_vec_retrieve;
    }
    t_old.get_cell_temperature(cell,m_vec_retrieve);
    m_vec_sum += m_vec_retrieve;
    t_old.set_cell_temperature(cell,m_vec_sum);
  }
  
  return;
}

