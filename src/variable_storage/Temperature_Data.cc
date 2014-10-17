/** @file   Temperature_Data.cc
  *   @author pmaginot
  *   @brief Implement the Temperature_Data class
  *   Store temperature unknowns
*/
#include "Temperature_Data.h"

Temperature_Data::Temperature_Data(const int n_cells, const Fem_Quadrature& fem_quad)
  /// initilaize range members
  : m_cells{ n_cells } ,     
    m_el_per_cell{fem_quad.get_number_of_interpolation_points() },
    m_t_length{ m_cells*m_el_per_cell} ,
    m_t(m_t_length,0.)    
  {
    fem_quad.get_dfem_interpolation_point_weights(m_dfem_w);
  }
  
/// Public accessor functions
double Temperature_Data::get_temperature(const int el, const int cell) const
{  
  return m_t[temperature_data_locator(el,cell)];
}

void Temperature_Data::set_temperature(const int el, const int cell, const double val)
{
  int loc = temperature_data_locator(el,cell);
  m_t[loc] = val;
  return ;
}

void Temperature_Data::get_cell_temperature(const int cell, Eigen::VectorXd& vec) const
{  
  int base = temperature_data_locator(0,cell);
  for(int i=0; i< m_el_per_cell; i++ )
    vec(i) = m_t[base + i];
  return; 
}

void Temperature_Data::set_cell_temperature(const int cell, const Eigen::VectorXd& vec)
{
  int loc = temperature_data_locator(0,cell);
  for(int i=0; i< m_el_per_cell ; i++)
    m_t[loc+i] = vec(i);
    
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

Temperature_Data& Temperature_Data::operator= (const Temperature_Data& t_data)
{
  if( (m_cells != t_data.m_cells) ||
      (m_el_per_cell != t_data.m_el_per_cell) ||
      (m_t_length != t_data.m_t_length) )
  {
    std::cerr << "Trying to copy non-identical Temperature_Data objects \n";
    exit(EXIT_FAILURE);
  }
  for(int i=0; i< m_t_length; i++)
    m_t[i] = t_data.m_t[i];
  
  return *this;
}

double Temperature_Data::calculate_average(void)
{
  double val = 0.;
  int cnt = 0;
  for(int c=0; c<m_cells; c++)
  {
    for(int el=0;el<m_el_per_cell;el++)
    {
      val += m_t[cnt]*m_dfem_w[el];
      cnt++;
    }
  }
  val /= (2.* double(m_cells));
  return val;
}

