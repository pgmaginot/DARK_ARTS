/** @file   Intensity_Data.cc
  *   @author pmaginot
  *   @brief Implement the Intensity_Data class
  *   Store the group intensity (I); angle integrated group intensity (phi)
*/
#include "Intensity_Data.h"

Intensity_Data::Intensity_Data(const Cell_Data& cell_data, const Angular_Quadrature& ang_quad,
  const Fem_Quadrature& fem_quad)
  /// initilaize range members
  : 
  m_cells{cell_data.get_total_number_of_cells() } , 
  m_groups{ang_quad.get_number_of_groups() } , 
  m_dir{ang_quad.get_number_of_dir() } , 
  m_leg{ ang_quad.get_number_of_leg_moments() }, 
  m_el_per_cell{fem_quad.get_number_of_interpolation_points() }  ,
  m_dir_div_2{m_dir/2},    
  m_offset{ m_dir_div_2*m_el_per_cell*m_dir*m_groups*m_cells },
  m_el_times_dir_div_2{ m_dir_div_2*m_el_per_cell }, 
  m_el_times_dir_div_2_times_grp{m_el_times_dir_div_2 * m_groups},
  m_el_times_l_mom{m_leg*m_el_per_cell},
  m_el_times_l_mom_times_group{m_el_times_l_mom*m_groups},    
  m_i_length{m_cells*m_groups*m_dir*m_el_per_cell},
  m_phi_length{m_cells*m_groups*m_leg*m_el_per_cell}
  {    
    m_i.resize(m_i_length,0.);
    m_phi.resize(m_phi_length,0.);
  }


void Intensity_Data::get_cell_intensity(const int cell, const int group, 
  const int dir, Eigen::VectorXd& loc_i_vec) const
{
  bool bad_input = intensity_range_check(0,cell,group,dir);
  if(bad_input)
  {
    std::cerr << "Error.  Accessing out of logical range intensity\n";
    exit(EXIT_FAILURE);
  }
  
  int val_loc = intensity_data_locator(0,cell,group,dir);
  bool bad_location = intensity_bounds_check(val_loc);
  
  if(bad_location)
  {
    std::cerr << "Error.  Intensity location out of possible range\n";
    exit(EXIT_FAILURE);
  }
  
  for(int i=0; i<m_el_per_cell; i++)
    loc_i_vec(i) = m_i[val_loc+i];
    
  return;
}

void Intensity_Data::get_cell_angle_integrated_intensity(const int cell,
  const int group, const int l_mom, Eigen::VectorXd& loc_phi_vec) const
{
  bool bad_input = angle_integrated_range_check(0,cell,group,l_mom);
  if(bad_input)
  {
    std::cerr << "Error.  Attempting to get out of logical range angle integrated intensity\n";
    exit(EXIT_FAILURE);
  }
  
  int val_loc = angle_integrated_data_locator(0,cell,group,l_mom);
  
  for(int i=0; i< m_el_per_cell; i++)
    loc_phi_vec(i) = m_phi[val_loc+i];

  return;
}


void Intensity_Data::set_cell_intensity(const int cell,
  const int group, const int dir, const Eigen::VectorXd& val) 
{  
  int loc = intensity_data_locator(0,cell,group,dir);
  bool bad_location = intensity_bounds_check(loc);
  
  if(bad_location)
  {
    std::cerr << "Error.  Trying to write to an intensity location out of range\n";
    exit(EXIT_FAILURE);
  }
  
  for(int i=0; i<m_el_per_cell ; i++)
    m_i[loc+i] = val(i);
    
  return;
}


void Intensity_Data::set_cell_angle_integrated_intensity(const int cell,
    const int group, const int l_mom, const Eigen::VectorXd& val) 
{
  bool bad_input = angle_integrated_range_check(0,cell,group,l_mom);
  if(bad_input)
  {
    std::cerr << "Error.  Attempting to set out of logical range angle integrated intensity\n";
    exit(EXIT_FAILURE);
  }
  
  int val_loc = angle_integrated_data_locator(0,cell,group,l_mom);
  
  for(int i=0; i< m_el_per_cell ; i++)
    m_phi[val_loc+i] = val(i);
    
  return;
}

  /* ***************************************************
  *
  *   Protected Functions
  *
  *************************************************** */

bool Intensity_Data::intensity_range_check(const int el, const int cell, 
  const int grp, const int dir) const
{
  bool is_bad = false;
  
  if( (el >= m_el_per_cell ) || (el < 0) )
    is_bad = true;
    
  if( (grp < 0) || (grp >= m_groups) )
    is_bad = true;
    
  if( (cell < 0) || (cell >= m_cells) )
    is_bad = true;
    
  if( (dir < 0) || (dir >= m_dir) )
    is_bad = true;
  
  return is_bad;
}

bool Intensity_Data::angle_integrated_range_check(const int el, const int cell, 
  const int grp, const int l_mom) const
{
  bool is_bad = false;
  
  if( (el >= m_el_per_cell ) || (el < 0) )
    is_bad = true;
    
  if( (grp < 0) || (grp >= m_groups) )
    is_bad = true;
    
  if( (cell < 0) || (cell >= m_cells) )
    is_bad = true;
    
  if( (l_mom < 0) || (l_mom >= m_leg) )
    is_bad = true;
  
  return is_bad;
}

/// This function controls the layout of intensity in memory!!
int Intensity_Data::intensity_data_locator(const int el, const int cell, const int group, const int dir) const
{
  int loc = -1;
  
  /**
    Arrange intensity data as follows:
    From closest together to farthest apart:
    element
    direction (by positive/negative)
    group
    cell
    
    Sweeps will do all the directions in a given cell for a group, all groups, then move to the next cell
    
    Upwinding values will be saved in a vector = N_dir/2 to avoid scanning through the intensity data
  
    This will hopefully minizmize data movement
    
  */
  if(dir < m_dir_div_2)
  {
    /// mu < 0
    loc = el + dir*m_el_per_cell + group*m_el_times_dir_div_2 + (m_cells-cell-1)*m_el_times_dir_div_2_times_grp;
  }
  else
  {
    /// mu > 0
    loc = m_offset + el + (dir-m_dir_div_2)*m_el_per_cell + group*m_el_times_dir_div_2 + cell*m_el_times_dir_div_2_times_grp;
  }
  
  bool bad_location = intensity_bounds_check(loc);
  
  if(bad_location)
  {
    std::cerr << "Error.  Intensity location out of possible range\n";
    exit(EXIT_FAILURE);
  }
  
  return loc;
}

/// This function controls the layout of angle_integrated intensities in memory!!
int Intensity_Data::angle_integrated_data_locator(const int el, const int cell, const int group, const int l_mom) const
{
  int loc = -1;
  
  /**
    Arrange data as:
    element
    moment
    group
    cell
  */
  loc = el + l_mom*m_el_per_cell + group*m_el_times_l_mom + cell*m_el_times_l_mom_times_group;
  
  bool bad_location = angle_integrated_bounds_check(loc);
  if(bad_location)
  {
    std::cerr << "Error.  Angle integrated intensity location out of possible range\n";
    exit(EXIT_FAILURE);
  }
  
  return loc;
}

bool Intensity_Data::intensity_bounds_check(const int loc) const
{
  bool is_bad_loc = false;
  if( (loc < 0) || (loc >= m_i_length) )
    is_bad_loc = true;  
  
  return is_bad_loc;
}

bool Intensity_Data::angle_integrated_bounds_check(const int loc) const
{
  bool is_bad_loc = false;
  if( (loc < 0) || (loc >= m_phi_length) )
    is_bad_loc = true;  
  
  return is_bad_loc;
}






