/** @file   Intensity_Data.cc
  *   @author pmaginot
  *   @brief Implement the Intensity_Data class
  *   Store the group intensity (I); angle integrated group intensity (phi)
*/
#include "Intensity_Data.h"

Intensity_Data::Intensity_Data(const Cell_Data& cell_data, const Angular_Quadrature& ang_quad,
    const Fem_Quadrature& fem_quad)
  /// initilaize range members
  : m_cells{cell_data.get_total_number_of_cells() } , 
    m_groups{ang_quad.get_number_of_groups() } , 
    m_dir{ang_quad.get_number_of_dir() } , 
    m_leg{ ang_quad.get_number_of_leg_moments() }, 
    m_el_per_cell{fem_quad.get_number_of_interpolation_points() }  
  {
    m_i_length = m_cells*m_groups*m_dir*m_el_per_cell;
    m_phi_length = m_cells*m_groups*m_leg*m_el_per_cell;
    
    m_i.resize(m_i_length,0.);
    m_phi.resize(m_phi_length,0.);
  }
/**
  we want to lay out the angular flux so that we minimize memory movement access
  
  for cell=0:1:<end_of_mesh
    for g=0:1:<n_groups
      for d=n_dir_div_2:1:<n_dir // that we are sweeping in this direction
  
  mu > 0  
  cell*n_groups*n_dir*n_dfem_interp + g*n_dir*n_dfem_interp + d*n_dfem_interp + el
  
  mu < 0
    
*/
double Intensity_Data::get_intensity(const int el, const int cell,
    const int group, const int dir) const
{
  bool bad_input = intensity_range_check(el,cell,group,dir);
  if(bad_input)
  {
    std::cerr << "Error.  Accessing out of logical range intensity\n";
    exit(EXIT_FAILURE);
  }
  
  int val_loc = intensity_data_locator(el,cell,group,dir);
  bool bad_location = intensity_bounds_check(val_loc);
  
  if(bad_location)
  {
    std::cerr << "Error.  Intensity location out of possible range\n";
    exit(EXIT_FAILURE);
  }
  
  return m_i[val_loc];
}

void Intensity_Data::set_intensity(const int el, const int cell,
    const int group, const int dir, const double val) 
{
  
  
  int loc = intensity_data_locator(el,cell,group,dir);
  bool bad_location = intensity_bounds_check(loc);
  
  if(bad_location)
  {
    std::cerr << "Error.  Trying to write to an intensity location out of range\n";
    exit(EXIT_FAILURE);
  }
  
  m_i[loc] = val;
  return;
}

double Intensity_Data::get_angle_integrated_intensity(const int el, const int cell,
    const int group, const int l_mom) const
{
  bool bad_input = angle_integrated_range_check(el,cell,group,l_mom);
  if(bad_input)
  {
    std::cerr << "Error.  Attempting to get out of logical range angle integrated intensity\n";
    exit(EXIT_FAILURE);
  }
  
  int val_loc = angle_integrated_data_locator(el,cell,group,l_mom);
  
  return m_phi[val_loc];
}

void Intensity_Data::set_angle_integrated_intensity(const int el, const int cell,
    const int group, const int l_mom, const double val) 
{
  bool bad_input = angle_integrated_range_check(el,cell,group,l_mom);
  if(bad_input)
  {
    std::cerr << "Error.  Attempting to set out of logical range angle integrated intensity\n";
    exit(EXIT_FAILURE);
  }
  
  int val_loc = angle_integrated_data_locator(el,cell,group,l_mom);
  m_phi[val_loc] = val;
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






