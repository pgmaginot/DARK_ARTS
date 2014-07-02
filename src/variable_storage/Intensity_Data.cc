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
/*
double Intensity_Data::get_intensity(const int el, const int cell,
    const int group, const int dir) const
{
  return m_i[0];
}

void Intensity_Data::get_cell_intensity(std::vector<double>& i_cell, const int cell,
    const int group, const int dir) const
{
  return m_i[0];
}

void Intensity_Data::set_intensity(const int el, const int cell,
    const int group, const int dir, const double val) 
{
  return m_i[0] = val;
}


double Intensity_Data::get_angle_integrated_intensity(const int el, const int cell,
    const int group, const int dir) const
{
  return m_phi[0];
}

void Intensity_Data::set_angle_integrated_intensity(const int el, const int cell,
    const int group, const double val) 
{
  return m_phi[0] = val;
}


*/


