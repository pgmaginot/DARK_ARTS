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
  m_el_per_cell{fem_quad.get_number_of_interpolation_points() }  ,
  m_dir_div_2{m_dir/2},    
  m_offset{ m_dir_div_2*m_el_per_cell*m_dir*m_groups*m_cells },
  m_el_times_dir_div_2{ m_dir_div_2*m_el_per_cell }, 
  m_el_times_dir_div_2_times_grp{m_el_times_dir_div_2 * m_groups}, 
  m_i_length{m_cells*m_groups*m_dir*m_el_per_cell},
  m_i(m_i_length,0.)
  {    
  }
  
Intensity_Data::Intensity_Data(const Cell_Data& cell_data, 
  const Angular_Quadrature& ang_quad,
  const Fem_Quadrature& fem_quad, 
  Materials& materials,
  const Input_Reader& input_reader)
  /// initilaize range members
  : 
  m_cells{cell_data.get_total_number_of_cells() } , 
  m_groups{ang_quad.get_number_of_groups() } , 
  m_dir{ang_quad.get_number_of_dir() } , 
  m_el_per_cell{fem_quad.get_number_of_interpolation_points() }  ,
  m_dir_div_2{m_dir/2},    
  m_offset{ m_dir_div_2*m_el_per_cell*m_dir*m_groups*m_cells },
  m_el_times_dir_div_2{ m_dir_div_2*m_el_per_cell }, 
  m_el_times_dir_div_2_times_grp{m_el_times_dir_div_2 * m_groups}, 
  m_i_length{m_cells*m_groups*m_dir*m_el_per_cell},
  m_i(m_i_length,0.)
  {    
    /// load the initial conditions
    if(input_reader.get_radiation_ic_type() == PLANCKIAN_IC)
    {
      /// loop over regions.  Each region could have a different initial condition
      std::vector<int> cell_per_reg;
      std::vector<double> temp_in_reg;
      input_reader.get_cells_per_region_vector(cell_per_reg);
      // input_reader.get_region_temperature(temp_in_reg);
      
      int cell_cnt = 0;
      double iso_emission = 0.;
      for(int reg = 0 ; reg < input_reader.get_n_regions() ; reg++)
      {
        int n_cell_reg = cell_per_reg[reg];
        for(int grp = 0; grp < m_groups ; grp++)
        {        
          /// assume isotropic planck emission
          if(m_groups > 1)
          {          
            iso_emission = materials.get_mf_planck(temp_in_reg[reg], grp);
          }
          else
          {
            iso_emission = materials.get_grey_planck(temp_in_reg[reg]);
          }          
          iso_emission /= ang_quad.get_sum_w();
          
          for(int dir=0; dir < m_dir ; dir++)
          {
            for(int cell = 0; cell < n_cell_reg ; cell++)
            {
              set_cell_intensity( (cell+cell_cnt) , grp, dir, iso_emission);
            }
          }
        }
        cell_cnt += n_cell_reg;
      }    
    }
    else
    {
      std::cerr << "Unknown intensity initial condition type\n";
      exit(EXIT_FAILURE);
    }
  }

double Intensity_Data::get_intensity(const int el, const int cell, const int group, const int dir) const
{
  bool bad_input = intensity_range_check(el,cell,group,dir);
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
  
  return m_i[val_loc];
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

void Intensity_Data::set_cell_intensity(const int cell,
  const int group, const int dir, const double val) 
{  
  int loc = intensity_data_locator(0,cell,group,dir);
  bool bad_location = intensity_bounds_check(loc);
  
  if(bad_location)
  {
    std::cerr << "Error.  Trying to write to an intensity location out of range\n";
    exit(EXIT_FAILURE);
  }
  
  for(int i=0; i<m_el_per_cell ; i++)
    m_i[loc+i] = val;
    
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


bool Intensity_Data::intensity_bounds_check(const int loc) const
{
  bool is_bad_loc = false;
  if( (loc < 0) || (loc >= m_i_length) )
    is_bad_loc = true;  
  
  return is_bad_loc;
}






