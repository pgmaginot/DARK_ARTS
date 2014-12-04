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
  m_offset{ m_dir_div_2*m_el_per_cell*m_groups*m_cells },
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
  m_offset{ m_dir_div_2*m_el_per_cell*m_groups*m_cells },
  m_el_times_dir_div_2{ m_dir_div_2*m_el_per_cell }, 
  m_el_times_dir_div_2_times_grp{m_el_times_dir_div_2 * m_groups}, 
  m_i_length{m_cells*m_groups*m_dir*m_el_per_cell},
  m_i(m_i_length,0.)
  {    
    /// load the initial conditions    
    RADIATION_IC_TYPE ic_type = input_reader.get_radiation_ic_type();
    if( ic_type == PLANCKIAN_IC)
    {
      /// loop over regions.  Each region could have a different initial condition
      std::vector<int> cell_per_reg;
      std::vector<double> temp_in_reg;
      
      input_reader.get_cells_per_region_vector(cell_per_reg);           
      
      int cell_cnt = 0;
      double iso_emission = 0.;
      for(int reg = 0 ; reg < input_reader.get_n_regions() ; reg++)
      {
        double temp_eval = input_reader.get_region_radiation_temperature(reg);
        int n_cell_reg = cell_per_reg[reg];
        
        for(int grp = 0; grp < m_groups ; grp++)
        {        
          /// assume isotropic planck emission
          if(m_groups > 1)
          {          
            iso_emission = materials.get_mf_planck(temp_eval, grp);
          }
          else
          {
            iso_emission = materials.get_grey_planck(temp_eval);
          }          
          /// planck function returns are per steradian
          
          for(int cell = 0; cell < n_cell_reg ; cell++)
          {
            for(int dir=0; dir < m_dir ; dir++)
            {
              set_cell_intensity( cell+cell_cnt , grp, dir, iso_emission);
            }
          }
        }
        cell_cnt += n_cell_reg;
      }    
    }
    else if(ic_type == MMS_RADIATION_IC)
    {
      std::vector<int> cell_per_reg;      
      input_reader.get_cells_per_region_vector(cell_per_reg);           
      
      Eigen::VectorXd local_i;
      local_i = Eigen::VectorXd::Zero(m_el_per_cell);
      std::vector<double> dfem_loc;
      fem_quad.get_dfem_interpolation_point(dfem_loc);
      
      /// get radiation space angle and temporal components
      MMS_Intensity i_ic(input_reader,ang_quad);
      
      int cell_cnt = 0;
      const double time = input_reader.get_t_start();
      for(int reg = 0 ; reg < input_reader.get_n_regions() ; reg++)
      {
        int n_cell_reg = cell_per_reg[reg];        
        for(int c = 0; c < n_cell_reg ; c++)
        {
          int cell = c + cell_cnt;       
          double xL = cell_data.get_cell_left_edge(cell);
          double dx = cell_data.get_cell_width(cell);
          for(int dir=0; dir < m_dir ; dir++)
          {
            for(int el = 0; el < m_el_per_cell ; el++)
            {
              double position = xL + dx/2.*( 1. + dfem_loc[el] );
              local_i(el) = i_ic.get_mms_intensity(position, time, dir);
            }
            set_cell_intensity( cell , 0, dir, local_i);
          }
        }
        cell_cnt += n_cell_reg;
      }    
    }
    else
      throw Dark_Arts_Exception( VARIABLE_STORAGE , "Unknown intensity initial condition type");
  }

double Intensity_Data::get_intensity(const int el, const int cell, const int group, const int dir) const
{ 
  int val_loc = intensity_data_locator(el,cell,group,dir);
 
  return m_i[val_loc];
}

void Intensity_Data::get_cell_intensity(const int cell, const int group, 
  const int dir, Eigen::VectorXd& loc_i_vec) const
{
  int val_loc = intensity_data_locator(0,cell,group,dir);
  
  for(int i=0; i<m_el_per_cell; i++)
    loc_i_vec(i) = m_i[val_loc+i];
    
  return;
}

void Intensity_Data::set_cell_intensity(const int cell,
  const int group, const int dir, const double val) 
{  
  int loc = intensity_data_locator(0,cell,group,dir);
  
  for(int i=0; i<m_el_per_cell ; i++)
    m_i[loc+i] = val;
    
  return;
}

void Intensity_Data::set_cell_intensity(const int cell,
  const int group, const int dir, const Eigen::VectorXd& val) 
{  
  int loc = intensity_data_locator(0,cell,group,dir);
  
  for(int i=0; i<m_el_per_cell ; i++)
    m_i[loc+i] = val(i);
    
  return;
}

  /* ***************************************************
  *
  *   Protected Functions
  *
  *************************************************** */

/// This function controls the layout of intensity in memory!!
int Intensity_Data::intensity_data_locator(const int el, const int cell, const int group, const int dir) const
{
  intensity_range_check(el,cell,group,dir);
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
  
  intensity_bounds_check(loc);
  return loc;
}

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
    
  if(is_bad)
  {
    std::stringstream err;
    err <<" Attemping to access intensity outside of logical bounds. \n" ;
    err << "Requested Element: " << el << " of cell: " << cell << " group: " << grp << " direction: " << dir << std::endl;
    err << "Max Element: " << m_el_per_cell << "Max Cell: " << m_cells << " Max Dir: " << m_dir << " Max group: " << m_groups << std::endl;
    throw Dark_Arts_Exception( VARIABLE_STORAGE , err.str());
  }
  
  return is_bad;
}

bool Intensity_Data::intensity_bounds_check(const int loc) const
{
  bool is_bad_loc = false;
  if( (loc < 0) || (loc >= m_i_length) )
  {
    is_bad_loc = true;  
    std::stringstream err; 
    err << "Calculated intensity memory index outside of possible values\n";
    err << "Requested element: " << loc << " Length of m_i: " << m_i_length << std::endl;
    throw Dark_Arts_Exception( VARIABLE_STORAGE , err.str() );
  }
  
  return is_bad_loc;
}






