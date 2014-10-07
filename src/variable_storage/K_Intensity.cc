/** @file   K_Intensity.cc
  *   @author pmaginot
  *   @brief Implement the K_Intensity class
  *   Store k quantities for DIRK time integration, for intensity data (no angular moment storage)
*/
#include "K_Intensity.h"

K_Intensity::K_Intensity(const int n_cells, const int n_stages, 
    const Fem_Quadrature& fem_quadrature, const Angular_Quadrature& angular_quadrature)
  /// initilaize range members
  :
  m_cells{ n_stages },
  m_el_per_cell{ fem_quadrature.get_number_of_interpolation_points()},
  m_n_stages{n_stages},
  m_n_dir{ angular_quadrature.get_number_of_dir() },
  m_n_grp{ angular_quadrature.get_number_of_groups() },
  m_k_length{m_n_grp*m_n_dir*m_n_stages*m_el_per_cell*m_cells},
  m_el_stage{ m_el_per_cell*m_n_stages},
  m_el_stage_dir_div_2{ m_el_stage*m_n_dir/2},
  m_el_stage_dir_div_2_grp{m_el_stage_dir_div_2*m_n_grp},
  m_offset{m_n_dir/2*m_n_grp*m_cells*m_n_stages},
  m_dir_div_2{m_n_dir/2},
  m_vec_sum{ Eigen::VectorXd::Zero(m_el_per_cell) },
  m_vec_retrieve{ Eigen::VectorXd::Zero(m_el_per_cell) }  
{  
  m_ki.resize(m_k_length,0.) ;
  m_rk_b.resize(m_el_per_cell,0.);
}
/// Public accessor functions
void K_Intensity::get_ki(const int cell, const int grp, const int dir, const int stage, Eigen::VectorXd& ki) const
{  
  /// find first element location
  int base_loc = ki_data_locator(cell, grp, dir, stage);
  for(int i=0;i<m_el_per_cell;i++)
    ki(i) = m_ki[base_loc + i];
    
  return;
}

/// Public functions to save values
void K_Intensity::set_ki(const int cell, const int grp, const int dir, const int stage, Eigen::VectorXd& ki)
{
  /// find first element location
  int base_loc = ki_data_locator(cell, grp, dir, stage);
  for(int i=0;i<m_el_per_cell;i++)
    m_ki[base_loc + i]=ki(i);   
   
  return ;
}
bool K_Intensity::ki_range_check(const int cell, const int grp, const int dir,const int stage) const
{
  bool is_bad = false;
  
  if( (stage >= m_n_stages ) || (stage < 0) )
    is_bad = true;
       
  if( (cell < 0) || (cell >= m_cells) )
    is_bad = true;
    
  if( (grp < 0 ) || (grp >= m_n_grp) )
    is_bad = true;
    
  if( (dir < 0) || ( dir > m_n_dir) ) 
    is_bad = true;
  
  return is_bad;
}
  
/**
  This is the function that controls k_I data layout!!
  Layout is critical to memory performance!
*/
int K_Intensity::ki_data_locator(const int cell, const int grp, const int dir,const int stage) const
{
  
  if( ki_range_check(cell,grp,dir,stage) )
  {
    std::cerr << "Attempting to access illogical k_intensity location\n";
    exit(EXIT_FAILURE);
  }
  
  /** we are sweeping as follows:   
      mu < 0 
      for(cell=1end...1)
        for(grp=1 ... g_max)
          for(d=1 ... n_dir/2)  
          
      mu > 0 
      for(cell=1...end)
        for(grp=1 ... g_max)
          for(d=n_dir/2+1 >>> n_dir)
          
      -----> RK stages will be the innermost loop!      
  */
  int loc_val = 0;
  if( dir < m_dir_div_2)
  {
    loc_val = stage*m_el_per_cell + dir*m_el_stage + grp*m_el_stage_dir_div_2 + (m_cells-cell-1)*m_el_stage_dir_div_2_grp;
  }
  else
  {
    loc_val = m_offset + stage*m_el_per_cell + (dir-m_dir_div_2)*m_el_stage + grp*m_el_stage_dir_div_2 + cell*m_el_stage_dir_div_2_grp;;
  }
  
  if( ki_bounds_check(loc_val) )
  {
    std::cerr << "Location out of bounds in k_intensity data\n";
    exit(EXIT_FAILURE);
  }
  
  return loc_val;
}
   
bool K_Intensity::ki_bounds_check(const int loc) const
{
  bool is_bad_loc = false;
  // if( (loc < 0) || (loc >= m_k_length) )
    // is_bad_loc = true;  
  
  return is_bad_loc;
}

void K_Intensity::advance_intensity(Intensity_Data& i_old, const double dt, const Time_Data* time_data)
{
  int start = 0;
  int end = m_cells;
  int incr = -1;
  int offset = 0;
  int dir = 0;
  time_data->get_b_dt_constants(m_rk_b,dt);
  for(int dir_set = 0; dir_set < 2; dir_set++)
  {
    if(dir_set ==0)
    {
      /// negative mu
      start = m_cells;
      end = -1;
      incr = -1;
      offset = 0;
    }
    else if(dir_set ==1)
    {
      /// positive mu
      start = 0;
      end = m_cells;
      incr = 1;
      offset = m_n_dir/2;
    }
    
    for(int c = start ; c != end ; c += incr)
    {
      for(int g = 0 ; g <= m_n_grp ; g++)
      {
        for(int d = 0; d < m_dir_div_2 ; d++)
        {
          dir = offset + d;
          m_vec_sum = Eigen::VectorXd::Zero(m_el_per_cell) ;
          for(int s = 0; s < m_n_stages ; s++)
          {
            get_ki(c,g,dir,s,m_vec_retrieve);
            m_vec_sum += m_rk_b[s]*m_vec_retrieve;
          }
          i_old.get_cell_intensity(c,g,dir,m_vec_retrieve);
          m_vec_sum += m_vec_retrieve;
          i_old.set_cell_intensity(c,g,dir,m_vec_sum);
        }
      }
    }
  }
  
  return;
}

