/** @file   Intensity_Moment_Data.cc
  *   @author pmaginot
  *   @brief Implement the Intensity_Moment_Data class
  *   Store angle integrated group intensity (phi)
*/
#include "Intensity_Moment_Data.h"

Intensity_Moment_Data::Intensity_Moment_Data(const Cell_Data& cell_data, const Angular_Quadrature& ang_quad,
  const Fem_Quadrature& fem_quad, const std::vector<double>& reference_phi_norm)
  : 
  m_cells{cell_data.get_total_number_of_cells() } , 
  m_groups{ang_quad.get_number_of_groups() } , 
  m_leg{ ang_quad.get_number_of_leg_moments() }, 
  m_el_per_cell{fem_quad.get_number_of_interpolation_points() }  ,
  m_el_times_l_mom{m_leg*m_el_per_cell},
  m_el_times_l_mom_times_group{m_el_times_l_mom*m_groups},   
  m_phi_length{m_cells*m_groups*m_leg*m_el_per_cell},
  m_norm_for_err(reference_phi_norm),
  m_small_ratio{1.0E-6},
  m_phi(m_phi_length,0.)
  {    }
    
/// Initializer used for ard_phi within Time_Marcher.  Initialize from I_old.
Intensity_Moment_Data::Intensity_Moment_Data(const Cell_Data& cell_data, const Angular_Quadrature& ang_quad,
  const Fem_Quadrature& fem_quad, const Intensity_Data& i_old)
  : 
  m_cells{cell_data.get_total_number_of_cells() } , 
  m_groups{ang_quad.get_number_of_groups() } , 
  m_leg{ ang_quad.get_number_of_leg_moments() }, 
  m_el_per_cell{fem_quad.get_number_of_interpolation_points() }  ,
  m_el_times_l_mom{m_leg*m_el_per_cell},
  m_el_times_l_mom_times_group{m_el_times_l_mom*m_groups},   
  m_phi_length{m_cells*m_groups*m_leg*m_el_per_cell},
  m_norm_for_err(m_groups,0.),
  m_small_ratio{1.0E-6},
  m_phi(m_phi_length,0.)
  { 
    /// calculate 0th angular moment of each group of i_old to estimate scale for phi error comparison
    calculate_reference_phi_norms(i_old, ang_quad, fem_quad);
  }
  
Intensity_Moment_Data::Intensity_Moment_Data(const Intensity_Moment_Data& intensity_moment)
  :
  m_cells{ intensity_moment.m_cells } , 
  m_groups{intensity_moment.m_groups } , 
  m_leg{ intensity_moment.m_leg }, 
  m_el_per_cell{intensity_moment.m_el_per_cell }  ,
  m_el_times_l_mom{intensity_moment.m_el_times_l_mom},
  m_el_times_l_mom_times_group{intensity_moment.m_el_times_l_mom_times_group},   
  m_phi_length{intensity_moment.m_phi_length},
  m_norm_for_err(intensity_moment.m_norm_for_err),
  m_small_ratio{intensity_moment.m_small_ratio},  
  m_phi(intensity_moment.m_phi)
{    
}

Intensity_Moment_Data& Intensity_Moment_Data::operator= (const Intensity_Moment_Data& intensity_moment)
{
  if( (m_phi_length != intensity_moment.m_phi_length) ||
      (m_groups != intensity_moment.m_groups) ||
      (m_cells != intensity_moment.m_cells) ||
      (m_leg != intensity_moment.m_leg) 
  )
  {
    std::cerr << "Error assigning to Intensity_Moment_Data object, objects not the same size\n";
    exit(EXIT_FAILURE);
  }
  
  for(int i=0; i<m_phi_length; i++)
    m_phi[i] = intensity_moment.m_phi[i];
    
  return *this;
}

double Intensity_Moment_Data::get_angle_integrated_intensity(const int el, const int cell,
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


  
void Intensity_Moment_Data::get_cell_angle_integrated_intensity(const int cell,
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

void Intensity_Moment_Data::set_cell_angle_integrated_intensity(const int cell,
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

void Intensity_Moment_Data::add_contribution(const int cell, const int grp, const int l_mom, Eigen::VectorXd& contrib)
{
  bool bad_input = angle_integrated_range_check(0,cell,grp,l_mom);
  if(bad_input)
  {
    std::cerr << "Error.  Attempting to set out of logical range angle integrated intensity\n";
    exit(EXIT_FAILURE);
  }
  
  int val_loc = angle_integrated_data_locator(0,cell,grp,l_mom);
  
  for(int i=0; i< m_el_per_cell ; i++)
    m_phi[val_loc+i] += contrib(i);
    
  return;
}


bool Intensity_Moment_Data::angle_integrated_range_check(const int el, const int cell, 
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

/// This function controls the layout of angle_integrated intensities in memory!!
int Intensity_Moment_Data::angle_integrated_data_locator(const int el, const int cell, const int group, const int l_mom) const
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


bool Intensity_Moment_Data::angle_integrated_bounds_check(const int loc) const
{
  bool is_bad_loc = false;
  if( (loc < 0) || (loc >= m_phi_length) )
    is_bad_loc = true;  
  
  return is_bad_loc;
}

void Intensity_Moment_Data::clear_angle_integrated_intensity(void)
{
  for(int i=0;i<m_phi_length;i++)
    m_phi[i] = 0.;
    
  return;
}

void Intensity_Moment_Data::normalized_difference(Intensity_Moment_Data& phi_compare, Err_Phi& err_phi) const
{
  /**
    Find the maximum, normalized difference between this iterate and another Intensity_Moment_Data object
    
    tol = small_num*\norm{\phi(t=0) 
    norm = avg?
    
    if( phi < tol)
    {
      loc_err = abs( phi_compare - phi)
    }
    else
    {
      loc_err = abs( (phi_compare - phi)/phi ) 
    }
    
  */
  double phi_new_this = 0.;
  double loc_err;
  for(int cell = 0; cell < m_cells; cell++)
  {
    for(int grp = 0; grp < m_groups ; grp++)
    {
      for(int l_mom = 0; l_mom < m_leg ; l_mom++)
      {
        for(int el = 0; el < m_el_per_cell ; el++)
        {
          phi_new_this = get_angle_integrated_intensity(el,cell,grp,l_mom);
          if(fabs(phi_new_this) > m_norm_for_err[grp])
          {
            /// will not be normalizing to zero
            loc_err = fabs( (phi_new_this - phi_compare.get_angle_integrated_intensity(el,cell,grp,l_mom) )/phi_new_this);
          }
          else
          {
            loc_err = fabs( phi_new_this - phi_compare.get_angle_integrated_intensity(el,cell,grp,l_mom) );
          }
          
          if(loc_err < err_phi.get_worst_err() )
          {
            err_phi.set_error(cell, grp, l_mom, loc_err);
          }          
        }
      }
    }
  }
  
  

  /// Err_Phi.err , Err_Phi.el_num , Err_Phi.cell_num, Err_Phi.group_num, Err_Phi.l_mom_num
}

void Intensity_Moment_Data::get_all_moments(
  std::vector<Eigen::VectorXd>& local_phi, const int cell, const int grp) const
{
  for(int l=0; l<m_leg; l++)
  {
    get_cell_angle_integrated_intensity(cell, grp, l, local_phi[l]);
  }
  
  return;
}

void Intensity_Moment_Data::get_phi_norm(std::vector<double>& norm_vec) const
{
  norm_vec = m_norm_for_err;
  return;
}

void Intensity_Moment_Data::calculate_reference_phi_norms(const Intensity_Data& i_old, 
  const Angular_Quadrature& ang_quad, 
  const Fem_Quadrature& fem_quad)
{
  std::vector<double> w_dfem_points;
  fem_quad.get_dfem_interpolation_point_weights(w_dfem_points);
  
  for(int cell = 0; cell < m_cells; cell++)
  {
    for(int grp = 0; grp < m_groups ; grp++)
    {
      for(int dir = 0; dir < ang_quad.get_number_of_dir() ; dir++)
      {
        for(int el = 0; el < m_el_per_cell ; el++)
        {
          m_norm_for_err[grp] += w_dfem_points[el]/2.*i_old.get_intensity(el,cell,grp,dir);
        }
      }
    }
  }
  
  for(int g=0; g<m_groups; g++)
    m_norm_for_err[g] = fabs(m_norm_for_err[g])*m_small_ratio/( double(m_cells) );

  return;
}