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
  m_sum_dfem_wts(fem_quad.get_sum_of_dfem_interpolation_weights()),
  m_n_dir(ang_quad.get_number_of_dir() ),
  m_phi(m_phi_length,0.)
  {  
    fem_quad.get_dfem_interpolation_point_weights(m_dfem_wts);
    m_sn_wts.resize(m_n_dir,0.);
    m_leg_moment_coeff_build.resize(m_n_dir*m_leg,0.);
    int cnt = 0;
    for(int dir = 0; dir < m_n_dir ; dir++)
    {
      m_sn_wts[dir] = ang_quad.get_w(dir);
      for(int l = 0; l< m_leg; l++)
      {
        m_leg_moment_coeff_build[cnt] = ang_quad.get_leg_moment_coeff_build(dir,l);
        cnt++;
      }
    }    
  }
    
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
  m_sum_dfem_wts(fem_quad.get_sum_of_dfem_interpolation_weights()),
  m_n_dir(ang_quad.get_number_of_dir() ),
  m_phi(m_phi_length,0.)
  { 
    fem_quad.get_dfem_interpolation_point_weights(m_dfem_wts);
    m_sn_wts.resize(m_n_dir,0.);
    m_leg_moment_coeff_build.resize(m_n_dir*m_leg,0.);
    int cnt = 0;
    for(int dir = 0; dir < m_n_dir ; dir++)
    {
      m_sn_wts[dir] = ang_quad.get_w(dir);
      for(int l = 0; l< m_leg; l++)
      {
        m_leg_moment_coeff_build[cnt] = ang_quad.get_leg_moment_coeff_build(dir,l);
        cnt++;
      }
    }    
    
    /// initialize phi to whatever i_old would dictate, then calculate norm (on the fly)
    update_phi_and_norms(i_old);
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
  m_sum_dfem_wts(intensity_moment.m_sum_dfem_wts),
  m_n_dir(intensity_moment.m_n_dir),
  m_sn_wts(intensity_moment.m_sn_wts),
  m_leg_moment_coeff_build(intensity_moment.m_leg_moment_coeff_build),  
  m_dfem_wts(intensity_moment.m_dfem_wts),
  m_phi(intensity_moment.m_phi)
{    
}

void Intensity_Moment_Data::mms_cheat(const double time_stage, const Cell_Data& cell_data, const std::vector<double> dfem_interp_points, const Input_Reader& input_reader, const Angular_Quadrature& angular_quadrature)
{
  MMS_Intensity mms( input_reader, angular_quadrature);
  Eigen::VectorXd phi_loc = Eigen::VectorXd::Zero(m_el_per_cell);
  for(int c  = 0 ; c < m_cells ; c++)
  {
    double dx = cell_data.get_cell_width(c);
    double xL = cell_data.get_cell_left_edge(c);
    for(int el = 0 ; el < m_el_per_cell ; el++)
    {
      double x = xL +dx/2.*(1. + dfem_interp_points[el]);
      phi_loc(el) = mms.get_mms_phi(x,time_stage);
    }
    set_cell_angle_integrated_intensity(c,0,0,phi_loc);
  }
}

Intensity_Moment_Data& Intensity_Moment_Data::operator= (const Intensity_Moment_Data& intensity_moment)
{
  if( (m_phi_length != intensity_moment.m_phi_length) ||
      (m_groups != intensity_moment.m_groups) ||
      (m_cells != intensity_moment.m_cells) ||
      (m_leg != intensity_moment.m_leg)       
    )
    throw Dark_Arts_Exception( VARIABLE_STORAGE , "Trying to assign different sized Intensity_Moment_Data object");
  
  for(int i=0; i<m_phi_length; i++)
    m_phi[i] = intensity_moment.m_phi[i];
    
  return *this;
}
  
void Intensity_Moment_Data::get_cell_angle_integrated_intensity(const int cell,
  const int group, const int l_mom, Eigen::VectorXd& loc_phi_vec) const
{
  int val_loc = angle_integrated_data_locator(0,cell,group,l_mom);
  
  for(int i=0; i< m_el_per_cell; i++)
    loc_phi_vec(i) = m_phi[val_loc+i];

  return;
}

void Intensity_Moment_Data::set_cell_angle_integrated_intensity(const int cell,
    const int group, const int l_mom, const Eigen::VectorXd& val) 
{
  int val_loc = angle_integrated_data_locator(0,cell,group,l_mom);
  
  for(int i=0; i< m_el_per_cell ; i++)
    m_phi[val_loc+i] = val(i);
    
  return;
}

/// take in a contribution to m_phi, contrib, and add to the correct elements
void Intensity_Moment_Data::add_contribution(const int cell, const int grp, const int l_mom, Eigen::VectorXd& contrib)
{
  int val_loc = angle_integrated_data_locator(0,cell,grp,l_mom);
  
  for(int i=0; i< m_el_per_cell ; i++)
    m_phi[val_loc+i] += contrib(i);
    
  return;
}

/// This function controls the layout of angle_integrated intensities in memory!!
int Intensity_Moment_Data::angle_integrated_data_locator(const int el, const int cell, const int group, const int l_mom) const
{
  int loc = -1;
  
  angle_integrated_range_check(el,cell,group,l_mom);
  
  /**
    Arrange data as:
    element
    moment
    group
    cell
  */
  loc = el + l_mom*m_el_per_cell + group*m_el_times_l_mom + cell*m_el_times_l_mom_times_group;
  
  angle_integrated_bounds_check(loc);
  
  return loc;
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
    
  if(is_bad)
  {
    std::stringstream err;
    err <<" Attemping to access intensity_moment_data outside of logical bounds.\n" ;
    err << "Requested Element: " << el << " of cell: " << cell << " group: " << grp << " Legendre moment: " << l_mom << std::endl;
    err << "Max Element: " << m_el_per_cell << "Max Cell: " << m_cells  << " Max group: " << m_groups << " Max L_mom " << m_leg<< std::endl;
    throw Dark_Arts_Exception( VARIABLE_STORAGE , err.str());
  }
  
  return is_bad;
}

bool Intensity_Moment_Data::angle_integrated_bounds_check(const int loc) const
{
  bool is_bad_loc = false;
  if( (loc < 0) || (loc >= m_phi_length) )
    is_bad_loc = true;  
    
  if(is_bad_loc)
  {
    std::stringstream err;
    err << "Calculated memory idex is out of range.  Calculated: " << loc << " . m_phi_length is: " << m_phi_length ; 
    throw Dark_Arts_Exception(VARIABLE_STORAGE , err.str() );
  }
  
  return is_bad_loc;
}

void Intensity_Moment_Data::clear_angle_integrated_intensity(void)
{
  for(int i=0;i<m_phi_length;i++)
    m_phi[i] = 0.;
    
  return;
}

void Intensity_Moment_Data::normalized_difference(const Intensity_Moment_Data& phi_compare, Err_Phi& err_phi) const
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
  
  Eigen::VectorXd phi_this = Eigen::VectorXd::Zero(m_el_per_cell);
  Eigen::VectorXd phi_other =  Eigen::VectorXd::Zero(m_el_per_cell);
  double loc_err;
  for(int cell = 0; cell < m_cells; cell++)
  {
    for(int grp = 0; grp < m_groups ; grp++)
    {
      for(int l_mom = 0; l_mom < m_leg ; l_mom++)
      {
        get_cell_angle_integrated_intensity(cell,grp,l_mom,phi_this);
        phi_compare.get_cell_angle_integrated_intensity(cell,grp,l_mom,phi_other);
        
        for(int el = 0; el < m_el_per_cell ; el++)
        {
          if(fabs(phi_this(el)) > m_norm_for_err[grp])
          {
            /// not dividing by a near zero
            loc_err = fabs( (phi_this(el) - phi_other(el) )/phi_this(el) );
          }
          else
          {
            loc_err = fabs( phi_this(el) - phi_other(el) );
          }
          
          if(loc_err > err_phi.get_worst_err() )
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

void Intensity_Moment_Data::update_phi_and_norms(const Intensity_Data& i_old)
{  
  Eigen::VectorXd i_local(m_el_per_cell);  
  Eigen::VectorXd phi_contrib_local(m_el_per_cell);  
  
  for(int cell = 0; cell < m_cells; cell++)
  {
    for(int grp = 0; grp < m_groups ; grp++)
    {
      for(int dir = 0; dir < m_n_dir ; dir++)
      {
        i_old.get_cell_intensity(cell,grp,dir,i_local);  
        /// calculate numerical average of scalar flux in this cell, and keep a running tally
        for(int el = 0; el < m_el_per_cell ; el++)
        {
          m_norm_for_err[grp] += m_dfem_wts[el]*m_sn_wts[dir]*i_local(el);
        }
        /// now save flux moments of solution
        for(int l = 0; l< m_leg; l++)
        {
          phi_contrib_local = m_sn_wts[dir]*m_leg_moment_coeff_build[m_leg*dir + l]* i_local;
          add_contribution(cell,grp,l,phi_contrib_local );
        }        
      }
    }
  }
  
  const double div = m_small_ratio/( m_sum_dfem_wts * double(m_cells) );
  /// numerical average of phi across all cells
  for(int g=0; g<m_groups; g++)
    m_norm_for_err[g] = fabs(m_norm_for_err[g])*div;

  return;
}

void Intensity_Moment_Data::add_from_array_pointer( double *update_array , const int cell, const int group)
{
  int val_loc = angle_integrated_data_locator(0,cell,group,0);
  
  for(int i=0; i< m_el_per_cell ; i++)
    m_phi[val_loc+i] += update_array[i];
    
  return;
}