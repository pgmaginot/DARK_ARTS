#include "WGRS_FP_Sweeps.h"

WGRS_FP_Sweeps::WGRS_FP_Sweeps(const Input_Reader& input_reader, 
  const Fem_Quadrature& fem_quadrature, 
  const Cell_Data& cell_data, 
  Materials& materials, 
  const Angular_Quadrature& angular_quadrature, 
  const int n_stages, 
  const Temperature_Data& t_old, 
  const Intensity_Data& i_old,
  const K_Temperature& kt, 
  K_Intensity& ki,
  const Temperature_Data& t_star,
  std::vector<double>& phi_ref_norm)
  : 
  V_WGRS(input_reader, fem_quadrature,cell_data, materials,angular_quadrature, n_stages,t_old, i_old, kt, ki,t_star, phi_ref_norm),
  m_n_groups{ angular_quadrature.get_number_of_groups()  },
  m_wg_tolerance{ input_reader.get_within_group_solve_tolerance()  },
  m_max_sweeps{ input_reader.get_max_number_sweeps() },
  /// initialize to zero
  m_phi_old(cell_data, angular_quadrature,fem_quadrature,phi_ref_norm)
{
 
}
void WGRS_FP_Sweeps::solve(Intensity_Moment_Data& phi_new)
{
  /** phi_new is what we will return
   * m_phi_old is a variable local to WGRS_FP_SWEEPS, we need to clear phi_old to start
   */
  m_phi_old.clear_angle_integrated_intensity();
    
   /** perform sweeps until convergence is reaced */
  for(int n_sweep = 0; n_sweep < m_max_sweeps ; n_sweep++)
  {
    /// perform a transport sweep, updating phi_new
    m_transport_sweep.set_sweep_type(false,false);
    m_transport_sweep.sweep_mesh(m_phi_old , phi_new);
    
    /// clear out the error object
    m_err.clear();
    
    /// calculate the error, normalize to phi_new
    phi_new.normalized_difference(m_phi_old,m_err);
    
    if( m_err.get_worst_err() < m_wg_tolerance )
    {
      /// converged !  stop the iteration
      break;
    }
    /// we have redefined the = operator to only copy over the values of flux moments
    m_phi_old = phi_new;
  }
  
  return;
}


void WGRS_FP_Sweeps::set_time_data( const double dt, const double time_stage, const std::vector<double>& rk_a_of_stage_i, const int stage )
{
  m_transport_sweep.set_time_data(dt,time_stage,rk_a_of_stage_i, stage);
  
  return;
}