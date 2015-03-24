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
  m_n_groups( angular_quadrature.get_number_of_groups()  ),
  m_wg_tolerance( input_reader.get_within_group_solve_tolerance()  ),
  m_max_sweeps( input_reader.get_max_number_sweeps() ),
  /// initialize to zero
  m_phi_old(cell_data, angular_quadrature,fem_quadrature,phi_ref_norm)
{
 
}

void WGRS_FP_Sweeps::kill_petsc_objects()
{
  return;
}

int WGRS_FP_Sweeps::solve(Intensity_Moment_Data& phi_new, bool& update_success)
{
  /** phi_new is what we will return
   * m_phi_old is a variable local to WGRS_FP_SWEEPS, we need to clear phi_old to start
   */
    
   /** perform sweeps until convergence is reaced */   
  int n_sweep = 0;
  for(n_sweep = 0; n_sweep < m_max_sweeps ; n_sweep++)
  {
    /// perform a transport sweep, updating phi_new
    const bool is_k_i = false;
    const bool is_krylov = false;
    
    m_transport_sweep.set_sweep_type(is_k_i,is_krylov);
    /// on first sweep, m_phi_old is a zero vector
    m_transport_sweep.sweep_mesh(m_phi_old , phi_new);
    
    /// calculate the error, normalize to phi_new
    m_err.clear();
    phi_new.normalized_difference(m_phi_old,m_err);    
    if( (m_err.get_worst_err() < m_wg_tolerance) )
    {
      /// converged !  stop the iteration
      
      break;
    }
    else
    {
      m_phi_old = phi_new;
    }
  }
  
  if( (m_err.get_worst_err() < m_wg_tolerance) )
    update_success = true;
  else
    update_success = false;
  
  
  return n_sweep;
}
