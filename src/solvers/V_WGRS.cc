#include "V_WGRS.h"

V_WGRS::V_WGRS(const Input_Reader& input_reader,
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
  m_transport_sweep(fem_quadrature, cell_data, materials,angular_quadrature,n_stages, t_old, i_old, kt, ki,t_star,input_reader) 
{
  switch(input_reader.get_phi_norm_type() )
  {
    case L1:
    {
      m_convergence_calculator = std::make_shared<L1_Phi_Error_Calculator>(fem_quadrature, cell_data, angular_quadrature);
      break;
    }
    case L1_RHO:
    {
      m_convergence_calculator = std::make_shared<L1_Rho_Phi_Error_Calculator>(fem_quadrature, cell_data, angular_quadrature);
      break;
    } 
    case POINTWISE:
    {
      m_convergence_calculator = std::make_shared<Pointwise_Phi_Error_Calculator>(fem_quadrature, cell_data, angular_quadrature);
      break;   
    }
    case INVALID_ERR_NORM_TYPE:
    {
      break;
    }
  }
}

void V_WGRS::set_time_data( const double dt, const double time_stage, const std::vector<double>& rk_a_of_stage_i, const int stage )
{
  m_transport_sweep.set_time_data(dt,time_stage,rk_a_of_stage_i, stage);
  
  return;
}

void V_WGRS::set_ard_phi_ptr(Intensity_Moment_Data* ard_phi_ptr)
{
  m_transport_sweep.set_ard_phi_ptr(ard_phi_ptr);
  return;
}

void V_WGRS::sweep_for_k_i(K_Intensity& k_i, Intensity_Moment_Data& ard_phi)
{
  /// set other indicators (possibly some matrix creators) to note that we are getting k_I and not sweeping
  /// ---- done in transport_sweep::sweep_mesh() when Solution_Saver is changed
  
  /// sweep the mesh and form k_I on the fly
  /// since we only have accees to one phi at this point (ard_phi, transport sweep needs to be modified have non co
  m_transport_sweep.set_sweep_type(false,true);
  
  
  m_transport_sweep.sweep_mesh(ard_phi, ard_phi);
  
  return;
}