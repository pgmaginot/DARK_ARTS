#include "V_Intensity_Update.h"

V_Intensity_Update::V_Intensity_Update(const Input_Reader& input_reader, 
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
{
  /// declare WGRS
   WG_SOLVE_TYPE solver_type = input_reader.get_within_group_solve_type();
   if(solver_type == FP_SWEEPS)
   {
    m_within_group_radiation_solver = std::shared_ptr<V_WGRS> (new WGRS_FP_Sweeps( input_reader, fem_quadrature, 
      cell_data, materials, angular_quadrature, n_stages,t_old, i_old, kt, ki, t_star, phi_ref_norm));
   }
   else
   {
    throw Dark_Arts_Exception( SUPPORT_OBJECT ,"Desired within group radiation solver not coded");
   }
}

void V_Intensity_Update::set_time_data( const double dt, const double time_stage, const std::vector<double>& rk_a_of_stage_i, const int stage )
{
  m_within_group_radiation_solver->set_time_data(dt,time_stage,rk_a_of_stage_i,stage);
  return;
}

void V_Intensity_Update::calculate_k_i(K_Intensity& k_i, Intensity_Moment_Data& ard_phi)
{
  /// every V_WGRS object has a transport sweep object, and that's how we calculate k_i
  m_within_group_radiation_solver->sweep_for_k_i(k_i,ard_phi);
  return;
}