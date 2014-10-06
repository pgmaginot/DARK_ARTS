#include "V_Intensity_Update.h"

V_Intensity_Update::V_Intensity_Update(const Input_Reader& input_reader, const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* materials, 
  Angular_Quadrature& angular_quadrature, const int n_stages, 
    const Temperature_Data* const t_old, const Temperature_Data* const t_star, 
    const Intensity_Data* const i_old,
    const K_Temperature* const kt, const K_Intensity* const ki)
{
  /// declare WGRS
   WG_SOLVE_TYPE solver_type = input_reader.get_within_group_solve_type();
   if(solver_type == FP_SWEEPS)
   {
    m_within_group_radiation_solver = std::shared_ptr<V_WGRS> (new WGRS_FP_Sweeps( input_reader, fem_quadrature, 
      cell_data, materials, angular_quadrature, n_stages,t_old, t_star, i_old, kt, ki));
   }
   else
   {
    std::cerr << "No other within group radiation solver coded\n";
    exit(EXIT_FAILURE);
   }
}

void V_Intensity_Update::set_time_data( const double dt, const double time_stage, const std::vector<double>& rk_a_of_stage_i, const int stage )
{
  m_within_group_radiation_solver->set_time_data(dt,time_stage,rk_a_of_stage_i,stage);
  return;
}