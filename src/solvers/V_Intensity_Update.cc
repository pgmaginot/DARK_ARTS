#include "V_Intensity_Update.h"

V_Intensity_Update::V_Intensity_Update(const Input_Reader& input_reader, const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* materials, 
  const Angular_Quadrature& angular_quadrature, const int n_stages)
    :
    m_dt{-1.},
    m_stage{-1}
{
  /// declare WGRS
   WG_SOLVE_TYPE solver_type = input_reader.get_within_group_solve_type();
   if(solver_type == FP_SWEEPS)
   {
    m_within_group_radiation_solver = std::shared_ptr<V_WGRS> (new WGRS_FP_Sweeps( input_reader, fem_quadrature, 
      cell_data, materials, angular_quadrature, n_stages));
   }
   else
   {
    std::cerr << "No other within group radiation solver coded\n";
    exit(EXIT_FAILURE);
   }
}

