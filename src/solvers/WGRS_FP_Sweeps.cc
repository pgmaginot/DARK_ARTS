#include "WGRS_FP_Sweeps.h"

WGRS_FP_Sweeps::WGRS_FP_Sweeps(const Input_Reader& input_reader, const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* materials, 
    const Angular_Quadrature& angular_quadrature, const int n_stages)
  : 
  V_WGRS(input_reader, fem_quadrature,cell_data, materials,angular_quadrature, n_stages),
  m_n_groups{ angular_quadrature.get_number_of_groups()  },
  m_wg_tolerance{ input_reader.get_within_group_solve_tolerance()  }
{

}


void WGRS_FP_Sweeps::solve(const Temperature_Data& t_star, Intensity_Moment_Data& phi)
{
  return;
}
