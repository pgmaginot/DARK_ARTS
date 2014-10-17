
#include "V_MF_ARD_Solver.h"

V_MF_ARD_Solver::V_MF_ARD_Solver(std::shared_ptr<V_WGRS> wgrs, const Input_Reader& input_reader)  
:
m_wgrs{wgrs},
m_ard_phi_tolerance{ input_reader.get_between_group_solve_tolerance() }
{

}