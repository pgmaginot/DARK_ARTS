#include "MF_ARD_Solver_FP_No_LMFGA.h"


MF_ARD_Solver_FP_No_LMFGA::MF_ARD_Solver_FP_No_LMFGA(const Input_Reader& input_reader, 
    const Fem_Quadrature& fem_quad, 
    const Cell_Data& cell_data, 
    const Angular_Quadrature& ang_quad, 
    std::shared_ptr<V_WGRS> wgrs, 
    std::vector<double>&  phi_ref_norm)
:
V_MF_ARD_Solver(wgrs, input_reader),
m_max_iterations{ input_reader.get_max_ard_iterations() },
m_ard_old(cell_data, ang_quad, fem_quad, phi_ref_norm )
{
  
}
  
int MF_ARD_Solver_FP_No_LMFGA::solve_ard_problem(Intensity_Moment_Data& ard_phi_new)
{
  int inners = 0;
  m_ard_old.clear_angle_integrated_intensity();
  for(int iter = 0; iter < m_max_iterations ; iter++)
  {
    /// set which phi we will be using for ard
    m_wgrs->set_ard_phi_ptr(&m_ard_old);
    /// get a new ard_phi based on the old value
    bool wg_success = false;
    inners += m_wgrs->solve(ard_phi_new,wg_success);
    /// get normalized change
    ard_phi_new.normalized_difference(m_ard_old,m_ard_err);
    /// check if change indciates convergence
  }
  return inners;
}
