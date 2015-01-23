#ifndef MF_ARD_Solver_FP_No_LMFGA_h
#define MF_ARD_Solver_FP_No_LMFGA_h

#include "V_MF_ARD_Solver.h"

class MF_ARD_Solver_FP_No_LMFGA : public V_MF_ARD_Solver
{
public:
  MF_ARD_Solver_FP_No_LMFGA(const Input_Reader& input_reader, 
    const Fem_Quadrature& fem_quadrature, 
    const Cell_Data& cell_data, 
    const Angular_Quadrature& angular_quadrature, 
    std::shared_ptr<V_WGRS> wgrs, 
    std::vector<double>&  phi_ref_norm);
  
  virtual ~MF_ARD_Solver_FP_No_LMFGA(){}
  
  int solve_ard_problem(Intensity_Moment_Data& ard_phi_new) override;

protected:
  const int m_max_iterations;
  Intensity_Moment_Data m_ard_old;
  Err_Phi m_ard_err;


};

#endif