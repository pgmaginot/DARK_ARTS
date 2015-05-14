#ifndef V_MF_ARD_Solver_h
#define V_MF_ARD_Solver_h

#include "Input_Reader.h"
#include "Fem_Quadrature.h"
#include "Cell_Data.h"
#include "Angular_Quadrature.h"
#include <memory>
#include "V_WGRS.h"

class V_MF_ARD_Solver
{
public:
  V_MF_ARD_Solver(std::shared_ptr<V_WGRS> wgrs, const Input_Reader& input_reader, const Fem_Quadrature fem_quadrature, 
    const Cell_Data& cell_data, const Angular_Quadrature& angular_quadrature) ;
  virtual ~V_MF_ARD_Solver(){}
  
  virtual int solve_ard_problem(Intensity_Moment_Data& ard_phi_new) = 0;
  
protected:
  std::shared_ptr<V_WGRS> m_wgrs;
  const double m_ard_phi_tolerance;
  
  std::shared_ptr<V_Phi_Error_Calculator> m_convergence_calculator;
};

#endif