#ifndef V_MF_ARD_Solver_h
#define V_MF_ARD_Solver_h

#include "Input_Reader.h"
#include <memory>
#include "V_WGRS.h"

class V_MF_ARD_Solver
{
public:
  V_MF_ARD_Solver(std::shared_ptr<V_WGRS> wgrs, const Input_Reader& input_reader) ;
  virtual ~V_MF_ARD_Solver(){}
  
  virtual void solve_ard_problem(Intensity_Moment_Data& ard_phi_new) = 0;
  
protected:
  std::shared_ptr<V_WGRS> m_wgrs;
  const double m_ard_phi_tolerance;
};

#endif