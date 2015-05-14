#include "V_MF_ARD_Solver.h"

V_MF_ARD_Solver::V_MF_ARD_Solver(std::shared_ptr<V_WGRS> wgrs, const Input_Reader& input_reader, const Fem_Quadrature fem_quadrature, 
    const Cell_Data& cell_data, const Angular_Quadrature& angular_quadrature)  
:
m_wgrs(wgrs),
m_ard_phi_tolerance( input_reader.get_between_group_solve_tolerance() )
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