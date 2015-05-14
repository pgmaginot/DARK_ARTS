#ifndef V_Phi_Error_Calculator_h
#define V_Phi_Error_Calculator_h

#include "Intensity_Moment_Data.h"
#include "Cell_Data.h"
#include "Fem_Quadrature.h"
#include "Input_Reader.h"
#include <algorithm>

/** @file   V_Phi_Error_Calculator.h
  *   @author pmaginot
  *   @brief Base class to implement interface for various measures to assess convergence of phi: L1, L2 , rho normalized L1, rho normalized L2, pointwise convergence check
 */

class V_Phi_Error_Calculator
{
public:
  V_Phi_Error_Calculator( const Fem_Quadrature& fem_quadrature, 
    const Cell_Data& cell_data,
    const Angular_Quadrature& angular_quadrature);
    
  virtual ~V_Phi_Error_Calculator(){}
  
  virtual double calculate_phi_error_norm(const Intensity_Moment_Data& phi_new , const Intensity_Moment_Data& phi_old, const int iter) = 0;
protected:
  const int m_n_cell;
  const int m_n_groups;
  const int m_n_l_mom;
  const int m_n_el;
};

#endif