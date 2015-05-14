#ifndef L1_Rho_Phi_Error_Calculatorr_h
#define L1_Rho_Phi_Error_Calculator_h

#include "V_Phi_Error_Calculator.h"
#include "Intensity_Moment_Data.h"
#include "Cell_Data.h"
#include "Fem_Quadrature.h"
#include "Input_Reader.h"

/** @file   V_Phi_Error_Calculator.h
  *   @author pmaginot
  *   @brief Base class to implement interface for various measures to assess convergence of phi: L1, L2 , rho normalized L1, rho normalized L2, pointwise convergence check
 */

class L1_Rho_Phi_Error_Calculator : public V_Phi_Error_Calculator
{
public:
  L1_Rho_Phi_Error_Calculator( const Fem_Quadrature& fem_quadrature, 
    const Cell_Data& cell_data,
    const Angular_Quadrature& angular_quadrature);
    
  virtual ~L1_Rho_Phi_Error_Calculator(){}
  
  double calculate_phi_error_norm(const Intensity_Moment_Data& phi_new , const Intensity_Moment_Data& phi_old, const int iter) override;
private:
  double m_norm_diff_old;
  double m_norm_diff_new;
};

#endif