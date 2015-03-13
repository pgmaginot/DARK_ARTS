#include <iostream>
#include <vector>
#include <algorithm>
#include <string>

#include "MIP_Kappa_Calculator.h"
#include "Dark_Arts_Exception.h"

/**
  Goal of this unit test is to check source moment formation, and reaction matrix for SLXS Lobatto scheme
  -This is a MMS problem with a spatially varying sigma_a, constant cv , zero sig_s
*/ 

int main(int argc, char** argv)
{
  int val = 0;
  const double tol = 1.0E-8;
  
  try{
    const double dx_l = 1.2;
    const double dx_r = 2.1;
    
    const double d_1 = 1.3;
    const double d_2 = 2.9;
    
    const int p_ord = 3;
    
    const double z_mip = 4.;
    
    const double kappa_left_edge = std::max(0.25, d_1/dx_l*(4.*3*(3.+1.)));
    const double kappa_center_edge = std::max(0.25, 4./2.*(3.*(3.+1.))*d_1/dx_l + 4./2.*(3.*(3.+1.))*d_2/dx_r );
    const double kappa_right_edge = std::max(0.25, d_2/dx_r*(4.*3*(3.+1.)));;
       
    MIP_Kappa_Calculator kappa_calculator(p_ord,z_mip);    
    
    if( fabs(kappa_calculator.calculate_boundary_kappa(dx_l,d_1) - kappa_left_edge) > tol )
      throw Dark_Arts_Exception(MIP , "Not Calculating left boundary kappa correctly");
      
    if( fabs(kappa_calculator.calculate_boundary_kappa(dx_r,d_2) - kappa_right_edge) > tol )
      throw Dark_Arts_Exception(MIP , "Not Calculating right boundary kappa correctly");
      
    if( fabs(kappa_calculator.calculate_interior_edge_kappa(dx_l,dx_r,d_1,d_2) - kappa_center_edge) > tol )
      throw Dark_Arts_Exception(MIP , "Not Calculating interior edge kappa correctly");
      
  }
  catch(const Dark_Arts_Exception& da)
  {
    val = -1;
    da.testing_message();
  }
  
  return val;
}
