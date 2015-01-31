#include <iostream>
#include <vector>
#include <math.h>
#include <string>

#include "Quadrule_New.h"
#include "Dark_Arts_Exception.h"

int main(int argc, char** argv)
{
  int val = 0;
  
  const int n_p = 3;
  
  std::vector<double> expected_roots(n_p,0.);
  std::vector<double> expected_weights(n_p,0.);
  
  expected_roots[0] = -0.7745966692;
  expected_roots[1] = 0.0000000000;
  expected_roots[2] = 0.7745966692;
  
  expected_weights[0] = 0.5555555556;
  expected_weights[1] = 0.8888888889;
  expected_weights[2] = 0.5555555556;
  
  Quadrule_New quad_fun;
  
  std::vector<double> calculated_roots(n_p,0.);
  std::vector<double> calculated_weights(n_p,0.);
  
  double tol = 1.0E-6;
  
  try{
    std::cout << "legendre_ek_compute \n";
    for(int i=0; i < n_p; i++)
    {
      calculated_roots[i] = 0.;
      calculated_weights[i] = 0.;
    }
    quad_fun.legendre_ek_compute( n_p , calculated_roots , calculated_weights);
    
    
    for(int i=0; i < n_p ; i++)
    {
      std::cout << "Calculated x["<<i<<"]: " << calculated_roots[i] << " w["<<i<<"]: " << calculated_weights[i] <<std::endl;
      if( (fabs( (calculated_weights[i] - expected_weights[i] )) > tol ) )
      {
        std::stringstream err;
        err << "Gauss Quadrature: Too large of difference in calculated weight " << i << " actual weight should be: " << expected_weights[i] ;
        throw Dark_Arts_Exception(SUPPORT_OBJECT , err.str() );
        break;
      }
      if( (fabs( (calculated_roots[i] - expected_roots[i] )) > tol) )
      {
        std::stringstream err;
        err << "Gauss Quadrature: Too large of difference in calculated root " << i << " actual root should be: " << expected_roots[i] ;
        throw Dark_Arts_Exception(SUPPORT_OBJECT , err.str() );
        break;
      }
    }
  }
  catch(const Dark_Arts_Exception& da_exception)
  {
    da_exception.testing_message();
    val = -1;
  }
  try{
    std::cout << "legendre_set\n";
    for(int i=0; i < n_p; i++)
    {
      calculated_roots[i] = 0.;
      calculated_weights[i] = 0.;
    }
    quad_fun.legendre_set( n_p , calculated_roots , calculated_weights);
    
    for(int i=0; i < n_p ; i++)
    {
      std::cout << "Calculated x["<<i<<"]: " << calculated_roots[i] << " w["<<i<<"]: " << calculated_weights[i] <<std::endl;
      if( (fabs( (calculated_weights[i] - expected_weights[i] )) > tol) )
      {
        std::stringstream err;
        err << "Gauss Quadrature: Too large of difference in calculated weight " << i << " actual weight should be: " << expected_weights[i] ;
        throw Dark_Arts_Exception(SUPPORT_OBJECT , err.str() );
        break;
      }
      if( (fabs( (calculated_roots[i] - expected_roots[i] ) ) > tol ) )
      {
        std::stringstream err;
        err << "Gauss Quadrature: Too large of difference in calculated root " << i << " actual root should be: " << expected_roots[i] ;
        throw Dark_Arts_Exception(SUPPORT_OBJECT , err.str() );
        break;
      }
    }
  }
  catch(const Dark_Arts_Exception& da_exception)
  {
    da_exception.testing_message();
    val = -1;
  }  
  
  // Return 0 if tests passed, something else if failing
  return val;
}
