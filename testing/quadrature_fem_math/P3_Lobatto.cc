#include <iostream>
#include <vector>
#include <math.h>
#include <string>

#include "Quadrule_New.h"
#include "Dark_Arts_Exception.h"

int main(int argc, char** argv)
{
  int val = 0;
  
  const int n_p = 4;
  
  std::vector<double> expected_roots(n_p,0.);
  std::vector<double> expected_weights(n_p,0.);
  
  expected_roots[0] = -1.;
  expected_roots[1] = -0.4472135955;
  expected_roots[2] = 0.4472135955;
  expected_roots[3] = 1.;
  
  expected_weights[0] = 0.1666666667;
  expected_weights[1] = 0.8333333333;
  expected_weights[2] = 0.8333333333;
  expected_weights[3] = 0.1666666667;
  
  Quadrule_New quad_fun;
  
  std::vector<double> calculated_roots(n_p,0.);
  std::vector<double> calculated_weights(n_p,0.);
  
  double tol = 1.0E-6;
  try{
    std::cout << "lobatto_compute \n";
    for(int i=0; i < n_p; i++)
    {
      calculated_roots[i] = 0.;
      calculated_weights[i] = 0.;
    }
    quad_fun.lobatto_compute( n_p , calculated_roots , calculated_weights);
    for(int i=0; i < n_p ; i++)
    {
      std::cout << "Calculated x["<<i<<"]: " << calculated_roots[i] << " w["<<i<<"]: " << calculated_weights[i] <<std::endl;
      if( (fabs( (calculated_weights[i] - expected_weights[i] )/calculated_weights[i]) > tol) )
      {
        std::stringstream err;
        err << "Lobatto Quadrature: Too large of difference in calculated weight " << i << " actual weight should be: " << expected_weights[i] ;
        throw Dark_Arts_Exception(SUPPORT_OBJECT , err.str() );
        break;
      }
      if( (fabs( (calculated_roots[i] - expected_roots[i] )/calculated_roots[i] ) > tol) )
      {
        std::stringstream err;
        err << "Lobatto Quadrature: Too large of difference in calculated root " << i << " actual root should be: " << expected_roots[i] ;
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
    quad_fun.lobatto_set( n_p , calculated_roots , calculated_weights);
    
    for(int i=0; i < n_p ; i++)
    {
      std::cout << "Calculated x["<<i<<"]: " << calculated_roots[i] << " w["<<i<<"]: " << calculated_weights[i] <<std::endl;
      if( (fabs( (calculated_weights[i] - expected_weights[i] )/calculated_weights[i]) > tol) )
      {
        std::stringstream err;
        err << "Lobatto Quadrature: Too large of difference in calculated weight " << i << " actual weight should be: " << expected_weights[i] ;
        throw Dark_Arts_Exception(SUPPORT_OBJECT , err.str() );
        break;
      }
      if( (fabs( (calculated_roots[i] - expected_roots[i] )/calculated_roots[i]) > tol) )
      {
        std::stringstream err;
        err << "Lobatto Quadrature: Too large of difference in calculated root " << i << " actual root should be: " << expected_roots[i] ;
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
