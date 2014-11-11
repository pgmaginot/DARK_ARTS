#include <iostream>
#include <vector>
#include <math.h>
#include <string>

#include "Legendre_Poly_Evaluation.h"
#include "Dark_Arts_Exception.h"

class Test
{
public:
  Test(){}
  ~Test(){}
  
  void get_expected_legendre_polynomials(const double x , std::vector<double>& expected)
  {
    expected[0] = 1.;
    expected[1] = x;
    expected[2] = (3.*x*x - 1.)/2.;
    expected[3] = (5.*pow(x,3.) - 3.*x)/2.;
    expected[4] = (35.*pow(x,4.) - 30.*x*x + 3.)/8.;
    return;
  }
};

int main(int argc, char** argv)
{
  int val = 0;
  
  /** The legendre polynomials are:
  
    \f{eqnarray}
    {
      P_0(x) &=& 1 \\
      P_1(x) &=& x \\
      P_2(x) &=& \frac{1}{2} \left( 3x^2 - 1 \right) \\
      P_3(x) &=& \frac{1}{2} \left( 5x^3 - 3x \right) \\
      P_4(x) &=& \frac{1}{8} \left( 35x^4 - 30x^2 + 3 \right)
    \f}    
  */
  
  Legendre_Poly_Evaluation leg_poly;
  Test test;
  
  const int n_leg_moment = 4;
  std::vector<double> calculated_leg_evals(n_leg_moment+1,0.);
  
  std::vector<double> expected_leg_evals(n_leg_moment+1, 0.);
  
  const int offset = 0;
  
  /// try at a variety of x_eval points, first being -1
  double x_eval = -1.;
  leg_poly.get_evaluated_legendre_polynomials(x_eval,n_leg_moment , offset , calculated_leg_evals);
  test.get_expected_legendre_polynomials(x_eval, expected_leg_evals);
  try{
    std::cout << "x_eval= " << x_eval << std::endl;
    for(int i=0; i<= n_leg_moment ; i++)
    {
      std::cout << "For P_"<< i << " Expected poly of: " << expected_leg_evals[i] << " Got: " << calculated_leg_evals[i] << std::endl ;
      if( fabs(expected_leg_evals[i] - calculated_leg_evals[i] ) > 1.E-4)
      {
        std::stringstream err;
        err << "For x= " << x_eval << " difference in Legendre polynomial detected\n";
        throw Dark_Arts_Exception(SUPPORT_OBJECT , err.str() );
      }
    
    }
  }
  catch(const Dark_Arts_Exception& da_exception)
  {
    da_exception.testing_message();
    val = -1;
  } 
  
  x_eval = -0.2;
  leg_poly.get_evaluated_legendre_polynomials(x_eval,n_leg_moment , offset , calculated_leg_evals);
  test.get_expected_legendre_polynomials(x_eval, expected_leg_evals);
  try{
    std::cout << "x_eval= " << x_eval << std::endl;
    for(int i=0; i<= n_leg_moment ; i++)
    {
      std::cout << "For P_"<< i << " Expected poly of: " << expected_leg_evals[i] << " Got: " << calculated_leg_evals[i] << std::endl ;
      if( fabs(expected_leg_evals[i] - calculated_leg_evals[i] ) > 1.E-4)
      {
        std::stringstream err;
        err << "For x= " << x_eval << " difference in Legendre polynomial detected\n";
        throw Dark_Arts_Exception(SUPPORT_OBJECT , err.str() );
      }
    
    }
  }
  catch(const Dark_Arts_Exception& da_exception)
  {
    da_exception.testing_message();
    val = -1;
  }
  
  x_eval = 0.;
  leg_poly.get_evaluated_legendre_polynomials(x_eval,n_leg_moment , offset , calculated_leg_evals);
  test.get_expected_legendre_polynomials(x_eval, expected_leg_evals);
  try{
    std::cout << "x_eval= " << x_eval << std::endl;
    for(int i=0; i<= n_leg_moment ; i++)
    {
      std::cout << "For P_"<< i << " Expected poly of: " << expected_leg_evals[i] << " Got: " << calculated_leg_evals[i] << std::endl ;
      if( fabs(expected_leg_evals[i] - calculated_leg_evals[i] ) > 1.E-4)
      {
        std::stringstream err;
        err << "For x= " << x_eval << " difference in Legendre polynomial detected\n";
        throw Dark_Arts_Exception(SUPPORT_OBJECT , err.str() );
      }
    
    }
  }
  catch(const Dark_Arts_Exception& da_exception)
  {
    da_exception.testing_message();
    val = -1;
  }
  
  x_eval = 0.4;
  leg_poly.get_evaluated_legendre_polynomials(x_eval,n_leg_moment , offset , calculated_leg_evals);
  test.get_expected_legendre_polynomials(x_eval, expected_leg_evals);
  try{
    std::cout << "x_eval= " << x_eval << std::endl;
    for(int i=0; i<= n_leg_moment ; i++)
    {
      std::cout << "For P_"<< i << " Expected poly of: " << expected_leg_evals[i] << " Got: " << calculated_leg_evals[i] << std::endl ;
      if( fabs(expected_leg_evals[i] - calculated_leg_evals[i] ) > 1.E-4)
      {
        std::stringstream err;
        err << "For x= " << x_eval << " difference in Legendre polynomial detected\n";
        throw Dark_Arts_Exception(SUPPORT_OBJECT , err.str() );
      }
    
    }
  }
  catch(const Dark_Arts_Exception& da_exception)
  {
    da_exception.testing_message();
    val = -1;
  }
  
  x_eval = 1.;
  leg_poly.get_evaluated_legendre_polynomials(x_eval,n_leg_moment , offset , calculated_leg_evals);
  test.get_expected_legendre_polynomials(x_eval, expected_leg_evals);
  try{
    std::cout << "x_eval= " << x_eval << std::endl;
    for(int i=0; i<= n_leg_moment ; i++)
    {
      std::cout << "For P_"<< i << " Expected poly of: " << expected_leg_evals[i] << " Got: " << calculated_leg_evals[i] << std::endl ;
      if( fabs(expected_leg_evals[i] - calculated_leg_evals[i] ) > 1.E-4)
      {
        std::stringstream err;
        err << "For x= " << x_eval << " difference in Legendre polynomial detected\n";
        throw Dark_Arts_Exception(SUPPORT_OBJECT , err.str() );
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


