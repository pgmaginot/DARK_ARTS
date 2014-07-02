/** @file   Legendre_Poly_Evaluation.cc
  *   @author  pmaginot
  *   @brief Implement the Legendre_Poly_Evaluation class that gives access to the Legendre polynomial functions
 */
         
#include "Legendre_Poly_Evaluation.h"
// ##########################################################
// Public functions 
// ##########################################################

void Legendre_Poly_Evaluation::get_evaluated_legendre_polynomials(
  const double x, const int n, std::vector<double>& evals)
{
  if(n < 0)
  {
    std::cerr << "Trying to calculate negative legendre moment\n" ;
    exit(EXIT_FAILURE);
  }  
  else if(n==0)
  {
    evals[0] = 1.;
  }
  else if(n==1)
  {
    evals[0] = 1.;
    evals[1] = x;
  }
  else
  {
    double dbl_n = 2.;
    for(int k=2; k < (n+1); k++)
    {
      evals[k] = (2.*dbl_n + 1.)*x*evals[k-1] - dbl_n*evals[k-2];
      evals[k] /= dbl_n + 1.;
      dbl_n += 1.;
    }
  }
  
  return;
}

void Legendre_Poly_Evaluation::get_legendre_polynomials_explicitly(
  const double x, const int n, std::vector<double>& evals)
{
  return;
}

// ##########################################################
// Private functions 
// ##########################################################

double Legendre_Poly_Evaluation::evaluate_binomial_coefficient(const int n, const int k)
{
  const double dbl_n = double(n);
  double i_dbl= 1.;
  
  double ans = 1.;
  
  for(int i=0;i<k ;i++)
  {
    ans *= (dbl_n + 1. - i_dbl)/i_dbl;
    i_dbl += 1.;
  }
  
  return ans;
}