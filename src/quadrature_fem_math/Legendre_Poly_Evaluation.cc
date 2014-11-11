/** @file   Legendre_Poly_Evaluation.cc
  *   @author  pmaginot
  *   @brief Implement the Legendre_Poly_Evaluation class that gives access to the Legendre polynomial functions
 */
         
#include "Legendre_Poly_Evaluation.h"
// ##########################################################
// Public functions 
// ##########################################################

/**
  \fn get_evaluated_legendre_polynomials - evaluate Legendre polynomails at point x.  Results stored in a vector, \param evals
   with all the Legendre evaluations at a given point stored consecutively .   \param start is the offset in the vector
   Should be noted that get_legendre_polynomials is at the mercy of the calling function to avoid overwriting / messing up data
   \param n is the degree of the highest legendre polynomial considered
*/
void Legendre_Poly_Evaluation::get_evaluated_legendre_polynomials(
  const double x, const int n, const int start, std::vector<double>& evals)
{
  try{
    if(x < -1.)
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Trying to evalaute Legendre polynomial for x<-1");
      
    if(x > 1.)
      throw Dark_Arts_Exception(SUPPORT_OBJECT , "Trying to evaluate Legendre polynomial for x > 1");
      
    if(n < 0)
    {
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Trying to calculate negative legendre moment" );
    }  
    else if(n==0)
    {
      evals[start] = 1.;
    }
    else if(n==1)
    {
      evals[start] = 1.;
      evals[start+1] = x;
    }
    else
    {
      evals[start] = 1.;
      evals[start+1] = x;
      double dbl_n = 1.;
      for(int k=2; k < (n+1); k++)
      {
        evals[start+k] = (2.*dbl_n + 1.)*x*evals[start + k-1] - dbl_n*evals[ start + k-2];
        evals[start+k] /= dbl_n + 1.;
        dbl_n += 1.;
      }
    }
  }
  catch(const Dark_Arts_Exception& da_exception)
  {
    da_exception.message() ;
  }
  
  return;
}
