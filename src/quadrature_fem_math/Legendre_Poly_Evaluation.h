/** 
  @file   Legendre_Poly_Evaluation.h  
  @author pmaginot
  @brief implement evaluations of the legendre polynomials via recursion
  
  \f[
      \left( \begin{array}{c} n \\ k \end{array} \right) = 
        \prod_{i=1}^{k}{ \frac{n+1-i}{i} }
  \f]
*/
#ifndef Legendre_Poly_Evaluation_h
#define Legendre_Poly_Evaluation_h

#include <vector>
#include <stdlib.h>
#include <iostream>


class Legendre_Poly_Evaluation
{
/** Allow the radiative transfer code to call these quadratures, no others
 others may require the modification and inclusion of other routines from the orignial quadrule.hpp and quadrule.cpp */
public:
  Legendre_Poly_Evaluation(){}
  
  ~Legendre_Poly_Evaluation(){} 
   
   void get_evaluated_legendre_polynomials(const double x, const int n, std::vector<double>& evals);
   
   void get_legendre_polynomials_explicitly(const double x, const int n, std::vector<double>& evals);
private:
  double evaluate_binomial_coefficient(const int n, const int k);
  
};

#endif
