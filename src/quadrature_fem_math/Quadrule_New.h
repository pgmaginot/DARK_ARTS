/** 
  @file   Quadrule_New.h  
  @author The internet and pmaginot
  @brief implementations of various quadrature rules and stepping routines
    This file and Quadrule_New.cc were originally taken from QUADRULE:
      http://people.sc.fsu.edu/~jburkardt/cpp_src/quadrule/quadrule.html
    
    They have been modified to become the Quadrule_New class, so that the functions are not
    global in the global namespace 
*/
#ifndef Quadrule_New_h
#define Quadrule_New_h

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
#include <vector>
#include <stdlib.h>


class Quadrule_New
{
/** Allow the radiative transfer code to call these quadratures, no others
 others may require the modification and inclusion of other routines from the orignial quadrule.hpp and quadrule.cpp */
public:
  Quadrule_New(){}
  
  ~Quadrule_New(){}
  
  void legendre_dr_compute ( const int order, std::vector<double>& xtab, std::vector<double> weight );
  void legendre_ek_compute ( const int n, std::vector<double>& x, std::vector<double>& w );
  void legendre_set ( const int order, std::vector<double>& xtab, std::vector<double>& weight );

  void lobatto_compute ( const int n, std::vector<double>& x, std::vector<double>& w );
  void lobatto_set ( const int order, std::vector<double>& xtab, std::vector<double>& weight );

  void ncc_compute ( const int order, std::vector<double>& xtab, std::vector<double>& weight );  
  void ncc_set ( const int order, std::vector<double>& xtab, std::vector<double>& weight );
  
private:
  
  void imtqlx (  const int n, std::vector<double>& d, std::vector<double>& e, 
    std::vector<double>& z );

  void nc_compute ( const int order, const double a, const double b,
    const std::vector<double>& xtab, std::vector<double>& weight );

  void ncc_compute_points ( const int n, std::vector<double>& x );
  
  void ncc_compute_weights ( const int n, std::vector<double>& w );

  double r8_epsilon (void);
  
  double r8_max ( const double x, const double y );

  double r8_sign ( const double x );

  void r8vec_reverse ( const int n, std::vector<double>& x );
};

#endif
