#ifndef Angular_Quadrature_h
#define Angular_Quadrature_h

#include "Inputs_Allowed.h"
#include "Input_Reader.h"
#include "Quadrule_New.h"
#include "Legendre_Poly_Evaluation.h"

#include <vector>
#include <stdlib.h>
#include <iostream>
#include <sstream>

/** @file   Angular_Quadrature.hm_
  *   @author pmaginot
  *   @brief Angular_Quadrature class, contain group info, disrete ordinates info, and legendre evaluations of angular quadrature
 */

class Angular_Quadrature
{
public:
  /// Only able to initialize if given an Input_Reader object
  /// constructor defined in Fem_Quadrature.cc
  Angular_Quadrature(const Input_Reader& input_reader, const Quadrule_New& quad_fun );
  ~Angular_Quadrature(){}
    
  int get_number_of_dir(void) const;
  int get_number_of_groups(void) const;
  int get_number_of_leg_moments(void) const;
  double get_leg_poly(const int dir, const int mom) const;
  double get_leg_moment_coeff(const int dir, const int mom) const;
  double get_mu(const int dir) const;
  double get_w(const int dir) const;
  double get_sum_w(void) const;
  
  double calculate_boundary_conditions(const int dir, const int grp, const double time) const;
  
protected:  
  /// number of directions
  const int m_n_dir;
  /// number of groups
  const int m_n_groups;  
  /// number of angle integrated legendre moments to keep
  const int m_n_legendre_moments;
  

  /// discrete ordinates values
  std::vector<double> m_mu;
  
  /// quadrature weights
  std::vector<double> m_w;
  
  double m_sum_w;
  
  /// legendre polynomials evaluated at discrete ordinates
  std::vector<double> m_legendre_poly;  

};

#endif