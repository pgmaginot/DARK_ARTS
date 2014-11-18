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
  virtual ~Angular_Quadrature(){}
    
  int get_number_of_dir(void) const;
  int get_number_of_groups(void) const;
  int get_number_of_leg_moments(void) const;
  double get_leg_poly(const int dir, const int mom) const;
  double get_leg_moment_coeff_build(const int dir, const int mom) const;
  double get_mu(const int dir) const;
  double get_w(const int dir) const;
  double get_sum_w(void) const;
  
  double get_group_low_bound(const int grp) const;
  double get_group_upper_bound(const int grp) const;
  
  bool has_left_reflection(void) const;  
  
  double most_glance_mu(void) const;
  double most_normal_mu(void) const;  
protected:  
  /// number of directions
  const int m_n_dir;
  /// number of groups
  const int m_n_groups;  
  /// number of angle integrated legendre moments to keep
  const int m_n_legendre_moments;
  
  const bool m_left_reflecting_boundary;
  
  double m_sum_w;
  
  double m_mu_most_glancing;
  double m_mu_most_normal;

  /// discrete ordinates values
  std::vector<double> m_mu;
  
  /// quadrature weights
  std::vector<double> m_w;
  

  
  /// legendre polynomials evaluated at discrete ordinates
  std::vector<double> m_legendre_poly;  
  std::vector<double> m_legendre_poly_alone;

  /// frequency group bounds
  std::vector<double> m_grp_e_min;
  std::vector<double> m_grp_e_max;
};

#endif