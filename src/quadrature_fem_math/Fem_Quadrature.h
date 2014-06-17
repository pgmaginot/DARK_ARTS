#ifndef Fem_Quadrature_h
#define Fem_Quadrature_h

#include "Inputs_Allowed.h"
#include "Input_Reader.h"

#include <vector>
#include <stdlib.h>
#include <iostream>
#include <sstream>

/** @file   Fem_Quadrature.h
  *   @author pmaginot
  *   @brief Implement the Fem_Quadrature class that owns all things quadrature related
 */

class Fem_Quadrature
{
public:
  /// Only able to initialize if given an Input_Reader object
  /// constructor defined in Fem_Quadrature.cc
  Fem_Quadrature(Input_Reader&  input_reader);
  ~Fem_Quadrature(){}
  
  int get_number_of_integration_points(void);
  
protected:

/* ****************************************************
*
*     Protected Functions
*
  **************************************************** */
  void evaluate_lagrange_func(const std::vector<double>& interp_points, 
    const std::vector<double>& eval_points, std::vector<double>& func_evals);

  void evaluate_lagrange_func_derivatives(const std::vector<double>& interp_points, 
    const std::vector<double>& eval_points, std::vector<double>& deriv_evals);

/* ****************************************************
*
*     Protected Variables
*
  **************************************************** */
  const int m_xs_extra_points = 10;
  
  int m_n_integration_points = -1;
  int m_n_interpolation_points = -1;
  int m_n_xs_interpolation_points = -1;
  
  /// Lagrange basis function interpolation points and associated weights
  std::vector<double> m_dfem_interpolation_points;
  std::vector<double> m_dfem_interpolation_weights;
  
  /// quadrature that reaction and gradient matrices will be evaluated with
  std::vector<double> m_integration_points;
  std::vector<double> m_integration_weights;
  
  /// quadrature points that will be used in forming xs at integration points
  std::vector<double> m_xs_eval_points;
  std::vector<double> m_xs_eval_weights;
  
  /// m_n_integration_points times m_n_interpolation_points length vector
  std::vector<double> m_basis_at_integration_points;
  
  /// m_n_integration_points times m_n_interpolation_points length vector
  std::vector<double> m_d_basis_d_s_at_integration_points;
  
};

#endif