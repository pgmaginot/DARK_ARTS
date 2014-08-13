#ifndef Fem_Quadrature_h
#define Fem_Quadrature_h

#include "Inputs_Allowed.h"
#include "Quadrule_New.h"
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
  Fem_Quadrature(const Input_Reader& input_reader, const Quadrule_New& quad_fun);
  ~Fem_Quadrature(){}
  
  int get_number_of_integration_points(void) const;
  
  int get_number_of_interpolation_points(void) const ;
  
  int get_number_of_xs_point(void) const;
  
  void get_xs_eval_points(std::vector<double>& xs_eval_pts) const;
  void get_dfem_at_xs_eval_points(std::vector<double>& dfem_at_xs_pts) const;
  void get_dfem_integration_points(std::vector<double>& dfem_integration_pts) const;
  void get_dfem_at_edges(std::vector<double>& dfem_at_left_edge,std::vector<double>& dfem_at_right_edge) const;
  
private:

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
  
  /// number of quadrature points that will be used to evaluate DFEM matrices/integrals (L,M,R matrices)
  int m_n_integration_points = -1;
  /// number of DFEM interpolation points.  Better equal trial space degree + 1
  int m_n_interpolation_points = -1;
  /// number of points to evaluate opacity at within each spatial cell
  int m_n_xs_evaluation_points = -1;
  /** number of opacity interpolation points: 
    relevant only when OPACITY_SPATIAL_TREATMENT = INTERPOLATING
  */
  int m_n_xs_interpolation_points = -1;
  
  /// Lagrange basis function interpolation points and associated weights
  std::vector<double> m_dfem_interpolation_points;
  std::vector<double> m_dfem_interpolation_weights;
  
  /// quadrature that reaction and gradient matrices will be evaluated with
  std::vector<double> m_integration_points;
  std::vector<double> m_integration_weights;
  
  /** local quadrature points where opacity will be evaluated in each cell
    if SLXS- equal to DFEM integration points
    if MOMENT_PRESERVING- quadrature points that willl be used to form moments
    if INTERPOLATING- these points become the cross section interpolation points */
  std::vector<double> m_xs_eval_points;
  std::vector<double> m_xs_eval_weights;
  
  /// used only if OPACITY_TREATMENT = INTERPOLATING
  std::vector<double> m_xs_poly_at_integration_points;
  
  /// m_n_integration_points times m_n_interpolation_points length vector
  std::vector<double> m_basis_at_integration_points;
  
  /// m_n_integration_points times m_n_interpolation_points length vector
  std::vector<double> m_d_basis_d_s_at_integration_points;
  
  /** Vector of DFEM basis functions at cross section evaluation points,
    * m_n_xs_evaluation_points times m_n_interpolation points length,
    * used to get temperature at the cross section evaluation points
  */
  std::vector<double> m_basis_at_xs_points;
  
  std::vector<double> m_dfem_at_left_edge;
  std::vector<double> m_dfem_at_right_edge;
  
  std::vector<double> m_d_dfem_d_s_at_left_edge;
  std::vector<double> m_d_dfem_d_s_at_right_edge;
};

#endif