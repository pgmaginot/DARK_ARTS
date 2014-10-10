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
  
  /// accept/return vector of local quadrature points where xs will be evaluated
  void get_xs_eval_points(std::vector<double>& xs_eval_pts) const;
  void get_xs_eval_weights(std::vector<double>& xs_eval_wts) const;
  
  /// accept/return dfem basis functions at the xs_evaluation points (used to find temperature at xs evalaution points)
  void get_dfem_at_xs_eval_points(std::vector<double>& dfem_at_xs_pts) const;
  
  /// accept/return vector of dfem basis functions at dfem integration points
  void get_dfem_integration_points(std::vector<double>& dfem_integration_pts) const;
  
  /// accept/return vectors of dfem basis functions at the left and right edges
  void get_dfem_at_edges(std::vector<double>& dfem_at_left_edge,std::vector<double>& dfem_at_right_edge) const;
  
  /// accept/return a vector of the interpolatory cross section polynomials evaluated at the dfem integration points
  void get_xs_at_dfem_integration_points(std::vector<double>& xs_at_dfem_integration_pts) const;
  
  void get_dfem_at_integration_points(std::vector<double>& dfem_at_int_pts) const;
  
  void get_dfem_derivatives_at_integration_points(std::vector<double>& dfem_deriv_at_int_pts) const;

  MATRIX_INTEGRATION get_integration_type(void) const;
  
  void get_integration_weights(std::vector<double>& m_integration_weights) const;
  
  void get_dfem_at_source_points(std::vector<double>& dfem_at_source_quad) const;

  void get_source_weights(std::vector<double>& source_weights) const;
  
  void get_source_points(std::vector<double>& source_pts) const;
  
  int get_number_of_source_points(void) const;
  
  void get_dfem_interpolation_point_weights(std::vector<double>& dfem_weights) const;
private:

  void evaluate_lagrange_func(const std::vector<double>& interp_points, 
    const std::vector<double>& eval_points, std::vector<double>& func_evals);

  void evaluate_lagrange_func_derivatives(const std::vector<double>& interp_points, 
    const std::vector<double>& eval_points, std::vector<double>& deriv_evals);

  
  const int m_n_interpolation_points;
  
  const MATRIX_INTEGRATION m_int_method;
  
  const int m_n_source_points;
  
  const int m_xs_extra_points = 10;  
  
  
  /// number of quadrature points that will be used to evaluate DFEM matrices/integrals (L,M,R matrices)
  int m_n_integration_points = -1;
  /// number of DFEM interpolation points.  Better equal trial space degree + 1
  
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
  
  /// driving source moment evaluation quadrature
  std::vector<double> m_source_points;
  std::vector<double> m_source_weights;
  
  /// dfem basis functions at driving source moment quadrature points
  std::vector<double> m_dfem_at_source_moments;
};

#endif