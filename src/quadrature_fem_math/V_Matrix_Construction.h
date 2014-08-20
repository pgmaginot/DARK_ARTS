#ifndef V_Matrix_Construction_h
#define V_Matrix_Construction_h

#include "Inputs_Allowed.h"
#include "Fem_Quadrature.h"
#include "Eigen/Dense"

#include <vector>
#include <stdlib.h>

/** @file   V_Matrix_Construction.h
  *   @author pmaginot
  *   @brief Provide a base class that constructs mass, reaction, and gradient matrices, as well as upwind vectors
  *     Concrete cases for  SELF_LUMPING , TRAD_LUMPING , and EXACT integration techniques
 */

class V_Matrix_Construction
{
public:
  /// Only able to initialize if given an Input_Reader object
  /// constructor defined in Fem_Quadrature.cc
  V_Matrix_Construction(const Fem_Quadrature& fem_quadrature);
  virtual ~V_Matrix_Construction(){}
  
  virtual void construct_mass_matrix(Eigen::MatrixXd& mass_mat) = 0;
  
  virtual void construct_reaction_matrix(Eigen::MatrixXd& rx_mat, std::vector<double>& xs) = 0;
  
  void construct_pos_gradient_matrix(Eigen::MatrixXd& l_mat);
  
  void construct_neg_gradient_matrix(Eigen::MatrixXd& l_mat);
  
  void construct_left_upwind_vector(Eigen::VectorXd& f_mu_pos);
  
  void construction_right_upwind_vector(Eigen::VectorXd& f_mu_neg);
  
protected:

  

/* ****************************************************
*
*     Protected Functions
*
  **************************************************** */


/* ****************************************************
*
*     Protected Variables
*
  **************************************************** */
  
  /**
    Store the evaluated basis functions and quadrature rules from
    Fem_Quadrature objects
  */
  
  int m_n_quad_pts = -1;
  int m_n_basis_pts = -1;
  
  std::vector<double> m_integration_weights;
  
  std::vector<double> m_basis_at_quad;
  
private:
  /// only used in the gradient matrix and upwind contributions, don't need to copy to derived matrix construction types
  std::vector<double> m_basis_deriv_at_quad;
  
  std::vector<double> m_basis_at_left_edge;
  std::vector<double> m_basis_at_right_edge;
};

#endif