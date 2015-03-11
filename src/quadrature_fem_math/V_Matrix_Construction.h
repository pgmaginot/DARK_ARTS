#ifndef V_Matrix_Construction_h
#define V_Matrix_Construction_h

#include "Inputs_Allowed.h"
#include "Fem_Quadrature.h"
#include "Materials.h"
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
  V_Matrix_Construction(const Fem_Quadrature& fem_quadrature, Materials& materials);
  virtual ~V_Matrix_Construction(){}
  
  
  
  /**
    These are the only callable material properties reaction matrix constructors
    \fn construct_r_cv , \fn construct_r_sigma_a , \fn construct_r_sigma_s
       
    The matrices returned will already include the Jacobian of the transformation
  */
   
  /// construct most of \f$ \mathbf{M} \f$, except for the \f$ \frac{\Delta x}{2} \f$ constant multiplication
  virtual void construct_dimensionless_mass_matrix(Eigen::MatrixXd& mass_mat) = 0;
    
  void construct_r_cv(Eigen::MatrixXd& r_cv);
  
  void construct_r_sigma_a(Eigen::MatrixXd& r_sig_a, const int grp);
  
  void construct_r_sigma_s( std::vector<Eigen::MatrixXd>& r_sig_s, const int grp, const int l_mom);
  void construct_r_sigma_s( Eigen::MatrixXd& r_sig_s, const int grp, const int l_mom);
  /// Construct driving source moments
  
  void construct_temperature_source_moments(Eigen::VectorXd& s_t, const double time);
  
  void construct_radiation_source_moments(Eigen::VectorXd& s_i, const double time, const int dir, const int grp);
 
  void construct_pos_gradient_matrix(Eigen::MatrixXd& l_pos);
  void construct_neg_gradient_matrix(Eigen::MatrixXd& l_neg);
  void construct_pos_upwind_vector(Eigen::VectorXd& f_pos);
  void construct_neg_upwind_vector(Eigen::VectorXd& f_neg);
  
protected:
  
  /** construct a reaction matrix, called only from any of the following:
      \fn construct_r_cv , 
      \fn construct_r_sigma_a,
      \fn construct_r_sigma_s
  */
  virtual void construct_reaction_matrix(Eigen::MatrixXd& rx_mat, std::vector<double>& xs) = 0;
  
  /// quadrature integration and scaling by Jacobian of driving source moments
  void construct_source_moments(Eigen::VectorXd& source_mom, std::vector<double>& source_evals);
  
  /// Calculate gradient quantities (without mu multiplied through)

  
  /// pointer to the Materials class object where all properties of given materials reside
  Materials& m_materials;  
  
  /**
    Store the evaluated basis functions and quadrature rules from
    Fem_Quadrature objects
  */
  /// number of DFEM basis functions (1 basis function per interpoaltion/basis point
  const int m_n_basis_pts ;

  /// DFEM matrix quadrature
  const int m_n_quad_pts ;

  
  /// source moment quadrature (may be different / more exact than matrix quadrature
  const int m_n_source_quad_pts;
  std::vector<double> m_source_weights;
  std::vector<double> m_dfem_at_source_quad;
  
  /// a vector of driving source evaluations at the source moment quadrature points
  std::vector<double> m_source_evals;
  
  /// a vector of material properties evalauted at the dfem integration points
  std::vector<double> m_xs_evals;
  
  /// quadrature weights for DFEM matrix formation
  std::vector<double> m_integration_weights;
  
  /// DFEM basis functions evalauted at quadrature points
  std::vector<double> m_basis_at_quad;
  
  /// derivatives of 
  std::vector<double> m_basis_deriv_at_quad;
  
  std::vector<double> m_basis_at_left_edge;
  std::vector<double> m_basis_at_right_edge;

};

#endif