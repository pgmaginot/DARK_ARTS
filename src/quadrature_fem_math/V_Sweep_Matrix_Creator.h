#ifndef V_Sweep_Matrix_Creator_h
#define V_Sweep_Matrix_Creator_h

#include "Inputs_Allowed.h"
#include "Fem_Quadrature.h"
#include "Eigen/Dense"
#include "V_Matrix_Construction.h"
#include "Matrix_Construction_Exact.h"
#include "Matrix_Construction_Self_Lumping.h"
#include "Matrix_Construction_Trad_Lumping.h"

#include "Cell_Data.h"
#include "Materials.h"
#include "Temperature_Data.h"
#include "K_Intensity.h"
#include "K_Temperature.h"

#include <vector>
#include <memory>

/** @file   V_Matrix_Construction.h
  *   @author pmaginot
  *   @brief Provide a base class that constructs mass, reaction, and gradient matrices, as well as upwind vectors
  *     Concrete cases for  SELF_LUMPING , TRAD_LUMPING , and EXACT integration techniques
 */

class V_Sweep_Matrix_Creator
{
public:
  /// Only able to initialize if given an Input_Reader object
  /// constructor defined in Fem_Quadrature.cc
  V_Sweep_Matrix_Creator(const Fem_Quadrature& fem_quadrature);
  virtual ~V_Sweep_Matrix_Creator(){}
  
  virtual void construct_r_sig_t(const int cell, const int grp, Eigen::MatrixXd& r_sig_t) = 0;
  
  virtual void construct_r_sig_s(const int cell, const int grp, const int l_mom, Eigen::MatrixXd& r_sig_s) = 0;
  
  virtual void construct_s_i(const int cell,const int grp, const int l_mom, Eigen::VectorXd& s_i) = 0;
  
  void construct_l_matrix(const double mu, Eigen::MatrixXd& l_matrix);
  
  void construct_f_vector(const double mu, Eigen::VectorXd& f_vector);
  
  /// data that changes every time we update radiation intensity
  // void set_thermal_iteration_data()
  
  // /** data that changes once per stage (per time step)
    // stage number, a_vector, time
  // */
  // void set_stage_data()
  
  // /** data that only changes once per time step
    // dt
  // */
  // void set_timestep_data()
  
  
private:


  const MATRIX_INTEGRATION m_matrix_type; 
  
  const int m_np;
  
  Eigen::MatrixXd m_no_mu_pos_l_matrix;
  Eigen::MatrixXd m_no_mu_neg_l_matrix;
  Eigen::VectorXd m_no_mu_pos_f_vector;
  Eigen::VectorXd m_no_mu_neg_f_vector;
  
  /// object that can build reaction matrices
  std::shared_ptr<V_Matrix_Construction> m_mtrx_builder;
  
  const Temperature_Data* m_t_star;
  const Temperature_Data* m_t_old;
  
  const K_Temperature* m_k_t;
  const K_Intensity* m_k_i;
  
  const Cell_Data* m_cell_data;
  
  /// can be const because when we call Materials::m_
  Materials* m_material;
};

#endif