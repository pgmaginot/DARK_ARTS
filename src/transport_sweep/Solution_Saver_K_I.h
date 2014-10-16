#ifndef Solution_Saver_K_I_h
#define Solution_Saver_K_I_h

#include "V_Solution_Saver.h"
#include <memory>
#include "V_Sweep_Matrix_Creator.h"
#include "K_Intensity.h"

/** @file   Solution_Saver_K_I.h
  *   @author pmaginot
  *   @brief calculate K_I of the local cell and save.  Do not overrwrite flux angular moment
 */

class Solution_Saver_K_I: public V_Solution_Saver 
{
public:
  Solution_Saver_K_I(const Fem_Quadrature& fem_quadrature,std::shared_ptr<V_Sweep_Matrix_Creator> matrix_creator_ptr,
  const Angular_Quadrature& angular_quadrature, K_Intensity& k_i_ref, const double c);
    
  virtual ~Solution_Saver_K_I(){}
  
  void save_local_solution(Intensity_Moment_Data& phi_new, 
  const Eigen::VectorXd& local_intensity, 
  Psi_In& psi_in,
  const int cell, 
  const int grp, 
  const int dir) override;
protected:
  
  std::shared_ptr<V_Sweep_Matrix_Creator> m_sweep_matrix_ptr;
  K_Intensity& m_k_i_ref;
  
  Eigen::VectorXd m_local_ki;
  Eigen::VectorXd m_scratch_vec;
  Eigen::MatrixXd m_scratch_mat;
  
  const int m_zero;
  const double m_c;
  int m_stage;
  
};

#endif