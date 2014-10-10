#ifndef V_Solution_Saver_h
#define V_Solution_Saver_h

#include "Fem_Quadrature.h"
#include "Angular_Quadrature.h"
#include "Intensity_Moment_Data.h"
#include "Psi_In.h"
#include "Eigen/Dense"

/** @file   V_Solution_Saver.h
  *   @author pmaginot
  *   @brief provide a function interface to save local intensity solution during a transport sweep.  Two options:
  * 1) save angle integrated moments only, Solution_Saver_Moments objects
  * 2) calculate and save k_I, Solution_Saver_K_I objects  
 */

class V_Solution_Saver
{
public:
  V_Solution_Saver(const Fem_Quadrature& fem_quadrature, const Angular_Quadrature& angular_quadrature);
    
  virtual ~V_Solution_Saver(){}

  virtual void save_local_solution(Intensity_Moment_Data& phi_new, 
  const Eigen::VectorXd& local_intensity, 
  Psi_In& psi_in,
  const int cell, 
  const int grp, 
  const int dir) = 0;
protected:
  const int m_np;
  const int m_n_dir_div_2;
  const int m_n_l_mom;
  
  double m_outflow;
  
  const Angular_Quadrature& m_quad_ref; 
  
  std::vector<double> m_dfem_at_left_bound;
  std::vector<double> m_dfem_at_right_bound;
  
  double calculate_outflow(const int dir, const Eigen::VectorXd& local_intensity);
};

#endif