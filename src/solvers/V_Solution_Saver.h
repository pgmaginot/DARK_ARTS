#ifndef V_Solution_Saver_h
#define V_Solution_Saver_h

#include "Fem_Quadrature.h"
#include "Fem_Quadrature.h"
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
  V_Solution_Saver(const Fem_Quadrature& fem_quadrature);
    
  virtual ~V_Solution_Saver(){}

  virtual void save_local_solution(Eigen::VectorXd& local_intensity, const int cell, const int grp, const int dir) = 0;
protected:
  const int m_n_dfem_pts;
  
};

#endif