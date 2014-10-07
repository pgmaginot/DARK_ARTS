#ifndef Solution_Saver_K_I_h
#define Solution_Saver_K_I_h

#include "V_Solution_Saver.h"

/** @file   Solution_Saver_K_I.h
  *   @author pmaginot
  *   @brief calculate K_I of the local cell and save.  Do not overrwrite flux angular moment
 */

class Solution_Saver_K_I: public V_Solution_Saver 
{
public:
  Solution_Saver_K_I(const Fem_Quadrature& fem_quadrature);
    
  virtual ~Solution_Saver_K_I(){}

  void save_local_solution(Eigen::VectorXd& local_intensity, const int cell, const int grp, const int dir) override;
protected:
  
};

#endif