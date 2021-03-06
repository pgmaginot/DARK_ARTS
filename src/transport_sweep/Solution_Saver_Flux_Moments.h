#ifndef Solution_Saver_Flux_Moments_h
#define Solution_Saver_Flux_Moments_h

#include "V_Solution_Saver.h"


/** @file   Solution_Saver_Flux_Moments.h
  *   @author pmaginot
  *   @brief save the angle integrated integrated flux moments of the local solution obtained during a flux sweep
 */

class Solution_Saver_Flux_Moments: public V_Solution_Saver 
{
public:
  Solution_Saver_Flux_Moments(const Fem_Quadrature& fem_quadrature, const Angular_Quadrature& anggular_quadrature);
    
  virtual ~Solution_Saver_Flux_Moments(){}

  void save_local_solution(Intensity_Moment_Data& phi_new, 
  const Eigen::VectorXd& local_intensity, 
  Psi_In& psi_in,
  const int cell, 
  const int grp, 
  const int dir) override;
protected:

  Eigen::VectorXd m_loc_val;
};

#endif