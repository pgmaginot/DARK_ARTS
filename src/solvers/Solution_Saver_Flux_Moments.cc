#include "Solution_Saver_Flux_Moments.h"

Solution_Saver_Flux_Moments::Solution_Saver_Flux_Moments(const Fem_Quadrature& fem_quadrature)
:
V_Solution_Saver(fem_quadrature)
{
 
}


void Solution_Saver_Flux_Moments::save_local_solution(Eigen::VectorXd& local_intensity, const int cell, const int grp, const int dir)
{

  return;
}


