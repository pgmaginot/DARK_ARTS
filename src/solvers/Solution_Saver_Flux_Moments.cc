#include "Solution_Saver_Flux_Moments.h"

Solution_Saver_Flux_Moments::Solution_Saver_Flux_Moments(const Fem_Quadrature& fem_quadrature, const Angular_Quadrature& angular_quadrature)
:
V_Solution_Saver(fem_quadrature,angular_quadrature),
m_loc_val{Eigen::VectorXd::Zero(m_np)}
{
 
}


void Solution_Saver_Flux_Moments::save_local_solution(Intensity_Moment_Data& phi_new, 
  const Eigen::VectorXd& local_intensity, 
  Psi_In& psi_in,
  const int cell, 
  const int grp, 
  const int dir)
{
  for(int l_mom=0;l_mom<m_n_l_mom;l_mom++)
  {
    m_loc_val = m_quad_ref.get_w(dir)*m_quad_ref.get_leg_poly(dir,l_mom)*local_intensity;
    phi_new.add_contribution(cell,grp,l_mom,m_loc_val  );
  }
  
  psi_in(grp,dir) = calculate_outflow(dir, local_intensity);
  
  return;
}


