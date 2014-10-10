#include "V_Solution_Saver.h"

V_Solution_Saver::V_Solution_Saver(const Fem_Quadrature& fem_quadrature, const Angular_Quadrature& angular_quadrature)
:
m_np{fem_quadrature.get_number_of_interpolation_points()},
m_n_dir_div_2{ angular_quadrature.get_number_of_dir() },
m_n_l_mom{angular_quadrature.get_number_of_leg_moments() },
m_outflow{0.},
m_stage{0},
m_quad_ref(angular_quadrature)
{
 fem_quadrature.get_dfem_at_edges(m_dfem_at_left_bound,m_dfem_at_right_bound);
}

double V_Solution_Saver::calculate_outflow(const int dir, const Eigen::VectorXd& local_intensity)
{
  m_outflow = 0.;
  if(dir < m_n_dir_div_2)
  {
    /// mu < 0 outflow is the left edge
    for(int i=0;i<m_np;i++)
    {
      m_outflow += m_dfem_at_left_bound[i]*local_intensity(i);
    }
  }
  {
    /// mu > 0 outlfow is the right edge
    for(int i=0;i<m_np;i++)
    {
      m_outflow += m_dfem_at_right_bound[i]*local_intensity(i);
    }
  }
  return m_outflow;
}

void  V_Solution_Saver::set_stage(const int stage)
{
  m_stage = stage;
  return;
}
