#include "Transport_Sweep.h"

Transport_Sweep::Transport_Sweep(const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* material, 
  const Angular_Quadrature& angular_quadrature, const Time_Stepper& time_stepper)
  :
  m_n_cells{ cell_data->get_total_number_of_cells() },
  m_n_groups{ angular_quadrature.get_number_of_groups()  },
  m_n_dir{ angular_quadrature.get_number_of_dir() },
  m_n_l_mom{ angular_quadrature.get_number_of_leg_moments() },
  m_np{ fem_quadrature.get_number_of_interpolation_points() }, 
  m_matrix_scratch{ Eigen::MatrixXd(m_np,m_np) },
  m_rhs_vec{ Eigen::VectorXd::Zero(m_np) },
  m_lhs_mat{ Eigen::MatrixXd(m_np,m_np) },
  m_local_soln{ Eigen::VectorXd::Zero(m_np)  }
{
 
}

void Transport_Sweep::sweep_mesh(const bool is_krylov, Intensity_Data& intensity_new, const Intensity_Data& intensity_old,
  Temperature_Data& t_new, const Temperature_Data& t_star, const Temperature_Data& t_n,
  const K_Temperature& k_t, const K_Intensity& k_i, const int stage, 
  const std::vector<double>& rk_a, const double time, const double dt)
{
  return;
}
