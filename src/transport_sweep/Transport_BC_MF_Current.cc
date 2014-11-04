#include "Transport_BC_MF_Current.h"



Transport_BC_MF_Current::Transport_BC_MF_Current(const Angular_Quadrature& angular_quadrature,
    BC_ANGLE_DEPENDENCE incident_dependence,
    const double t_start,
    const double t_end,
    const double energy_current,
    const BC_ENERGY_DEPENDENCE e_dependence)
:
  V_Transport_BC(),
  m_sn_w{ angular_quadrature.get_sum_w() },
  m_incident_dependence{ incident_dependence },
  m_bc_time_start{ t_start} , 
  m_bc_time_end{t_end},
  m_abs_glance_angle{ fabs(angular_quadrature.most_glance_mu() ) },
  m_abs_normal_angle{ fabs(angular_quadrature.most_normal_mu() ) },
  m_planck(angular_quadrature.get_number_of_groups() , 0. )
{
  /// need an energy current that we want in and a temperature to evalaute the planckian, if this indeeed a planckia energy distribtuion
  std::cerr << "This feature is not coded yet.  Need to add input varialbes that clarify energy distribution of incident intensity\n";
  exit(EXIT_FAILURE);
}

double Transport_BC_MF_Current::get_boundary_condition(const double mu, const int grp , const double time) 
{
  double val = 0.;
  
  
  return val;
}

