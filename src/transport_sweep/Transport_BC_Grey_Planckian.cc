#include "Transport_BC_Grey_Planckian.h"

Transport_BC_Grey_Planckian::Transport_BC_Grey_Planckian(Materials& materials, 
  const Angular_Quadrature& angular_quadrature,
  BC_ANGLE_DEPENDENCE incident_dependence,
  const double t_start,
  const double t_end,
  const double boundary_t)
:
  V_Transport_BC(),
  m_sn_w{ angular_quadrature.get_sum_w() },
  m_incident_dependence{ incident_dependence },
  m_bc_time_start{ t_start} , 
  m_bc_time_end{t_end},
  m_abs_glance_angle{ fabs(angular_quadrature.most_glance_mu() ) },
  m_abs_normal_angle{ fabs(angular_quadrature.most_normal_mu() ) },
  m_planck{ materials.get_grey_planck(boundary_t) }
{
  
}

double Transport_BC_Grey_Planckian::get_boundary_condition(const double mu, const int grp, const double time) 
{
  double val = 0.;
  if( (time > m_bc_time_start) && (time < m_bc_time_end) )
  {
    
    
    if(m_incident_dependence == BC_ISOTROPIC)
    {
      val = m_planck/m_sn_w;
    }
    else if(m_incident_dependence == BC_GLANCE)
    {
      /// if not the glancing angle, vacuum BC
      if( fabs( fabs(mu) - m_abs_glance_angle)  > 1.E-6)
      {
        val = 0.;
      }
      else
        val = m_planck;

    }
    else if(m_incident_dependence == BC_NORMAL)
    {
      /// if not the most normal angle, vacuum BC
      if( fabs( fabs(mu) - m_abs_normal_angle ) > 1.E-6)
      {
        val = 0.;
      }
      else
        val = m_planck;
      
    }
    else
    {
      std::cerr << "Invalid angular dependence for a radiation boundary condition\n";
      exit(EXIT_FAILURE);
    }
  }
  
  return val;
}

