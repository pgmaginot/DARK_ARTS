#include "Transport_BC_Grey_Current.h"

Transport_BC_Grey_Current::Transport_BC_Grey_Current(const Angular_Quadrature& angular_quadrature,
    BC_ANGLE_DEPENDENCE incident_dependence,
    const double t_start,
    const double t_end,
    const double energy_current)
:
  V_Transport_BC(),
  m_incident_dependence{ incident_dependence },
  m_bc_time_start{ t_start} , 
  m_bc_time_end{t_end},
  m_abs_glance_angle{ fabs(angular_quadrature.most_glance_mu() ) },
  m_abs_normal_angle{ fabs(angular_quadrature.most_normal_mu() ) }
{
  double half_sum_mu_w = 0.;
  if(m_incident_dependence == BC_ISOTROPIC)
  {
    
    for(int d = 0 ; d < angular_quadrature.get_number_of_dir()/2 ; d++)
      half_sum_mu_w += fabs( angular_quadrature.get_mu(d)*angular_quadrature.get_w(d) ) ;      
    
  }
  else if(m_incident_dependence == BC_GLANCE)
  {
    half_sum_mu_w = fabs( angular_quadrature.get_mu(0)*angular_quadrature.get_w(0) );
  }
  else if(m_incident_dependence == BC_NORMAL)
  {
    half_sum_mu_w = fabs( angular_quadrature.get_mu( angular_quadrature.get_number_of_dir()/2 )*angular_quadrature.get_w( angular_quadrature.get_number_of_dir()/2 ) );
  }
  
  m_non_zero_value = energy_current/half_sum_mu_w;
  
}

double Transport_BC_Grey_Current::get_boundary_condition(const double mu, const int grp, const double time) 
{
  double val = 0.;
  if( (time > m_bc_time_start) && (time < m_bc_time_end) )
  {    
    if(m_incident_dependence == BC_ISOTROPIC)
    {
      val = m_non_zero_value;
    }
    else if(m_incident_dependence == BC_GLANCE)
    {
      if( fabs( fabs(mu) - m_abs_glance_angle ) < 1.E-6)
        val = m_non_zero_value;

    }
    else if(m_incident_dependence == BC_NORMAL)
    {
      /// if not the most normal angle, vacuum BC
      if( fabs( fabs(mu) - m_abs_normal_angle ) < 1.E-6)
        val = m_non_zero_value;
    }   
    else
    {
      throw Dark_Arts_Exception( SUPPORT_OBJECT ,  "Invalid angular dependence for Transport_BC_Grey_Current_Boundary condition");
    }
  }
  
  return val;
}

