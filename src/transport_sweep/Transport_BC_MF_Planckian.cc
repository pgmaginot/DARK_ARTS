#include "Transport_BC_MF_Planckian.h"



Transport_BC_MF_Planckian::Transport_BC_MF_Planckian(Materials& materials, 
  const Angular_Quadrature& angular_quadrature,
  const BC_ANGLE_DEPENDENCE incident_dependence,
  const double t_start,
  const double t_end,
  const double boundary_temp,
  const BC_ENERGY_DEPENDENCE e_dependence)
:
  V_Transport_BC(),
  m_incident_dependence{ incident_dependence },
  m_bc_time_start{ t_start} , 
  m_bc_time_end{t_end},
  m_abs_glance_angle{ fabs(angular_quadrature.most_glance_mu() ) },
  m_abs_normal_angle{ fabs(angular_quadrature.most_normal_mu() ) },
  m_planck(angular_quadrature.get_number_of_groups() , 0. )
{
  for(int g = 0; g < angular_quadrature.get_number_of_groups() ; g++ )
  {
    m_planck[g] = materials.get_mf_planck(boundary_temp,g);
    if(m_incident_dependence != BC_ISOTROPIC)
      m_planck[g] *= angular_quadrature.get_sum_w();
  }
  
}

double Transport_BC_MF_Planckian::get_boundary_condition(const double mu, const int grp , const double time) 
{
  double val = 0.;
  if( (time > m_bc_time_start) && (time < m_bc_time_end) )
  {   
    if(m_incident_dependence == BC_ISOTROPIC)
    {
      val = m_planck[grp];
    }
    else if(m_incident_dependence == BC_GLANCE)
    {
      /// if not the glancing angle, vacuum BC
      if( fabs( fabs(mu) - m_abs_glance_angle)  < 1.E-6)
        val = m_planck[grp];

    }
    else if(m_incident_dependence == BC_NORMAL)
    {
      /// if not the most normal angle, vacuum BC
      if( fabs( fabs(mu) - m_abs_normal_angle ) < 1.E-6)
        val = m_planck[grp];
      
    }
    else
    {
      throw Dark_Arts_Exception( SUPPORT_OBJECT ,  "Invalid angular dependence for MF_Planckian radiation boundary condition");
    }
  }
  
  return val;
}

