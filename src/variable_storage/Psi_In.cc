#include "Psi_In.h"

Psi_In::Psi_In(const int n_groups, const int n_dir)
:
  m_n_groups{n_groups},
  m_n_dir{n_dir}
{
  m_outflows.resize(m_n_groups*m_n_dir,0.);
}

double& Psi_In::operator()(const int group, const int dir)
{
  return m_outflows[location(group,dir)];
}

int Psi_In::location(const int group, const int dir) const
{
  if(group > m_n_groups)
    throw Dark_Arts_Exception(TRANSPORT, "Trying to access Psi_In object group that exceeds m_n_groups");
    
  if(dir > m_n_dir)
    throw Dark_Arts_Exception(TRANSPORT, "Trying to access Psi_In object direction that exceeds m_n_dir");
    
  return group*m_n_dir + dir;
}