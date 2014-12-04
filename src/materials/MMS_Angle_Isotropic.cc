#include "MMS_Angle_Isotropic.h"

MMS_Angle_Isotropic::MMS_Angle_Isotropic(const double sum_w) 
:
  V_MMS_Angle(),
  m_sum_w_inv(1./sum_w)
{
  
}

double MMS_Angle_Isotropic::get_angle_component(const int dir)
{ 
  return m_sum_w_inv;
}

double MMS_Angle_Isotropic::get_angle_component(const double mu)
{ 
  return m_sum_w_inv;
}