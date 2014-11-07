#include "Transport_BC_Reflective.h"

Transport_BC_Reflective::Transport_BC_Reflective()
:
  V_Transport_BC()
{
 
}

double Transport_BC_Reflective::get_boundary_condition(const double mu, const int grp , const double time) 
{
  throw Dark_Arts_Exception( SUPPORT_OBJECT , "Reflective condition should be being taken care of elsewhere");
  
  return 0.;
}

