#include "Transport_BC_Vacuum.h"

Transport_BC_Vacuum::Transport_BC_Vacuum()
:
  V_Transport_BC()
{
 
}

double Transport_BC_Vacuum::get_boundary_condition(const double mu, const int grp , const double time) 
{
  return 0.;
}

