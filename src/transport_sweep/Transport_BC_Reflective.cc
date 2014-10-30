#include "Transport_BC_Reflective.h"

Transport_BC_Reflective::Transport_BC_Reflective()
:
  V_Transport_BC()
{
 
}

double Transport_BC_Reflective::get_boundary_condition(const double mu, const int grp , const double time) 
{
  std::cerr << "Accessing reflective boundary condition creator\n";
  std::cerr << "By covnention, this should not be happening\n";
  exit(EXIT_FAILURE);
  return 0.;
}

