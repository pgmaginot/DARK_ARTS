/**
  A custom data structure that will hold the result of comparing two Intensity_Moment_Data objects
*/
#include "Err_Phi.h"
Err_Phi::Err_Phi(void)
:
error{0.},
cell{-1},
group{-1},
leg_moment{-1}
{

}



void Err_Phi::set_error(int c, int g, int l, double err)
{
  cell = c;
  group = g;
  leg_moment = l;
  error = err;
  return;
}

void Err_Phi::clear(void)
{
  cell = -1;
  group = -1;
  leg_moment = -1;
  error = 0.;
}

int Err_Phi::get_cell_with_worst_err(void) const
{
  return cell;
}

int Err_Phi::get_group_with_worst_err(void) const
{
  return group;
}

int Err_Phi::get_legendre_moment_with_worst_err(void) const
{
  return leg_moment;
}
  
double Err_Phi::get_worst_err(void) const
{
  return error;
}