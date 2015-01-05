/**
  A custom data structure that will hold the result of comparing two Intensity_Moment_Data objects
*/
#include "Err_Temperature.h"

Err_Temperature::Err_Temperature(const int n_el)
:
m_np( n_el ),
big_delta( Eigen::VectorXd::Zero(m_np)),
error(0.),
cell(-1),
element(-1),
small_number(1.)
{

}

void Err_Temperature::set_error(const int c, const int el, const double err, Eigen::VectorXd& delta)
{
  cell = c;
  element = el;
  error = err;
  big_delta = delta;
  
  return;
}

void Err_Temperature::clear(void)
{
  cell = -1;
  element = -1;
  error = 0.;
  big_delta = Eigen::VectorXd::Zero(m_np);
  
  return;
}

int Err_Temperature::get_cell_with_worst_err(void) const
{
  return cell;
}

int Err_Temperature::get_element_with_worst_err(void) const
{
  return element;
}
  
double Err_Temperature::get_worst_err(void) const
{
  return error;
}

void Err_Temperature::get_big_delta(Eigen::VectorXd& vec) const
{
  vec = big_delta;
  return;
}

double Err_Temperature::get_small_number(void) const
{
  return small_number;
}

void Err_Temperature::set_small_number(const double val)
{
  small_number = val;
  return;
}



