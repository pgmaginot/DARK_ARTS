/** @file   XS_Treatment_Interpolating.cc
  *   @author pmaginot
  *   @brief Implement the XS_Treatment_Interpolating class
  *   xs evaluation points are different from dfem interpolation points
  * evalaute at xs_interpolation points, then calcualte what xs would be at dfem integration points
*/
#include "XS_Treatment_Interpolating.h"

XS_Treatment_Interpolating::XS_Treatment_Interpolating(const Fem_Quadrature& fem_quad) 
{
  m_n_xs_evals = fem_quad.get_number_of_xs_point();
  m_n_integration_pts = fem_quad.get_number_of_integration_points();
  fem_quad.get_xs_at_dfem_integration_points(m_xs_at_dfem_integration_points);
}

XS_Treatment_Interpolating::~XS_Treatment_Interpolating()
{

}

void XS_Treatment_Interpolating::calculate_xs_at_integration_points(
  const std::vector<double>& xs_evaluations, std::vector<double>& xs_at_dfem_integration_points)
{  
  int cnt = 0;
  for(int xs_p=0; xs_p < m_n_xs_evals; xs_p++)
  {
    for(int dfem_int=0; dfem_int < m_n_integration_pts; dfem_int ++)
    {
      xs_at_dfem_integration_points[dfem_int] += xs_evaluations[xs_p]*m_xs_at_dfem_integration_points[cnt];
      cnt++;
    }
  }
  
  return;
}
