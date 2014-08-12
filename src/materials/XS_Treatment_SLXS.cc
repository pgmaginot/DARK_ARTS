/** @file   XS_Treatment_SLXS.cc
  *   @author pmaginot
  *   @brief Implement the XS_Treatment_SLXS class
  *   xs evaluation points are the dfem itnegration points
*/
#include "XS_Treatment_SLXS.h"

XS_Treatment_SLXS::XS_Treatment_SLXS(const Fem_Quadratre& fem_quad) 
{
  m_n_xs_evals = fem_quad.get_number_of_xs_point();
}

XS_Treatment_SLXS::~XS_Treatment_SLXS()
{

}

void XS_Treatment_SLXS::calculate_xs_at_integration_points(
  const std::vector<double>& xs_evaluations, std::vector<double>& xs_at_dfem_integration_points)
{  
  for(int i=0; i < m_n_xs_evals; i++)
    xs_at_dfem_integration_points[i] = xs_evaluations[i];
  
  return;
}
