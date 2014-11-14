/** @file   XS_Treatment_Moment_Preserving.cc
  *   @author pmaginot
  *   @brief Implement the XS_Treatment_Moment_Preserving class
  *   xs evaluation points and legendre polynomials are used to 
*/
#include "XS_Treatment_Moment_Preserving.h"

XS_Treatment_Moment_Preserving::XS_Treatment_Moment_Preserving(const Fem_Quadrature& fem_quad,
  const Input_Reader& input_reader) 
  :
  m_n_xs_evals{ fem_quad.get_number_of_xs_point() },
  m_n_dfem_int_pts{ fem_quad.get_number_of_integration_points() } ,
  m_n_leg_mom{ input_reader.get_opacity_degree() + 1}
{
  /// local quadrature locations to evaluate xs at
  fem_quad.get_xs_eval_points(m_xs_quad_pts);
  /// weights of local quadrature points to evaluate legendre moments of quadrature
  fem_quad.get_xs_eval_weights(m_xs_quad_wts); 
  /// dfem integration points legendre fit of cross section will be mapped to
  fem_quad.get_dfem_integration_points(m_dfem_int_pts);
  
  /// resize vector that stores evaluation of the legendre polynomials at the xs_evaluation_points
  m_leg_poly_at_xs_evals.resize(m_n_xs_evals*m_n_leg_mom , 0.);
  for(int p=0;p < m_n_xs_evals; p++)
  {
    m_leg_poly.get_evaluated_legendre_polynomials(m_xs_quad_pts[p], m_n_leg_mom-1 , p*m_n_leg_mom, m_leg_poly_at_xs_evals);
  }
  
  /// resize the vector that stores the evaluation of the legendre polynomials at the dfem_integration_points
  m_leg_poly_at_dfem_integration_points.resize(m_n_dfem_int_pts*m_n_leg_mom);
  for(int p=0;p < m_n_dfem_int_pts; p++)
  {
    m_leg_poly.get_evaluated_legendre_polynomials(m_dfem_int_pts[p], m_n_leg_mom-1 , p*m_n_leg_mom, m_leg_poly_at_dfem_integration_points);
  }
  
  /// resize coefficents storage vector
  m_leg_coeff_hold.resize(m_n_leg_mom,0.);
}


void XS_Treatment_Moment_Preserving::calculate_xs_at_integration_points(
  const std::vector<double>& xs_evaluations, std::vector<double>& xs_at_dfem_integration_points)
{  
  /** From: http://www.idea.wsu.edu/Quantum/legendre_series.htm
    \f{eqnarray}{
        f(s) & \approx & \sum_{k=0}^{N}{a_k P_k(s) } \\
        a_k &=& \frac{2k+1}{2} \int_{-1}^1{P_k(s) f(s)~ds} \\
            &\approx & \frac{2k+1}{2} \sum_{q=1}^{N_q}{w_q P_k(s_q) f(s_q)} \\
        \sum_{q=1}^{N_q}{w_q} &=& 2
    \f}
  */
  
  /// clear out legednre polynomial coefficient matrix
  for(int n_leg = 0; n_leg < m_n_leg_mom ; n_leg++)
    m_leg_coeff_hold[n_leg] = 0.;
  
  /// loop over quadrature points on the outside (given the layout of m_leg_poly_at_xs_evals) 
  int cnt = 0;
  for(int p=0; p < m_n_xs_evals; p++)
  {
    double w = m_xs_quad_wts[p];
    double f = xs_evaluations[p];
    for(int n_leg = 0; n_leg < m_n_leg_mom; n_leg++)
    {
      m_leg_coeff_hold[n_leg] += w*m_leg_poly_at_xs_evals[cnt]*f;
      cnt++;
    }
  }
  /// normalize legendre moment coefficients
  for(int n_leg = 0; n_leg < m_n_leg_mom; n_leg++)
    m_leg_coeff_hold[n_leg] *= (2.* double(n_leg) + 1.)/2.;
  
  /// we now have the legendre moment coefficients, map onto the dfem integration points
  cnt = 0;
  for(int p=0; p < m_n_dfem_int_pts; p++)
  {
    xs_at_dfem_integration_points[p] = 0.;
    for(int n_leg = 0; n_leg < m_n_leg_mom; n_leg++)
    {
      xs_at_dfem_integration_points[p] += m_leg_coeff_hold[n_leg]*m_leg_poly_at_dfem_integration_points[cnt];
      cnt++;
    }
  }
  
  
  return;
}
