#ifndef XS_Treatment_Moment_Preserving_h
#define XS_Treatment_Moment_Preserving_h

#include "V_XS_Treatment.h"
#include "Fem_Quadrature.h"
#include "Legendre_Poly_Evaluation.h"
#include "Input_Reader.h"
#include <vector>
#include <stdlib.h>
#include <math.h>

class XS_Treatment_Moment_Preserving: public V_XS_Treatment
{
public:
  XS_Treatment_Moment_Preserving(const Fem_Quadrature& fem_quad, const Input_Reader& input_reader);
  virtual ~XS_Treatment_Moment_Preserving(){}

  void calculate_xs_at_integration_points( const std::vector<double>& xs_evaluations,
    std::vector<double>& xs_at_dfem_integration_points) override;
    
private:
  const int m_n_xs_evals;
  const int m_n_dfem_int_pts;
  const int m_n_leg_mom;
  
  
  Legendre_Poly_Evaluation m_leg_poly;
  
  std::vector<double> m_leg_coeff_hold;
  
  std::vector<double> m_xs_quad;
  
  std::vector<double> m_dfem_int_pts;
  std::vector<double> m_leg_poly_at_dfem_integration_points;
  std::vector<double> m_leg_poly_at_xs_evals;
  std::vector<double> m_xs_quad_pts;
  std::vector<double> m_xs_quad_wts;
};

#endif