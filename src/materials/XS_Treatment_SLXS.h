#ifndef XS_Treatment_SLXS_h
#define XS_Treatment_SLXS_h

#include "V_XS_Treatment.h"
#include "Fem_Quadrature.h"
#include <vector>
#include <stdlib.h>
#include <math.h>

class XS_Treatment_SLXS: public V_XS_Treatment
{
public:
  XS_Treatment_SLXS(const Fem_Quadrature& fem_quad);
  virtual ~XS_Treatment_SLXS();

  void calculate_xs_at_integration_points( const std::vector<double>& xs_evaluations,
    std::vector<double>& xs_at_dfem_integration_points) override;
    
private:
  int m_n_xs_evals = -1;
};

#endif