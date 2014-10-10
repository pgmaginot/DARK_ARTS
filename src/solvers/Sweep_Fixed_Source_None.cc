#include "Sweep_Fixed_Source_None.h"

Sweep_Fixed_Source_None::Sweep_Fixed_Source_None(const Fem_Quadrature& fem_quadrature)
:
V_Sweep_Fixed_Source(fem_quadrature)
{
 
}

void Sweep_Fixed_Source_None::get_source(Eigen::VectorXd& source_vec)
{
  for(int i=0; i< m_n_dfem_pts ; i++)
    source_vec(i) = 0.;
    
  return;
}

