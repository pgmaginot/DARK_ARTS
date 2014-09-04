/** @file   Sweep_Matrix_Creator_MF.cc
  *   @author pmaginot
  *   @brief Implement the Sweep_Matrix_Creator_MF class (make pseudo matrices and sources for transport sweep)
  *    Creates the "funny" matrices that arise in the Planck/temperature linearization when accounting for spatially varying material properties
 */

#include "Sweep_Matrix_Creator_MF.h"

Sweep_Matrix_Creator_MF::Sweep_Matrix_Creator_MF(const Fem_Quadrature& fem_quadrature, Materials* const materials,
    Cell_Data* const cell_data, const int n_stages)
:
  V_Sweep_Matrix_Creator( fem_quadrature, materials, cell_data, n_stages )
{  
  
}

/**
  \f[
    \text{m_r_sig_t} = \bar{\bar{\mathbf R}}_{\sigma_{t,i}} = \mathbf{R}_{\sigma_t} + \frac{1}{c \Delta t a_{ii} } \mathbf{M}
  \f]
*/
void Sweep_Matrix_Creator_MF::construct_r_sig_t(void)
{
  return;
}

void Sweep_Matrix_Creator_MF::construct_r_sig_s(void)
{
  return;
}

void Sweep_Matrix_Creator_MF::construct_s_i(void)
{
  return;
}
