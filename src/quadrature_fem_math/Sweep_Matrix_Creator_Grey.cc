/** @file   Sweep_Matrix_Creator_Grey.cc
  *   @author pmaginot
  *   @brief Implement the Sweep_Matrix_Creator_Grey class (make pseudo matrices and sources for transport sweep)
  *    Creates the "funny" matrices that arise in the Planck/temperature linearization when accounting for spatially varying material properties
 */

#include "Sweep_Matrix_Creator_Grey.h"

Sweep_Matrix_Creator_Grey::Sweep_Matrix_Creator_Grey(const Fem_Quadrature& fem_quadrature, Materials* const materials,
    Cell_Data* const cell_data, const int n_stages)
:
  V_Sweep_Matrix_Creator( fem_quadrature, materials, cell_data, n_stages )
{  
  
}
void Sweep_Matrix_Creator_Grey::construct_r_sig_t(void)
{
  /**
    \f[
      \bar{\bar{\mathbf R}}_{\sigma_{t,i}} = \mathbf{R}_{\sigma_t} + \frac{1}{c \Delta t a_{ii} }\mathbf{M}    
    \f]
  */
  m_material->get_sigma_a(0,m_temp_mat_vec);  
  m_mtrx_builder->construct_reaction_matrix(m_r_sig_a,m_temp_mat_vec);
  m_r_sig_a *= dx;
  
  r_sig_t = 1./(m_dt*m_c*m_rk_a[stage])*m_mass;
  
  return;
}

void Sweep_Matrix_Creator_Grey::construct_r_sig_s(void)
{
  return;
}

void Sweep_Matrix_Creator_Grey::construct_s_i(void)
{
  return;
}
