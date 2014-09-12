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


  /// calculate \f$ \mathbf{R}_{C_v}^{-1} \f$, \f$ \mathbf{M} \f$, get \f$ \vec{T}^*,~\vec{T}_n \f$
void Sweep_Matrix_Creator_Grey::update_cell_dependencies(const int cell)
{
  /// set cell number
  m_cell_num = cell;
  /// set cell width
  m_dx = m_cell_data->get_cell_width(m_cell_num);
  
  /// calculate \f$ \mathbf{M} \f$ by getting generic mass matrix and multiplying by cell width
  m_dx_div_2_mass = m_dx/2.*m_mass;
  
  /// get \f$ \vec{T}^n \f$
  m_t_old->get_cell_temperature(m_cell_num,m_t_old_vec) ;
  
  /// get \f$ \vec{T}^* \f$
  m_t_star->get_cell_temperature(m_cell_num,m_t_star_vec) ;
  
  /// get Planck vector since it won't change
  m_materials->get_grey_planck(m_t_star_vec, m_planck_vec);
  
  /// load \f$ \mathbf{R}_{C_v} \f$ here
  
 /// then invert it and store in a temporary matrix

  /// put this into 
  
  return;
}

  /**
    get \f$ \vec{\widehat{B}}_g, ~\mathbf{D}_G^* \f$
    calculate \f$ \bar{\bar{\mathbf{R}}}_{\sigma_{t,g}} \f$ , the isotropic components of \f$ \vec{S}_I \f$, and 
    if grey \f$ \mathbf{R}_{\sigma_{s,0}} + \bar{\bar{\nu}}\mathbf{R}_{\sigma_a} \f$ 
  */  
void Sweep_Matrix_Creator_Grey::update_group_dependencies(const int grp)
{
  return;
}