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

    /// calculate \f$ \mathbf{R}_{C_v}^{-1} \f$, \f$ \mathbf{M} \f$, get \f$ \vec{T}^*,~\vec{T}_n \f$
void Sweep_Matrix_Creator_MF::update_cell_dependencies(const int cell)
{
  return;
}
  
  /**
    get \f$ \vec{\widehat{B}}_g, ~\mathbf{D}_G^* \f$
    calculate \f$ \bar{\bar{\mathbf{R}}}_{\sigma_{t,g}} \f$ , the isotropic components of \f$ \vec{S}_I \f$, and 
    if grey \f$ \mathbf{R}_{\sigma_{s,0}} + \bar{\bar{\nu}}\mathbf{R}_{\sigma_a} \f$ 
  */  

  void Sweep_Matrix_Creator_MF::update_group_dependencies(const int grp)
{
  return;
}
