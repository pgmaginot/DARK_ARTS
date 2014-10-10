/** @file   Sweep_Matrix_Creator_MF.cc
  *   @author pmaginot
  *   @brief Implement the Sweep_Matrix_Creator_MF class (make pseudo matrices and sources for transport sweep)
  *    Creates the "funny" matrices that arise in the Planck/temperature linearization when accounting for spatially varying material properties
 */

#include "Sweep_Matrix_Creator_MF.h"

Sweep_Matrix_Creator_MF::Sweep_Matrix_Creator_MF(const Fem_Quadrature& fem_quadrature, 
  Materials& materials,
  const int n_stages, const double sn_w, 
  const int n_l_mom,
  const Temperature_Data& t_old,
  const Intensity_Data& i_old,
  const K_Temperature& kt, 
  const K_Intensity& ki,
  const Temperature_Data& t_star)
:
  V_Sweep_Matrix_Creator( fem_quadrature, materials, n_stages , sn_w, n_l_mom, t_old, i_old, kt, ki,t_star)
{  

}

    /// calculate \f$ \mathbf{R}_{C_v}^{-1} \f$, \f$ \mathbf{M} \f$, get \f$ \vec{T}^*,~\vec{T}_n \f$
void Sweep_Matrix_Creator_MF::update_cell_dependencies(const int cell)
{
  /// set cell number
  m_cell_num = cell;
  
  /// get \f$ \vec{T}^n \f$
  m_t_old_ref.get_cell_temperature(m_cell_num,m_t_old_vec) ;
  
  /// get \f$ \vec{T}^* \f$
  m_t_star_ref.get_cell_temperature(m_cell_num,m_t_star_vec) ;
  
  /// populate Materials object with local temperature and position to evalaute material properties
  m_materials.calculate_local_temp_and_position(cell,m_t_star_vec);
  
  /// get Planck vector since it won't change
  m_materials.get_grey_planck(m_t_star_vec, m_planck_vec);
  
  /// set cell width
  m_dx = m_materials.get_cell_width();
  
  /// calculate \f$ \mathbf{M} \f$ by getting generic mass matrix and multiplying by cell width
  m_dx_div_2_mass = m_dx/2.*m_mass;
  
  m_mass_inv = m_dx_div_2_mass.inverse();
    
  /// load \f$ \mathbf{R}_{C_v} \f$ here
  m_mtrx_builder->construct_r_cv(m_r_cv);
  
 /// then invert it and store in a temporary matrix
  m_coefficient = m_r_cv.inverse(); 
  /// put this back into m_r_cv
  m_r_cv = m_coefficient;
  
  /// we assume that whatever intensity we are updating in the sweep has already calculated an absoprtion rate density
  // m_ard_vec
  
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

void Sweep_Matrix_Creator_MF::update_direction_dependencies(const int dir)
{
  return;
}

