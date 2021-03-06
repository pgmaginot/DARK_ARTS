#ifndef Sweep_Matrix_Creator_MF_h
#define Sweep_Matrix_Creator_MF_h

#include "V_Sweep_Matrix_Creator.h"

/** @file   Sweep_Matrix_Creator_MF.h
  *   @author pmaginot
  *   @brief Provide a concrete class that matrices and sources for MF radiative transfer sweeps
 */

class Sweep_Matrix_Creator_MF : public V_Sweep_Matrix_Creator
{
public:
  /// Only able to initialize if given an Input_Reader object
  /// constructor defined in Fem_Quadrature.cc
  Sweep_Matrix_Creator_MF(const Fem_Quadrature& fem_quadrature, 
    Materials& materials,
    const int n_stages, 
    const double sn_w, 
    const int n_l_mom,    
    const int n_groups,
    const Temperature_Data& t_old, 
    const Intensity_Data& i_old,
    const K_Temperature& kt, 
    const K_Intensity& ki,
    const Temperature_Data& t_star);
  ~Sweep_Matrix_Creator_MF(){}
  
    /// calculate \f$ \mathbf{R}_{C_v}^{-1} \f$, \f$ \mathbf{M} \f$, get \f$ \vec{T}^*,~\vec{T}_n \f$
  void update_cell_dependencies(const int cell) override;
  
  /**
    get \f$ \vec{\widehat{B}}_g, ~\mathbf{D}_G^* \f$
    calculate \f$ \bar{\bar{\mathbf{R}}}_{\sigma_{t,g}} \f$ , the isotropic components of \f$ \vec{S}_I \f$, and 
    if grey \f$ \mathbf{R}_{\sigma_{s,0}} + \bar{\bar{\nu}}\mathbf{R}_{\sigma_a} \f$ 
  */  
  void update_group_dependencies(const int grp) override;
   
  void update_direction_dependencies(const int dir) override;
  
private:
  const int m_n_groups;
  int m_group_num;
  Eigen::MatrixXd m_spectrum;
  Eigen::VectorXd m_phi_vec;
  Eigen::VectorXd m_sigma_planck;
  Eigen::VectorXd m_local_ard;
  Eigen::VectorXd m_group_independent_xi;
  
};

#endif