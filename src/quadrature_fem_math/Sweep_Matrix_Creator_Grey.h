#ifndef Sweep_Matrix_Creator_Grey_h
#define Sweep_Matrix_Creator_Grey_h

#include "V_Sweep_Matrix_Creator.h"

/** @file   Sweep_Matrix_Creator_Grey.h
  *   @author pmaginot
  *   @brief Provide a base class that constructs mass, reaction, and gradient matrices, as well as upwind vectors
  *     Concrete cases for  SELF_LUMPING , TRAD_LUMPING , and EXACT integration techniques
 */

class Sweep_Matrix_Creator_Grey : public V_Sweep_Matrix_Creator
{
public:
  /// Only able to initialize if given an Input_Reader object
  /// constructor defined in Fem_Quadrature.cc
  Sweep_Matrix_Creator_Grey(const Fem_Quadrature& fem_quadrature, Materials* const materials,
    const int n_stages, const double sn_w);
  virtual ~Sweep_Matrix_Creator_Grey(){}
  
  /// calculate \f$ \mathbf{R}_{C_v}^{-1} \f$, \f$ \mathbf{M} \f$, get \f$ \vec{T}^*,~\vec{T}_n \f$
  void update_cell_dependencies(const int cell) override;
  
  /**
    get \f$ \vec{\widehat{B}}_g, ~\mathbf{D}_G^* \f$
    calculate \f$ \bar{\bar{\mathbf{R}}}_{\sigma_{t,g}} \f$ , the isotropic components of \f$ \vec{S}_I \f$, and 
    if grey \f$ \mathbf{R}_{\sigma_{s,0}} + \bar{\bar{\nu}}\mathbf{R}_{\sigma_a} \f$ 
  */  
  void update_group_dependencies(const int grp) override;
   
  void get_s_i(Eigen::VectorXd& s_i, const int dir) override;
  
private:
  
};

#endif