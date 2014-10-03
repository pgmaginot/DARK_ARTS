#ifndef V_Sweep_Fixed_Source_h
#define V_Sweep_Fixed_Source_h

#include "Fem_Quadrature.h"
#include "Eigen/Dense"

/** @file   V_Sweep_Fixed_Source.h
  *   @author pmaginot
  *   @brief provide function interfaces for within group radiation solvers
  *   concrete V_WGRS objects WGRS_FP_Sweeps, WGRS_FP_DSA, WGRS_Krylov_Sweeps, WGRS_Krylov_DSA solve the within group
  *     scattering problem for a given iterate of absorption/re-emission source
  *  Angular_Quadrature is no longer const because the getting boundary conditions may modify some local variables (for the time bein)
  * at some point, (after boundary conditions are initialized, this probaly can be changed back to being a const object
 */

class V_Sweep_Fixed_Source
{
public:
  V_Sweep_Fixed_Source(const Fem_Quadrature& fem_quadrature);
    
  virtual ~V_Sweep_Fixed_Source(){}

  virtual void get_source(Eigen::VectorXd& source_vec, const int dir) = 0;
protected:
  const int m_n_dfem_pts;
  
};

#endif