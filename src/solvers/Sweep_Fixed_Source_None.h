#ifndef Sweep_Fixed_Source_None_h
#define Sweep_Fixed_Source_None_h

#include "V_Sweep_Fixed_Source.h"

/** @file   Sweep_Fixed_Source_None.h
  *   @author pmaginot
  *   @brief provide function interfaces for within group radiation solvers
  *   concrete V_WGRS objects WGRS_FP_Sweeps, WGRS_FP_DSA, WGRS_Krylov_Sweeps, WGRS_Krylov_DSA solve the within group
  *     scattering problem for a given iterate of absorption/re-emission source
  *  Angular_Quadrature is no longer const because the getting boundary conditions may modify some local variables (for the time bein)
  * at some point, (after boundary conditions are initialized, this probaly can be changed back to being a const object
 */

class Sweep_Fixed_Source_None : public V_Sweep_Fixed_Source
{
public:
  Sweep_Fixed_Source_None(const Fem_Quadrature& fem_quadrature);
    
  virtual ~Sweep_Fixed_Source_None(){}

  void get_source(Eigen::VectorXd& source_vec, const int dir) override;
protected:
  
};

#endif