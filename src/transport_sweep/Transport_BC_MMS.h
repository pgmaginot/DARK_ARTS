#ifndef Transport_BC_MMS_h
#define Transport_BC_MMS_h

#include "V_Transport_BC.h"
#include "Angular_Quadrature.h"
#include "MMS_Intensity.h"

/** @file   V_Transport_BC.h
  *   @author pmaginot
  *   @brief provide function interfaces for within group radiation solvers
  *   concrete V_WGRS objects WGRS_FP_Sweeps, WGRS_FP_DSA, WGRS_Krylov_Sweeps, WGRS_Krylov_DSA solve the within group
  *     scattering problem for a given iterate of absorption/re-emission source
  *  Angular_Quadrature is no longer const because the getting boundary conditions may modify some local variables (for the time bein)
  * at some point, (after boundary conditions are initialized, this probaly can be changed back to being a const object
 */

class Transport_BC_MMS : public V_Transport_BC
{
public:
  Transport_BC_MMS(const Angular_Quadrature& angular_quadrature,
    const Input_Reader& input_reader, const double x_b);
    
  virtual ~Transport_BC_MMS(){}

  double get_boundary_condition(const double mu, const int grp, const double time) override;
protected:
  const double m_x_boundary;
  MMS_Intensity m_intensity;
};

#endif