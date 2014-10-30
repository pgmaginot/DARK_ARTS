#ifndef Transport_BC_Grey_Planckian_h
#define Transport_BC_Grey_Planckian_h

#include "V_Transport_BC.h"
#include "Materials.h"
#include "Angular_Quadrature.h"


/** @file   V_Transport_BC.h
  *   @author pmaginot
  *   @brief provide function interfaces for within group radiation solvers
  *   concrete V_WGRS objects WGRS_FP_Sweeps, WGRS_FP_DSA, WGRS_Krylov_Sweeps, WGRS_Krylov_DSA solve the within group
  *     scattering problem for a given iterate of absorption/re-emission source
  *  Angular_Quadrature is no longer const because the getting boundary conditions may modify some local variables (for the time bein)
  * at some point, (after boundary conditions are initialized, this probaly can be changed back to being a const object
 */

class Transport_BC_Grey_Planckian : public V_Transport_BC
{
public:
  Transport_BC_Grey_Planckian(Materials& materials, 
    const Angular_Quadrature& angular_quadrature,
    BC_ANGLE_DEPENDENCE incident_dependence,
    const double t_start,
    const double t_end,
  const double boundary_t);
    
  virtual ~Transport_BC_Grey_Planckian(){}

  double get_boundary_condition(const double mu, const int grp, const double time) override;
protected:
  const double m_sn_w;
  const BC_ANGLE_DEPENDENCE m_incident_dependence;
  const double m_bc_time_start;
  const double m_bc_time_end;
  const double m_abs_glance_angle;
  const double m_abs_normal_angle;  
  const double m_planck;
};

#endif