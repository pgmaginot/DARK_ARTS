#ifndef V_WGRS_h
#define V_WGRS_h

#include "Input_Reader.h"
// #include "Fem_Quadrature.h"
// #include "Cell_Data.h"
// #include "Materials.h"
// #include "Angular_Quadrature.h"
#include "Transport_Sweep.h"
// #include "Intensity_Moment_Data.h"
// #include "Temperature_Data.h"

/** @file   V_WGRS.h
  *   @author pmaginot
  *   @brief provide function interfaces for within group radiation solvers
  *   concrete V_WGRS objects WGRS_FP_Sweeps, WGRS_FP_DSA, WGRS_Krylov_Sweeps, WGRS_Krylov_DSA solve the within group
  *     scattering problem for a given iterate of absorption/re-emission source
  *  Angular_Quadrature is no longer const because the getting boundary conditions may modify some local variables (for the time bein)
  * at some point, (after boundary conditions are initialized, this probaly can be changed back to being a const object
 */

class V_WGRS
{
public:
  V_WGRS(const Input_Reader& input_reader,
    const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* materials, 
    Angular_Quadrature& angular_quadrature, const int n_stages, 
    const Temperature_Data* const t_old, const Temperature_Data* const t_star, 
    const Intensity_Data* const i_old,
    const K_Temperature* const kt, const K_Intensity* const ki);
    
  virtual ~V_WGRS(){}

  virtual void solve(const Temperature_Data* const t_star, Intensity_Moment_Data& phi_new) = 0;

  void set_ard_phi_ptr(Intensity_Moment_Data* ard_phi_ptr);
  
  virtual void set_time_data( const double dt, const double time_stage, const std::vector<double>& rk_a_of_stage_i, const int stage ) = 0;
protected:
  std::shared_ptr<Transport_Sweep> m_transport_sweep;
  
};

#endif