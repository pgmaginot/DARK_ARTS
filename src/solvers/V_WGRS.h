#ifndef V_WGRS_h
#define V_WGRS_h

#include "Input_Reader.h"
#include "Transport_Sweep.h"

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
    const Fem_Quadrature& fem_quadrature, 
    const Cell_Data& cell_data, 
    Materials& materials, 
    const Angular_Quadrature& angular_quadrature, 
    const int n_stages, 
    const Temperature_Data& t_old, 
    const Intensity_Data& i_old,
    const K_Temperature& kt, 
    K_Intensity& ki,
    const Temperature_Data& t_star,
    std::vector<double>& phi_ref_norm);
    
  virtual ~V_WGRS(){}

  virtual void solve(Intensity_Moment_Data& phi_new) = 0;

  void set_ard_phi_ptr(Intensity_Moment_Data* ard_phi_ptr);
  
  virtual void set_time_data( const double dt, const double time_stage, const std::vector<double>& rk_a_of_stage_i, const int stage ) = 0;
  
  void sweep_for_k_i(K_Intensity& k_i, Intensity_Moment_Data& ard_phi);
protected:
  Transport_Sweep m_transport_sweep;
  
  std::vector<double>& m_phi_ref_norm;
};

#endif