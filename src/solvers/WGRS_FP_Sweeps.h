#ifndef WGRS_FP_Sweeps_h
#define WGRS_FP_Sweeps_h

#include "V_WGRS.h"

/** @file   WGRS_FP_Sweeps.h
  *   @author pmaginot
  *   @brief Solve for within group scattering (and scattering like if this is grey) using fixed point iteration and sweeps only
 */

class WGRS_FP_Sweeps : public V_WGRS
{
public:
  WGRS_FP_Sweeps(const Input_Reader& input_reader,
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
    
  virtual ~WGRS_FP_Sweeps(){}

  void kill_petsc_objects() override;
  
  int solve(Intensity_Moment_Data& phi_new, bool& update_success) override;
protected:
  
  
  const int m_n_groups;
  const double m_wg_tolerance;
  const int m_max_sweeps;
  
  Intensity_Moment_Data m_phi_old;
  Err_Phi m_err;

};

#endif