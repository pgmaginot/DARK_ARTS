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
    const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* materials, 
    Angular_Quadrature& angular_quadrature, const int n_stages);
    
  virtual ~WGRS_FP_Sweeps(){}

  void solve(const Temperature_Data& t_star, Intensity_Moment_Data& phi_new) override;
protected:
  
  
  const int m_n_groups;
  const double m_wg_tolerance;
  const int m_max_sweeps;
  
  Intensity_Moment_Data m_phi_old;
  Err_Phi m_err;

};

#endif