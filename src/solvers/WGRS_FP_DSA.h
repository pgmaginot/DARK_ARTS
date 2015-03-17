#ifndef WGRS_FP_DSA_h
#define WGRS_FP_DSA_h

#include "V_WGRS.h"
#include "Diffusion_Operator.h"



/** @file   WGRS_FP_DSA.h
  *   @author pmaginot
  *   @brief Solve for within group scattering (and scattering like if this is grey) using fixed point iteration and sweeps only
 */

class WGRS_FP_DSA : public V_WGRS
{
public:
  WGRS_FP_DSA(const Input_Reader& input_reader,
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
    
  virtual ~WGRS_FP_DSA(){}
  
  void kill_petsc_objects() override;

  int solve(Intensity_Moment_Data& phi_new) override;
  
  void set_time_data( const double dt, const double time_stage, const std::vector<double>& rk_a_of_stage_i, const int stage );
protected:  
  const int m_n_groups;
  const double m_wg_tolerance;
  const int m_max_sweeps;
  
  Intensity_Moment_Data m_phi_old;
  Err_Phi m_err;
  
  /// this object creates a MIP diffusion matrix and inverts for a given RHS
  Diffusion_Operator m_diffusion_operator;

};

#endif