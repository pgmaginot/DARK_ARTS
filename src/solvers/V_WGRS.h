#ifndef V_WGRS_h
#define V_WGRS_h

#include "Input_Reader.h"
#include "Fem_Quadrature.h"
#include "Cell_Data.h"
#include "Materials.h"
#include "Angular_Quadrature.h"
#include "Transport_Sweep.h"
#include "Intensity_Moment_Data.h"
#include "Temperature_Data.h"

/** @file   V_WGRS.h
  *   @author pmaginot
  *   @brief provide function interfaces for within group radiation solvers
  *   concrete V_WGRS objects WGRS_FP_Sweeps, WGRS_FP_DSA, WGRS_Krylov_Sweeps, WGRS_Krylov_DSA solve the within group
  *     scattering problem for a given iterate of absorption/re-emission source
 */

class V_WGRS
{
public:
  V_WGRS(const Input_Reader& input_reader,
    const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* materials, 
    const Angular_Quadrature& angular_quadrature, const int n_stages);
    
  virtual ~V_WGRS(){}

  virtual void solve(const Temperature_Data& t_star, Intensity_Moment_Data& phi_new) = 0;

  void set_ard_phi_ptr(Intensity_Moment_Data* ard_phi_ptr);
protected:
  std::shared_ptr<Transport_Sweep> m_transport_sweep;
  
};

#endif