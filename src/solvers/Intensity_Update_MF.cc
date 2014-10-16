#include "Intensity_Update_MF.h"

Intensity_Update_MF::Intensity_Update_MF(const Input_Reader& input_reader, 
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
  std::vector<double>& phi_ref_norm)
  :
  V_Intensity_Update(input_reader, fem_quadrature,cell_data,materials, angular_quadrature,n_stages ,t_old, i_old, kt, ki,t_star, phi_ref_norm)
{
  /// initialize the ARD solver
  

}


void Intensity_Update_MF::update_intensity(Intensity_Moment_Data& phi)
{
  /// clear out phi
  
  /** there are options for solving the MF absorption re-emission problem
    1) fixed point, no diffusion operator acceleration
    2) fixed point, with diffusion operator accleration (LMFGA)
    3) krylov with LMFG acceleration
  */
  /// use m_within_group_radiation_solver to get a new iterate for phi given an existing iterate for ard_phi
  

  return;
}
