#include "V_WGRS.h"

V_WGRS::V_WGRS(const Input_Reader& input_reader,
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
  m_transport_sweep(fem_quadrature, cell_data, materials,angular_quadrature,n_stages, t_old, i_old, kt, ki,t_star,input_reader) ,
  m_phi_ref_norm(phi_ref_norm)
{

}

void V_WGRS::set_ard_phi_ptr(Intensity_Moment_Data* ard_phi_ptr)
{
  m_transport_sweep.set_ard_phi_ptr(ard_phi_ptr);
  return;
}

void V_WGRS::sweep_for_k_i(K_Intensity& k_i, Intensity_Moment_Data& ard_phi)
{
  /// set other indicators (possibly some matrix creators) to note that we are getting k_I and not sweeping
  /// ---- done in transport_sweep::sweep_mesh() when Solution_Saver is changed
  
  /// sweep the mesh and form k_I on the fly
  /// since we only have accees to one phi at this point (ard_phi, transport sweep needs to be modified have non co
  m_transport_sweep.sweep_mesh(ard_phi, ard_phi , false, true);
  
  return;
}