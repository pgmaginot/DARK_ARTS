#include "V_WGRS.h"

V_WGRS::V_WGRS(const Input_Reader& input_reader,
  const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* materials, 
  Angular_Quadrature& angular_quadrature, const int n_stages, 
    const Temperature_Data* const t_old, 
    const Intensity_Data* const i_old,
    const K_Temperature* const kt, const K_Intensity* const ki)
{
  m_transport_sweep = std::shared_ptr<Transport_Sweep> ( new Transport_Sweep(fem_quadrature, cell_data, materials,angular_quadrature,n_stages,
    t_old, i_old, kt, ki) );
}

void V_WGRS::set_ard_phi_ptr(Intensity_Moment_Data* ard_phi_ptr)
{
  m_transport_sweep->set_ard_phi_ptr(ard_phi_ptr);
  return;
}

void V_WGRS::sweep_for_k_i(const Temperature_Data* t_star, K_Intensity& k_i, Intensity_Moment_Data& ard_phi)
{
  /// set t_star
  m_transport_sweep->set_t_star(t_star);
  /// set other indicators (possibly some matrix creators) to note that we are getting k_I and not sweeping
  /// ---- done in transport_sweep::sweep_mesh() when Solution_Saver is changed
  
  /// sweep the mesh and form k_I on the fly
  /// since we only have accees to one phi at this point (ard_phi, transport sweep needs to be modified have non co
  m_transport_sweep->sweep_mesh(ard_phi, ard_phi , false, true);
  
  return;
}