#include "V_WGRS.h"

V_WGRS::V_WGRS(const Input_Reader& input_reader,
  const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* materials, 
  Angular_Quadrature& angular_quadrature, const int n_stages, 
    const Temperature_Data* const t_old, const Temperature_Data* const t_star, 
    const Intensity_Data* const i_old,
    const K_Temperature* const kt, const K_Intensity* const ki)
{
  m_transport_sweep = std::shared_ptr<Transport_Sweep> ( new Transport_Sweep(fem_quadrature, cell_data, materials,angular_quadrature,n_stages,
    t_old, t_star, i_old, kt, ki) );
}

void V_WGRS::set_ard_phi_ptr(Intensity_Moment_Data* ard_phi_ptr)
{
  m_transport_sweep->set_ard_phi_ptr(ard_phi_ptr);
  return;
}

