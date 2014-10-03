#include "Intensity_Update_Grey.h"

Intensity_Update_Grey::Intensity_Update_Grey(const Input_Reader& input_reader, 
  const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, 
  Materials* materials, Angular_Quadrature& angular_quadrature, const int n_stages, 
    const Temperature_Data* const t_old, const Temperature_Data* const t_star, 
    const Intensity_Data* const i_old,
    const K_Temperature* const kt, const K_Intensity* const ki)
  :
  V_Intensity_Update(input_reader, fem_quadrature,cell_data,materials,angular_quadrature, n_stages,t_old, t_star, i_old, kt, ki)
{
 
}

void Intensity_Update_Grey::update_intensity(const Temperature_Data* const t_star, Intensity_Moment_Data& phi)
{
  /**
    For the grey problem, the entire solve is carried out by the V_WGRS solver
  */  
  m_within_group_radiation_solver->solve(t_star, phi);
  
  return;
}