#include "Intensity_Update_Grey.h"

Intensity_Update_Grey::Intensity_Update_Grey(const Input_Reader& input_reader, 
  const Fem_Quadrature& fem_quadrature, 
  const Cell_Data& cell_data, 
  Materials& materials, 
  const Angular_Quadrature& angular_quadrature, 
  const int n_stages, 
  const Temperature_Data& t_old, 
  const Intensity_Data& i_old,
  const K_Temperature& kt, 
  K_Intensity& ki,
  const Temperature_Data& t_star)
  :
  V_Intensity_Update(input_reader, fem_quadrature,cell_data,materials,angular_quadrature, n_stages,t_old, i_old, kt, ki, t_star)
{
 
}

void Intensity_Update_Grey::update_intensity(Intensity_Moment_Data& phi)
{
  /**
    For the grey problem, the entire solve is carried out by the V_WGRS solver
  */  
  m_within_group_radiation_solver->solve(phi);
  
  return;
}