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
  const Temperature_Data& t_star,
  std::vector<double>& phi_ref_norm)
  :
  V_Intensity_Update(input_reader, fem_quadrature,cell_data,materials,angular_quadrature, n_stages,t_old, i_old, kt, ki, t_star, phi_ref_norm)
{
 
}

void Intensity_Update_Grey::kill_petsc_objects()
{
  m_within_group_radiation_solver->kill_petsc_objects();
}

int Intensity_Update_Grey::update_intensity(Intensity_Moment_Data& phi, bool& update_success)
{
  /**
    For the grey problem, the entire solve is carried out by the V_WGRS solver, so this looks particluarly emptry.  The MF case is not empty!
  */  
  int inners = m_within_group_radiation_solver->solve(phi,update_success);
  
  return inners;
}