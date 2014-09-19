#include "Intensity_Update_Grey.h"

Intensity_Update_Grey::Intensity_Update_Grey(const Input_Reader& input_reader,
  Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* materials, 
  const Angular_Quadrature& angular_quadrature, const int n_stages)
  :
  V_Intensity_Update(fem_quadrature,cell_data,materials,angular_quadrature, n_stages)
{

}

void Intensity_Update_Grey::update_intensity(const Temperature_Data& t_star, Intensity_Moment_Data& phi)
{
  /**
    Given a temperature iterate, t_star, find the radiation (angle integrated intensity), phi that would be a result of that temperature (planck and material properties)
    Use a Radiation_Solver to get phi.  Possible Radiation_Solvers are:
    1.  Fixed point, sweeps only
    2.  Fixed point, sweeps + DSA (DSA incorporates fission term)
    3.  Krylov, sweep preconditioning only
    4.  Krylov, sweep + DSA preconditioning
  */
  return;
}