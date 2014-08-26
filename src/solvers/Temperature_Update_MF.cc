#include "Temperature_Update_MF.h"

Temperature_Update_MF::Temperature_Update_MF(const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* material, 
    const Angular_Quadrature& angular_quadrature)
  :
  V_Temperature_Update(fem_quadrature, cell_data,material, angular_quadrature)
{

}

void Temperature_Update_MF::update_temperature(const Intensity_Data& intensity, 
  Temperature_Data& t_new, const Temperature_Data& t_star, const Temperature_Data& t_n,
  const K_Temperature& k_t, const int stage, const double time, const double dt)
{
  return;
}