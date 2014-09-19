#include "Intensity_Update_MF.h"

Intensity_Update_MF::Intensity_Update_MF(const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* materials, 
    const Angular_Quadrature& angular_quadrature, const int n_stages)
  :
  V_Intensity_Update(fem_quadrature,cell_data,materials, angular_quadrature,n_stages )
{

}


void Intensity_Update_MF::update_intensity(const Temperature_Data& t_star, Intensity_Moment_Data& phi)
{
  return;
}
