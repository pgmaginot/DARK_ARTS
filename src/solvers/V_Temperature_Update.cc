#include "V_Temperature_Update.h"

V_Temperature_Update::V_Temperature_Update(const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* material)
  :
  m_np(fem_quadrature.get_number_of_interpolation_points() )
{
  std::cout << "Matrix that is saved: \n" << m_r_sig_a ;

}

void V_Temperature_Update::update_temperature(const Intensity_Data& intensity, 
  Temperature_Data& t_star, Temperature_Data& t_n) 
{
  return;
}