#include "V_Phi_Error_Calculator.h"

V_Phi_Error_Calculator::V_Phi_Error_Calculator( const Fem_Quadrature& fem_quadrature, 
    const Cell_Data& cell_data,
    const Angular_Quadrature& angular_quadrature)
    :
  m_n_cell(cell_data.get_total_number_of_cells() ),
  m_n_groups(angular_quadrature.get_number_of_groups() ),
  m_n_l_mom( angular_quadrature.get_number_of_leg_moments() ),
  m_n_el(fem_quadrature.get_number_of_interpolation_points() )
{
  
}
