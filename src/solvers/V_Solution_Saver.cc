#include "V_Solution_Saver.h"

V_Solution_Saver::V_Solution_Saver(const Fem_Quadrature& fem_quadrature)
:
m_n_dfem_pts{ fem_quadrature.get_number_of_interpolation_points() }
{
 
}


