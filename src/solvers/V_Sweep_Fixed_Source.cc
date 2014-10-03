#include "V_Sweep_Fixed_Source.h"

V_Sweep_Fixed_Source::V_Sweep_Fixed_Source(const Fem_Quadrature& fem_quadrature)
:
m_n_dfem_pts{ fem_quadrature.get_number_of_interpolation_points() }
{
 
}


