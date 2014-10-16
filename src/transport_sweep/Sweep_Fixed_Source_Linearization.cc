#include "Sweep_Fixed_Source_Linearization.h"

Sweep_Fixed_Source_Linearization::Sweep_Fixed_Source_Linearization(const Fem_Quadrature& fem_quadrature, std::shared_ptr<V_Sweep_Matrix_Creator> sweep_matrix_ptr)
:
V_Sweep_Fixed_Source(fem_quadrature),
 m_sweep_matrix_ptr{sweep_matrix_ptr}
{

}

void Sweep_Fixed_Source_Linearization::get_source(Eigen::VectorXd& source_vec)
{
  m_sweep_matrix_ptr->get_s_i(source_vec);
  
  return;
}

