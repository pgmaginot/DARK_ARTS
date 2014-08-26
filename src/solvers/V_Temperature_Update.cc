#include "V_Temperature_Update.h"

V_Temperature_Update::V_Temperature_Update(const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* material, 
  const Angular_Quadrature& angular_quadrature)
  :
  m_material(material),
  m_cell_data(cell_data),
  m_np(fem_quadrature.get_number_of_interpolation_points() ) , 
  m_matrix_type(fem_quadrature.get_integration_type() ),
  m_n_cells(cell_data->get_total_number_of_cells() ),
  m_sn_w(angular_quadrature.get_sum_w() )
{
 if(m_matrix_type ==  EXACT)
 {
    m_mtrx_builder = std::shared_ptr<V_Matrix_Construction> (new    
      Matrix_Construction_Exact(fem_quadrature) );
 }
 else if(m_matrix_type == TRAD_LUMPING)
 {
    m_mtrx_builder = std::shared_ptr<V_Matrix_Construction> (new 
      Matrix_Construction_Trad_Lumping(fem_quadrature) );
 }
 else if(m_matrix_type == SELF_LUMPING)
 {
    m_mtrx_builder = std::shared_ptr<V_Matrix_Construction> (new 
      Matrix_Construction_Self_Lumping(fem_quadrature) );
 }
 
 m_temp_mat_vec.resize(m_np,0.);
}
