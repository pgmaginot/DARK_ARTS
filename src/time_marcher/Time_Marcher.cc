#include "Time_Marcher.h"

Time_Marcher::Time_Marcher(const Input_Reader&  input_reader, Angular_Quadrature& angular_quadrature,
    const Fem_Quadrature& fem_quadrature, Cell_Data* const cell_data, Materials* const materials,  
    Temperature_Data& t_old, Intensity_Data& i_old,
    const Time_Data& time_data)
    :
    m_n_stages{time_data.get_number_of_stages()},
    m_k_i( cell_data->get_total_number_of_cells() ,m_n_stages, fem_quadrature, angular_quadrature),
    m_k_t( cell_data->get_total_number_of_cells(), m_n_stages, fem_quadrature),
    m_t_star( cell_data->get_total_number_of_cells(), fem_quadrature)
{
  if( angular_quadrature.get_number_of_groups() > 1){
    m_intensity_update = std::shared_ptr<V_Intensity_Update> (new Intensity_Update_MF(input_reader, fem_quadrature, cell_data, materials, angular_quadrature, m_n_stages, &t_old, &m_t_star, &i_old, &m_k_t, &m_k_i) );
    m_temperature_update = std::shared_ptr<V_Temperature_Update> (new Temperature_Update_MF(fem_quadrature, cell_data, materials, angular_quadrature, m_n_stages) );
  }
  else{
    m_intensity_update = std::shared_ptr<V_Intensity_Update> (new Intensity_Update_Grey(input_reader,fem_quadrature, cell_data, materials, angular_quadrature, m_n_stages,&t_old, &m_t_star, &i_old, &m_k_t, &m_k_i ) );
    m_temperature_update = std::shared_ptr<V_Temperature_Update> (new Temperature_Update_Grey( fem_quadrature, cell_data, materials, angular_quadrature, m_n_stages ) );
  }
}