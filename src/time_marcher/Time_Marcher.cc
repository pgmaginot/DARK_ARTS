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

void Time_Marcher::solve(Intensity_Data& i_old, Temperature_Data& t_old, Time_Data& time_data)
{
  double time = time_data.get_t_start();
  double time_stage = time;
  
  int max_step = int( (time_data.get_t_end() - time_data.get_t_start() )/time_data.get_dt_min() );
  
  double dt = 0.;
  
  /// set here in order to get to compile on a Firday afternoon
  int max_thermal_iter = 100;
  
  std::vector<double> rk_a_of_stage_i(m_n_stages,0.);
  
  for(int t_step = 0; t_step < max_step; t_step++)
  {
    dt = time_data.get_dt(t_step);
    for(int stage = 0; stage < m_n_stages ; stage++)
    {
      time_stage = time + dt*time_data.get_c(stage);
      for(int i = 0; i< stage; i++)
        rk_a_of_stage_i[i] = time_data.get_a(stage,i);
        
      /// set time (of this stage), dt (of the whole time step), rk_a for this stage
      for(int therm_iter = 0; therm_iter < max_thermal_iter; therm_iter++)
      {
        /// converge the thermal linearization
        
      }    
      /// calculate k_I and k_T
      
    }
    /// advance to the next time step
    time = time + dt;
    advance_intensity(i_old,m_k_i);
    advance_temperature(t_old,m_k_t);
    
  }

  return;
}

void Time_Marcher::advance_intensity(Intensity_Data& i_old, K_Intensity& k_i)
{
  return;
}
 
void Time_Marcher::advance_temperature(Temperature_Data& t_old, K_Temperature& k_t)
{
  return;
}