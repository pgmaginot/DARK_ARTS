#include "Time_Marcher.h"

Time_Marcher::Time_Marcher(const Input_Reader&  input_reader, const Angular_Quadrature& angular_quadrature,
    const Fem_Quadrature& fem_quadrature, const Cell_Data& cell_data, Materials& materials,  
    Temperature_Data& t_old, Intensity_Data& i_old,
    const Time_Data& time_data)
    :
    m_n_stages{time_data.get_number_of_stages()},
    m_time_data( time_data),
    m_k_i( cell_data.get_total_number_of_cells() ,m_n_stages, fem_quadrature, angular_quadrature),
    m_k_t( cell_data.get_total_number_of_cells(), m_n_stages, fem_quadrature),
    m_t_star( cell_data.get_total_number_of_cells(), fem_quadrature),
    m_ard_phi( cell_data, angular_quadrature, fem_quadrature, i_old ),
    m_damping{1.},
    m_err_temperature( fem_quadrature.get_number_of_interpolation_points() )
{   
  std::vector<double> phi_ref_norm;
  m_ard_phi.get_phi_norm(phi_ref_norm);
  if( angular_quadrature.get_number_of_groups() > 1){
    m_intensity_update = std::shared_ptr<V_Intensity_Update> (new Intensity_Update_MF(
      input_reader, 
      fem_quadrature, 
      cell_data, 
      materials, 
      angular_quadrature,
      m_n_stages, 
      t_old, 
      i_old,
      m_k_t, 
      m_k_i, 
      m_t_star, 
      phi_ref_norm ) );
    m_temperature_update = std::shared_ptr<V_Temperature_Update> (new Temperature_Update_MF(fem_quadrature, cell_data, materials, angular_quadrature, m_n_stages) );
  }
  else{
    m_intensity_update = std::shared_ptr<V_Intensity_Update> (new Intensity_Update_Grey(
      input_reader,fem_quadrature, 
      cell_data,
      materials, 
      angular_quadrature,
      m_n_stages,
      t_old, 
      i_old,
      m_k_t, 
      m_k_i, 
      m_t_star, 
      phi_ref_norm ) );
    m_temperature_update = std::shared_ptr<V_Temperature_Update> (new Temperature_Update_Grey( fem_quadrature, cell_data, materials, angular_quadrature, m_n_stages ) );
  }
}

void Time_Marcher::solve(Intensity_Data& i_old, Temperature_Data& t_old, Time_Data& time_data)
{
  double time = time_data.get_t_start();
  
  int max_step = int( (time_data.get_t_end() - time_data.get_t_start() )/time_data.get_dt_min() );
  
  double dt = 0.;
  
  /// set here in order to get to compile on a Firday afternoon
  int max_thermal_iter = 100;
  
  std::vector<double> rk_a_of_stage_i(m_n_stages,0.);
  
  for(int t_step = 0; t_step < max_step; t_step++)
  {
    dt = time_data.get_dt(t_step);
    
    /// initial temperature iterate guess is t_old
    m_t_star = t_old;
    m_damping = 1.;
    for(int stage = 0; stage < m_n_stages ; stage++)
    {
      double time_stage = time + dt*time_data.get_c(stage);
      for(int i = 0; i< stage; i++)
        rk_a_of_stage_i[i] = time_data.get_a(stage,i);
        
      /// set time (of this stage), dt (of the whole time step), rk_a for this stage
      // m_intensity_update->set_time_data( dt, time_stage, rk_a_of_stage_i, stage );
      m_temperature_update->set_time_data( dt, time_stage, rk_a_of_stage_i, stage );
      
      for(int therm_iter = 0; therm_iter < max_thermal_iter; therm_iter++)
      {
        /// converge the thermal linearization
        /// first get an intensity given the temperature iterate
        /// Intensity_Update objects are linked to m_star at construction
        m_intensity_update->update_intensity(m_ard_phi);
          
        /// then update temperature given the new intensity
        /// give a damping coefficient to possibly control this Newton (like) iteration)
        /// automatically overrwrite m_t_star, delta / error info tracked in m_temperature_err
        m_temperature_update->update_temperature(m_ard_phi, m_t_star, t_old, m_k_t, m_damping, m_err_temperature );       

      }    
      /** calculate k_I and k_T
       * our intensity and update objects were initialized with const ptr to m_k_i and m_k_t respectively,
       * but let's just call the calculate k_i and k_t functions by passing references to the objects we want to change
       * this will then imply that the more frequently called update functions are not modifying the 
      */
      /// give the converged \f$ \Phi \f$ so that all we have to do is sweep once to get m_k_i
      m_intensity_update->calculate_k_i(m_k_i, m_ard_phi);
      // m_temperature_update->calculate_k_t(&m_t_star, m_k_t, m_ard_phi);
    }
    /// advance to the next time step, overwrite t_old
    time += dt;
    /// these are the only functions that change I_old and T_old
    m_k_i.advance_intensity(i_old,dt,m_time_data);
    m_k_t.advance_temperature(t_old,dt,m_time_data);
    
  }

  return;
}
