#include "Time_Marcher.h"

Time_Marcher::Time_Marcher(const Input_Reader&  input_reader, const Angular_Quadrature& angular_quadrature,
    const Fem_Quadrature& fem_quadrature, const Cell_Data& cell_data, Materials& materials,  
    Temperature_Data& t_old, Intensity_Data& i_old,
    const Time_Data& time_data, std::string input_filename)
    :
    m_n_stages(time_data.get_number_of_stages()),
    m_time_data( time_data),
    m_thermal_tolerance( input_reader.get_thermal_tolerance()  ),
    m_k_i( cell_data.get_total_number_of_cells() ,m_n_stages, fem_quadrature, angular_quadrature),
    m_k_t( cell_data.get_total_number_of_cells(), m_n_stages, fem_quadrature),
    m_t_star( cell_data.get_total_number_of_cells(), fem_quadrature),
    m_ard_phi( cell_data, angular_quadrature, fem_quadrature, i_old ),
    m_damping(1.),
    m_iters_before_damping(input_reader.get_iters_before_damping() ),
    m_damping_decrease_factor(input_reader.get_damping_factor() ),
    m_iteration_increase_factor(input_reader.get_iter_increase_factor() ),
    m_checkpoint_frequency(input_reader.get_restart_frequency() ),
    m_max_damps(input_reader.get_max_damp_iters() ),
    m_err_temperature( fem_quadrature.get_number_of_interpolation_points() ),
    m_status_generator(input_filename),
    m_output_generator(angular_quadrature,fem_quadrature, cell_data, input_filename)
{   
  try{
    std::vector<double> phi_ref_norm;
    m_ard_phi.get_phi_norm(phi_ref_norm);
    if( angular_quadrature.get_number_of_groups() > 1)
    {
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
  catch(const Dark_Arts_Exception& da_exception)
  {
    da_exception.message();
  }
}

void Time_Marcher::solve(Intensity_Data& i_old, Temperature_Data& t_old, Time_Data& time_data)
{
  double time = time_data.get_t_start();
  
  int max_step = int( (time_data.get_t_end() - time_data.get_t_start() )/time_data.get_dt_min() );
  
  double dt = 0.;
  
  /// set here in order to get to compile on a Firday afternoon
  int max_thermal_iter = 1000;
  int times_damped = 0;
  int inners = 0;
  
  
  std::vector<double> rk_a_of_stage_i(m_n_stages,0.);
  
  // m_err_temperature.set_small_number( 1.0E-6*t_old.calculate_average() );
  
  
  for(int t_step = 0; t_step < max_step; t_step++)
  {
    dt = time_data.get_dt(t_step,time);
    
    /// initial temperature iterate guess is t_old
    m_t_star = t_old;
    
    for(int stage = 0; stage < m_n_stages ; stage++)
    {
      m_damping = 1.;  
      
      double time_stage = time + dt*time_data.get_c(stage);
      for(int i = 0; i< stage; i++)
        rk_a_of_stage_i[i] = time_data.get_a(stage,i);
        
      // /// set time (of this stage), dt (of the whole time step), rk_a for this stage
      m_intensity_update->set_time_data( dt, time_stage, rk_a_of_stage_i, stage );
      m_temperature_update->set_time_data( dt, time_stage, rk_a_of_stage_i, stage );
      
      std::cout << "dt:  " << dt << std::endl;
      if(dt < 0.)
        throw Dark_Arts_Exception(TIME_MARCHER, "dt < 0 in time marcher ");
        
      if( (time_stage < time_data.get_t_start()) || (time_stage > time_data.get_t_end() ) )
        throw Dark_Arts_Exception(TIME_MARCHER, "time_stage outside of plausible range");
      
      for(int therm_iter = 0; therm_iter < max_thermal_iter; therm_iter++)
      {
        /// converge the thermal linearization
        /// first get an intensity given the temperature iterate
        /// Intensity_Update objects are linked to m_star at construction        
        inners = m_intensity_update->update_intensity(m_ard_phi);
        std::cout << " Time step: " << t_step << " Stage: " << stage << " Thermal iteration: " << therm_iter << " Number of Transport solves: " << inners << std::endl;
          
        /// then update temperature given the new intensity
        /// give a damping coefficient to possibly control this Newton (like) iteration)
        /// automatically overrwrite m_t_star, delta / error info tracked in m_temperature_err
        m_temperature_update->update_temperature(m_ard_phi, m_t_star, t_old, m_k_t, m_damping, m_err_temperature );       

        double norm_relative_change = m_err_temperature.get_worst_err();
        
        /// write to iteration status file
        m_status_generator.write_iteration_status(t_step, stage, dt , inners , norm_relative_change, m_damping);
        
        /// check convergence of temperature
        if( norm_relative_change < m_thermal_tolerance)
        {
          break;
        }     
        else
        {
          /// damp if necessary 
          if(therm_iter > m_iters_before_damping)
          {
              therm_iter = 0;
              m_iters_before_damping *= m_iteration_increase_factor;
              times_damped++;
              m_damping *= m_damping_decrease_factor;
              m_t_star = t_old;
              if(times_damped > m_max_damps) 
              {
                std::stringstream err;
                err << "Failing to converge thermal iteration.  Exiting now.  Saving Last iterates";
                m_output_generator.write_xml(false,t_step,i_old);
                m_output_generator.write_xml(false,t_step,t_old);
                m_output_generator.write_xml(false,t_step,m_ard_phi);
                throw Dark_Arts_Exception(TIME_MARCHER, err);
              }
          }          
        }       
      }    
      /** calculate k_I and k_T
       * our intensity and update objects were initialized with const ptr to m_k_i and m_k_t respectively,
       * but let's just call the calculate k_i and k_t functions by passing references to the objects we want to change
       * this will then imply that the more frequently called update functions are not modifying the 
      */
      
      /// give the converged \f$ \Phi \f$ so that all we have to do is sweep once to get m_k_i
      /// inorder to integrate spatial-temporal error, we must form Phi_stage and T_stage in these function
      m_intensity_update->calculate_k_i(m_k_i, m_ard_phi);
      m_temperature_update->calculate_k_t(m_t_star, m_k_t, m_ard_phi);
    }
    /// advance to the next time step, overwrite t_old
    time += dt;
    /// these are the only functions that change I_old and T_old
    m_k_i.advance_intensity(i_old,dt,m_time_data);
    m_k_t.advance_temperature(t_old,dt,m_time_data);
    
    if( (t_step % m_checkpoint_frequency) == 0)
    {
      m_output_generator.write_xml(false,t_step,i_old);
      m_output_generator.write_xml(false,t_step,t_old);
      m_output_generator.write_xml(false,t_step,m_ard_phi);    
    }
    
    /// check to see if we are at the end of the time marching scheme
    if( abs( (time - time_data.get_t_end() )/time) < 1.0E-6)
      break;
  }
  
  /// dump final solutions, always!
  m_output_generator.write_xml(true,0,i_old);
  m_output_generator.write_xml(true,0,t_old);
  m_output_generator.write_xml(true,0,m_ard_phi);
  
  return;
}
