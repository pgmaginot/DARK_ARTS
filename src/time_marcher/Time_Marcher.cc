#include "Time_Marcher.h"

Time_Marcher::Time_Marcher(const Input_Reader&  input_reader, const Angular_Quadrature& angular_quadrature,
    const Fem_Quadrature& fem_quadrature, const Cell_Data& cell_data, Materials& materials,  
    Temperature_Data& t_old, Intensity_Data& i_old,
    const Time_Data& time_data, std::string& stat_filename)
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
    m_damp_trigger_initial(m_iters_before_damping),
    m_damping_decrease_factor(input_reader.get_damping_factor() ),
    m_iteration_increase_factor(input_reader.get_iter_increase_factor() ),
    m_checkpoint_frequency(input_reader.get_restart_frequency() ),
    m_max_damps(input_reader.get_max_damp_iters() ),
    m_max_thermal_iter( input_reader.get_max_thermal_iteration() ),
    m_suppress_output( input_reader.get_output_suppression() ),
    m_err_temperature( fem_quadrature.get_number_of_interpolation_points() ),
    m_status_generator( stat_filename),
    m_output_generator(angular_quadrature,fem_quadrature, cell_data, input_reader),
    m_calculate_space_time_error( input_reader.record_space_time_error() ),
    m_calculate_final_space_error( input_reader.record_final_space_error() ),
    m_temperature_update(fem_quadrature, cell_data, materials, angular_quadrature, m_n_stages, t_old, m_ard_phi),
    m_input_reader(input_reader),
    m_cell_data(cell_data),
    m_angular_quadrature(angular_quadrature),
    recent_iteration_errors(5,0.)
{       
  try{
    std::vector<double> phi_ref_norm;
    m_ard_phi.get_phi_norm(phi_ref_norm);
    if( angular_quadrature.get_number_of_groups() > 1)
    {
      m_intensity_update = std::make_shared<Intensity_Update_MF>(
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
        phi_ref_norm ) ;
    }
    else{
      m_intensity_update = std::make_shared<Intensity_Update_Grey>(
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
        phi_ref_norm ) ;
    }
    
    std::string output_directory = input_reader.get_output_directory();
    output_directory += input_reader.get_filename_base_for_results();
        
    // std::cout << "L2 error should be in " << filename_base1 << std::endl;
    
    if(m_calculate_space_time_error)
      m_space_time_error_calculator = std::make_shared<Space_Time_Error_Calculator>(angular_quadrature,
        fem_quadrature, cell_data,  input_reader, time_data, output_directory);
      
    if(m_calculate_final_space_error)
    {
      m_final_space_error_calculator = std::make_shared<Final_Space_Error_Calculator>(angular_quadrature,
        fem_quadrature, cell_data,  input_reader, time_data,  output_directory);
    }
    
    fem_quadrature.get_dfem_interpolation_point(dfem_interp_points);
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
  
  double dt = time_data.get_dt_min();
  double time_stage = 0.;
  
  int times_damped = 0;
  int total_inners = 0;
  
  int total_thermals = 0;
  int inners = 0;
  int t_step = 0;
  
  std::vector<double> rk_a_of_stage_i(m_n_stages,0.);  
  
  /// use some extra things to see if we are just slowly converging before damping

  bool need_to_cut_dt = false;
  while( true )
  {
    if( (t_step % 10) == 0)
    {
      m_err_temperature.set_small_number( 1.0E-3*t_old.calculate_average() );  
      m_ard_phi.update_error_cutoff();
    }
    if( need_to_cut_dt)
    {
      dt *= m_damping_decrease_factor;
    }
    else
    {
      dt = time_data.get_dt(t_step,time,dt);
    }
    
    need_to_cut_dt = false;
    
    /// initial temperature iterate guess is t_old
    m_t_star = t_old;
    
    /// this variable can be used to cut time step    
    for( int stage = 0 ; stage < m_n_stages ; stage++)
    {    
      int stage_inners = 0;
      int stage_thermals = 0;
      // m_t_star.make_non_zero_guess();
      if(need_to_cut_dt)
      {
        /// get out of this loop, and out to the time step loop
        break;
      }
      
      /// variables to control thermal iteration
      /// damping factor, time damped, reset threshold to apply damping     
      bool converged_thermal = false;      
      m_damping = 1.; 
      times_damped = 0;
      m_iters_before_damping = m_damp_trigger_initial;
      
      time_stage = time + dt*time_data.get_c(stage);
      
      
      // std::cout << "Time: " << time << " dt_full: " << dt << " stage: " << stage << " Time stage: " << time_stage << std::endl;
            
      for(int i = 0; i<= stage; i++)
        rk_a_of_stage_i[i] = time_data.get_a(stage,i);
        
      // // /// set time (of this stage), dt (of the whole time step), rk_a for this stage
      m_intensity_update->set_time_data( dt, time_stage, rk_a_of_stage_i, stage );
      m_temperature_update.set_time_data( dt, time_stage, rk_a_of_stage_i, stage );
              
      // if( (time_stage < time_data.get_t_start()) || (time_stage > time_data.get_t_end() ) )
        // throw Dark_Arts_Exception(TIME_MARCHER, "time_stage outside of plausible range");
      
      for(int therm_iter = 0;  ; therm_iter++)
      {
        
        /// converge the thermal linearization
        /// first get an intensity given the temperature iterate
        /// Intensity_Update objects are linked to m_star at construction        
        bool intensity_update_success = false;
        m_ard_phi.clear_angle_integrated_intensity();
        inners = m_intensity_update->update_intensity(m_ard_phi,intensity_update_success);
        total_inners += inners;
        stage_inners += inners;
        // m_ard_phi.mms_cheat(time_stage,m_cell_data,dfem_interp_points,m_input_reader,m_angular_quadrature);
        
        double norm_relative_change = 0.;
        if(intensity_update_success)
        {
          
          /// then update temperature given the new intensity
          /// give a damping coefficient to possibly control this Newton (like) iteration)
          /// automatically overrwrite m_t_star, delta / error info tracked in m_temperature_err
          m_err_temperature.clear();
          m_temperature_update.update_temperature(m_t_star, m_k_t, m_damping, m_err_temperature ); 
          // m_t_star.mms_cheat(time_stage,m_cell_data,dfem_interp_points,m_input_reader);
          norm_relative_change = m_err_temperature.get_worst_err();
          cascade_thermal_err(norm_relative_change);
          std::cout << " Time step: " << t_step << " Stage: " << stage << " Thermal iteration: " << therm_iter <<
            " Number of Transport solves: " << inners << " Thermal error: " << norm_relative_change << std::endl;
          /// write to iteration status file
          // m_status_generator.write_iteration_status(t_step, stage, therm_iter, dt , inners , norm_relative_change, m_damping);
        }
        else
        {
          total_thermals += therm_iter;
          
          stage_thermals += therm_iter;
          /// restart this time step, with a smaller dt
          std::cout << "Intensity Update Needing to cut dt" << std::endl;
          std::cout << "Converged_thermal is: " << converged_thermal << std::endl;
          converged_thermal = false;
          need_to_cut_dt = true;
          break; /// break into the stage loop.  Really want to break into the time loop
        }
        /// check convergence of temperature
        if( norm_relative_change < m_thermal_tolerance)
        {
          total_thermals += therm_iter;
          stage_thermals += therm_iter;
          converged_thermal = true;
          break;
        }     
        else
        {
          /// damp if necessary 
          if( (therm_iter > m_iters_before_damping) )
          {
            if( determine_if_just_slowly_converging() )
            {
              m_iters_before_damping *= 2;
            }
            else
            {
              // std::cout << "Damping thermal iteration\n" ;    
              total_thermals += therm_iter;    
              stage_thermals +=therm_iter;
              therm_iter = 0;
              m_iters_before_damping = int( ceil( m_iteration_increase_factor*double(m_iters_before_damping) ) );
              // std::cout << "iters before damping again: " << m_iters_before_damping << std::endl;
              times_damped++;
              m_damping *= m_damping_decrease_factor;
              // std::cout << "New damping factor: " << m_damping << std::endl;
              m_t_star = t_old;
              if(times_damped > m_max_damps) 
              {
                std::cout << "Damped too many times" << std::endl;
                std::cout << "Cutting dt and restarting time step" << std::endl;
                need_to_cut_dt = true;
              }
            }
          }          
        } // check of thermal iteration convergence       
      } // bottom of thermal iteration (of this stage) loop    
      
      if(need_to_cut_dt)
      {
        /// exit stage loop, better have a test immediately after to get back to the top of the time loop
        break;
      }
      
      /** calculate k_I and k_T
       * our intensity and update objects were initialized with const ptr to m_k_i and m_k_t respectively,
       * but let's just call the calculate k_i and k_t functions by passing references to the objects we want to change
       * this will then imply that the more frequently called update functions are not modifying the 
      */
      /// give the converged \f$ \Phi \f$ so that all we have to do is sweep once to get m_k_i
      m_intensity_update->calculate_k_i(m_k_i, m_ard_phi);
      m_temperature_update.calculate_k_t(m_t_star, m_k_t);
      
      /// m_ard_phi and m_t_star are the radiation and temperature profiles at this time stage
      if(m_calculate_space_time_error)
        m_space_time_error_calculator->record_error(dt, stage, time_stage, m_ard_phi, m_t_star);
      
      m_status_generator.write_iteration_status(t_step, stage, stage_thermals, dt , stage_inners , 0.0 , m_damping);
    } // bottom of stage loop
    
    if(need_to_cut_dt)
      continue;

    /// advance to the next time step, overwrite t_old
    time += dt;
    /// these are the only functions that change I_old and T_old
    m_k_i.advance_intensity(i_old,dt,m_time_data);
    m_k_t.advance_temperature(t_old,dt,m_time_data);

      
    if(!m_suppress_output)
    {
      if( (t_step % m_checkpoint_frequency) == 0)
      {
        m_output_generator.write_xml(false,t_step,i_old);
        m_output_generator.write_xml(false,t_step,t_old);
        m_output_generator.write_xml(false,t_step,m_ard_phi);    
      }
    }
      
      /// check to see if we are at the end of the time marching scheme
    if( fabs( (time - time_data.get_t_end() )/time) < 1.0E-4)
      break;
      
    t_step++;
  } /// bottom of time step loop (outermost)
  
  std::cout << "Needed: " << total_thermals << " thermal iterations" << std::endl;
  std::cout << "Needed: " << total_inners << " transport iterations" << std::endl;
  
  m_intensity_update->kill_petsc_objects();
  
  /// call for end spatial error only
  /// we were fancy with space_time error and output the error during the destructor call
  if(m_calculate_final_space_error)
  {
    /// need to calculate m_ard_phi from i_old, since m_ard_phi is at the last time stage value
    m_ard_phi.clear_angle_integrated_intensity();
    m_ard_phi.update_phi_and_norms(i_old);
    
    /// t_step is off by 1
    t_step++;
    m_final_space_error_calculator->record_error(time,t_step,t_old, m_ard_phi, m_status_generator.get_total_thermals() , m_status_generator.get_total_sweeps() );
  }
  /// dump final solutions, always!
  if(!m_suppress_output)
  {
    m_output_generator.write_xml(true,0,i_old);
    m_output_generator.write_xml(true,0,t_old);
    m_output_generator.write_xml(true,0,m_ard_phi);
    
    if(m_input_reader.get_output_type() == DUMP)
    {
      m_output_generator.write_txt(true,0,i_old);
      m_output_generator.write_txt(true,0,t_old);
      m_output_generator.write_txt(true,0,m_ard_phi);    
    }
  }
  
  return;
}

void Time_Marcher::cascade_thermal_err(const double last_err)
{
  recent_iteration_errors[4] = recent_iteration_errors[3];
  recent_iteration_errors[3] = recent_iteration_errors[2];
  recent_iteration_errors[2] = recent_iteration_errors[1];
  recent_iteration_errors[1] = recent_iteration_errors[0];
  recent_iteration_errors[0] = last_err;
  
  return;
}

bool Time_Marcher::determine_if_just_slowly_converging(void)
{
  if(m_already_extended)
  {
    m_already_extended = false;
    return false;
  }
  else
  {
    bool reducing = false;
    for(int i = 0 ; i < 4; i++)
    {
      double reduction = (recent_iteration_errors[i+1]/recent_iteration_errors[i]);
      // std::cout << "reduction factor[" << i << "]: " << reduction << std::endl;
      
      reducing = (reduction > 1.0 );
      if(!reducing)
      {
        return false;
      }
    }  
  }
  m_already_extended = true;
  return true;
}
