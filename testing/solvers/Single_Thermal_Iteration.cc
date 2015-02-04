#include <iostream>
#include <vector>
#include <math.h>
#include <string>

#include "Input_Reader.h"
#include "Materials.h"
#include "Fem_Quadrature.h"
#include "Quadrule_New.h"
#include "Cell_Data.h"
#include "Angular_Quadrature.h"
#include "Time_Data.h"
#include "Intensity_Data.h"
#include "Intensity_Moment_Data.h"
#include "Temperature_Data.h"
#include "K_Intensity.h"
#include "K_Temperature.h"
#include "Intensity_Update_Grey.h"
#include "Output_Generator.h"
#include "Err_Temperature.h"
#include "L2_Error_Calculator.h"
#include "Temperature_Update.h"

#include "Dark_Arts_Exception.h"

/**
  Goal of this unit test is to check source moment formation, and reaction matrix for SLXS Lobatto scheme
  -This is a MMS problem with a spatially varying sigma_a, constant cv , zero sig_s
*/ 

int main(int argc, char** argv)
{
  int val = 0;  
  Input_Reader input_reader;    
  try
  {
    input_reader.read_xml(argv[1]);
  }
  catch(const Dark_Arts_Exception& da_exception )
  {
    da_exception.message() ;
    val = -1;
  }       
  
  Quadrule_New quad_fun;  
  Fem_Quadrature fem_quadrature( input_reader , quad_fun);  
  Cell_Data cell_data( input_reader );  
  Angular_Quadrature angular_quadrature( input_reader , quad_fun );   
  L2_Error_Calculator l2_error_calculator(angular_quadrature,fem_quadrature, cell_data,input_reader);
  const int n_cell = cell_data.get_total_number_of_cells();
  const int n_p = fem_quadrature.get_number_of_interpolation_points();  
  
  /// Create a Materials object that contains all opacity, heat capacity, and source objects
  Materials materials( input_reader, fem_quadrature , cell_data, angular_quadrature );  
  Time_Data time_data(input_reader);
  Temperature_Data t_old(fem_quadrature, input_reader, cell_data);  
  Intensity_Data i_old(cell_data, angular_quadrature, fem_quadrature, materials,input_reader);
  
  Intensity_Data i_never_change(cell_data, angular_quadrature, fem_quadrature, materials,input_reader);
  
  const int n_stages = time_data.get_number_of_stages();
  
  K_Temperature kt(n_cell , n_stages, fem_quadrature);
  K_Intensity ki(n_cell ,  n_stages, fem_quadrature, angular_quadrature);  
    
  Intensity_Moment_Data ard_phi(cell_data,angular_quadrature,fem_quadrature,i_old);  

  
  
  try{
    std::string input_filename = argv[1];
    unsigned int found = input_filename.find_last_of("/");
    std::string short_input_filename = input_filename.substr(found+1);  
    Output_Generator output_generator(angular_quadrature,fem_quadrature, cell_data,input_reader);
    
    Temperature_Data t_star( cell_data.get_total_number_of_cells(), fem_quadrature);
  
    std::vector<double> phi_ref_norm;
    ard_phi.get_phi_norm(phi_ref_norm);
        
    std::shared_ptr<V_Intensity_Update> intensity_update;
    intensity_update = std::make_shared<Intensity_Update_Grey>(
        input_reader,fem_quadrature, 
        cell_data,
        materials, 
        angular_quadrature,
        n_stages,
        t_old, 
        i_old,
        kt, 
        ki, 
        t_star, 
        phi_ref_norm ) ;
        
      Temperature_Update temperature_update(fem_quadrature, cell_data, materials, angular_quadrature, n_stages, t_old, ard_phi);
        
    const int stage = 0;

    double time = time_data.get_t_start();
    int t_step = 0;

    double dt = time_data.get_dt(t_step,time);
    
    t_star = t_old;
        
    double time_stage = time + dt*time_data.get_c(stage);
    std::vector<double> rk_a_of_stage_i(n_stages,0.);
    for(int i = 0; i <= stage; i++)
      rk_a_of_stage_i[i] = time_data.get_a(stage,i);
      
    intensity_update->set_time_data(dt,time_stage,rk_a_of_stage_i , stage);
    temperature_update.set_time_data(dt,time_stage,rk_a_of_stage_i , stage);
   
    Err_Temperature err_temperature(n_p);
    err_temperature.set_small_number(1.0E-4*t_old.calculate_average());
 
    int inners = intensity_update->update_intensity(ard_phi);
    temperature_update.update_temperature(t_star, kt, 1.0, err_temperature);
        
    if(inners == 0)
      throw Dark_Arts_Exception(TRANSPORT,"Expecting to take some transport iterations");
        
    intensity_update->calculate_k_i(ki, ard_phi);    
    temperature_update.calculate_k_t(t_star, kt);
       
    ki.advance_intensity(i_old, dt, time_data);
    kt.advance_temperature(t_old, dt, time_data);
  
    std::cout << "dt: " << dt << "time_stage: " << time_stage << "t_step: " << t_step << std::endl;
    std::cout << "Temperature L2 error: " << l2_error_calculator.calculate_l2_error(time_stage, t_old) << std::endl;
    std::cout << "Phi L2 error: " << l2_error_calculator.calculate_l2_error(time_stage , ard_phi) << std::endl;
    time = time + dt;         
  }
  catch(const Dark_Arts_Exception& da)
  {
    val = -1;
    da.testing_message();
  }
  
  

  
  return val;
}
