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
#include "Err_Temperature.h"

#include "Dark_Arts_Exception.h"

/**
  Goal of this unit test is to check source moment formation, and reaction matrix for SLXS Lobatto scheme
  -This is a MMS problem with a spatially varying sigma_a, constant cv , zero sig_s
*/ 

int main(int argc, char** argv)
{
  int val = 0;
  const double tol = 1.0E-10;
  
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
  
  const int n_cell = cell_data.get_total_number_of_cells();
  const int n_dir = angular_quadrature.get_number_of_dir();
  const int n_p = fem_quadrature.get_number_of_interpolation_points();
  
  /// Create a Materials object that contains all opacity, heat capacity, and source objects
  Materials materials( input_reader, fem_quadrature , cell_data, angular_quadrature );  
  Time_Data time_data(input_reader);
  Temperature_Data t_old(fem_quadrature, input_reader, cell_data);  
  Intensity_Data i_old(cell_data, angular_quadrature, fem_quadrature, materials,input_reader);
  Temperature_Data t_star(n_cell,fem_quadrature);
  t_star = t_old;
  
  K_Temperature kt(n_cell , time_data.get_number_of_stages(), fem_quadrature);
  K_Intensity ki(n_cell ,  time_data.get_number_of_stages(),fem_quadrature, angular_quadrature);  
    
  Intensity_Moment_Data phi_old(cell_data,angular_quadrature,fem_quadrature,i_old);
  
  std::vector<double> phi_ref_norm;
  phi_old.get_phi_norm(phi_ref_norm);
  
  Intensity_Update_Grey i_update(input_reader, fem_quadrature, cell_data, materials, 
    angular_quadrature, 1,  t_old, i_old, kt, ki, t_star, phi_ref_norm);
  
  Intensity_Moment_Data phi_comparison(phi_old);
    
  try{
    const double dt = time_data.get_dt_max();
    const int stage = 0;
    std::vector<double> rk_a(1,0.);
    rk_a[0] = time_data.get_a(stage,stage);
    const double time_eval = time_data.get_t_start() + dt;
    
    i_update.set_time_data( dt, time_eval, rk_a, stage );
    i_update.update_intensity(phi_old);
    
    /// verify that phi_old[time^(n+1)] = phi[time_old^n] for constant in time problem
    const int grp = 0;    
    Eigen::VectorXd phi_new_loc = Eigen::VectorXd::Zero(n_p);
    Eigen::VectorXd phi_old_loc = Eigen::VectorXd::Zero(n_p);
    for(int cell = 0 ; cell < n_cell; cell++)
    {
      phi_comparison.get_cell_angle_integrated_intensity(cell,grp,0,phi_new_loc);
      phi_old.get_cell_angle_integrated_intensity(cell,grp,0,phi_old_loc);
      std::cout << "Phi_new:\n" << phi_new_loc << "\nPhi_old:\n" << phi_old_loc << std::endl;
      for(int el = 0 ; el < n_p ; el++)
      {
        if(fabs(phi_old_loc(el) - phi_new_loc(el)) > tol)
          throw Dark_Arts_Exception(TRANSPORT, "Intensity Update not working correctly");
      }
    }
    i_update.calculate_k_i(ki, phi_old);
    
    Eigen::VectorXd k_i_loc = Eigen::VectorXd::Zero(n_p);
    
    for(int cell = 0 ; cell < n_cell; cell++)
    {
      for(int dir = 0; dir < n_dir ; dir++)
      {
        ki.get_ki(cell, grp, dir, stage, k_i_loc);        
        std::cout << "K_i for stage: " << stage << " cell: " << cell << " direction: " << dir << " group: " << grp << std::endl << k_i_loc << std::endl; 
        
        for(int el = 0; el < n_p ; el++)
        {
          if( fabs(k_i_loc(el)) > tol)          
            throw Dark_Arts_Exception(TRANSPORT, "Expecting zero k_i");   
        }
      }        
    }
    
    ki.advance_intensity(i_old,dt,time_data);
    
    Intensity_Data i_repeat(cell_data, angular_quadrature, fem_quadrature, materials,input_reader);
    Eigen::VectorXd i_old_loc = Eigen::VectorXd::Zero(n_p);
    Eigen::VectorXd i_new_loc = Eigen::VectorXd::Zero(n_p);
    for(int cell = 0; cell < n_cell ; cell++)
    {
      for(int dir = 0; dir < n_dir ; dir++)
      {
        i_old.get_cell_intensity(cell,grp,dir,i_old_loc);
        i_repeat.get_cell_intensity(cell,grp,dir,i_new_loc);
        for(int el = 0 ; el < n_p ; el++)
        {
          if(fabs(i_old_loc(el) - i_new_loc(el))>tol)
            throw Dark_Arts_Exception(TRANSPORT, "Not advancing constant in time intensity");
        }
      }
    }
    
  }
  catch(const Dark_Arts_Exception& da)
  {
    val = -1;
    da.testing_message();
  }
  
  

  
  return val;
}
