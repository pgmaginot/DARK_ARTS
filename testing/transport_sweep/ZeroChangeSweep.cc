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
#include "Transport_Sweep.h"
#include "Temperature_Update_Grey.h"
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

  Transport_Sweep transport_sweep(
    fem_quadrature, 
    cell_data, 
    materials, 
    angular_quadrature, 
    time_data.get_number_of_stages(), 
    t_old, 
    i_old,
    kt, 
    ki,
    t_star,
    input_reader);
    
  Intensity_Moment_Data phi_old(cell_data,angular_quadrature,fem_quadrature,i_old);
  Intensity_Moment_Data phi_new(phi_old);
  phi_new.clear_angle_integrated_intensity();
    
    
  try{
    const double dt = time_data.get_dt_max();
    const int stage = 0;
    std::vector<double> rk_a(1,0.);
    rk_a[0] = time_data.get_a(stage,stage);
    
    transport_sweep.set_time_data(dt , time_data.get_t_start() + dt , rk_a , stage);
  
    /// sweep and update phi only
    bool is_krylov = false;
    bool is_get_k_i = false;
    transport_sweep.set_sweep_type(is_krylov,is_get_k_i);
    /// phi_old is the exact phi, therefore there will be only one richardson iteration to reach convergence!!!
    transport_sweep.sweep_mesh(phi_old,phi_new);
    
    Err_Phi err_phi;
    phi_new.normalized_difference(phi_old,err_phi);
    std::cout << "Maximum normalized difference: " << err_phi.get_worst_err() << std::endl;
        
    if(fabs(err_phi.get_worst_err() ) > tol)
      throw Dark_Arts_Exception(TRANSPORT , "Not getting zero change in phi");
    
    /// perform the same sweep , but check that k_i is zero
    is_krylov = false;
    is_get_k_i = true;
    transport_sweep.set_sweep_type(is_krylov,is_get_k_i);
    transport_sweep.sweep_mesh(phi_old,phi_new);
    
    Eigen::VectorXd k_i_loc = Eigen::VectorXd::Zero(n_p);
    
    const int grp = 0;
    for(int cell = 0 ; cell < n_cell; cell++)
    {
      Eigen::VectorXd phi_new_loc = Eigen::VectorXd::Zero(n_p);
      Eigen::VectorXd phi_old_loc = Eigen::VectorXd::Zero(n_p);
      phi_new.get_cell_angle_integrated_intensity(cell,grp,0,phi_new_loc);
      phi_old.get_cell_angle_integrated_intensity(cell,grp,0,phi_old_loc);
      std::cout << "Phi_new:\n" << phi_new_loc << "\nPhi_old:\n" << phi_old_loc << std::endl;
      for(int dir = 0; dir < n_dir ; dir++)
      {
        ki.get_ki(cell, grp, dir, stage, k_i_loc);        
        std::cout << "K_i for stage: " << stage << " cell: " << cell << " direction: " << dir << " group: " << grp << std::endl << k_i_loc << std::endl; 
        
        for(int el = 0; el < n_p ; el++)
        {
          if( !finite(k_i_loc(el)) )
            throw Dark_Arts_Exception(TRANSPORT , "NAN Ki");
        }
      }        
    }
    
    std::cout << "Advancing intensity now " << std::endl;
    /// now test out that when we add k_i we don't change the solution
    Intensity_Data i_new(i_old);
    ki.advance_intensity(i_old, dt, time_data);
    Eigen::VectorXd i_old_loc = Eigen::VectorXd::Zero(n_p);
    Eigen::VectorXd i_new_loc = Eigen::VectorXd::Zero(n_p);
    for(int cell = 0; cell < n_cell ; cell++)
    {
      for(int dir = 0; dir < n_dir ; dir++)
      {
        i_old.get_cell_intensity(cell,grp,dir,i_old_loc);
        i_new.get_cell_intensity(cell,grp,dir,i_new_loc);
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
