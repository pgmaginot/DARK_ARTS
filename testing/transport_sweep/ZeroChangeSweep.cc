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

#include "Dark_Arts_Exception.h"

/**
  Goal of this unit test is to check source moment formation, and reaction matrix for SLXS Lobatto scheme
  -This is a MMS problem with a spatially varying sigma_a, constant cv , zero sig_s
*/ 

int main(int argc, char** argv)
{
  int val = 0;
  const double tol = 1.0E-6;
  
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
  const int n_dfem_p = fem_quadrature.get_number_of_interpolation_points();
  const double sn_w = angular_quadrature.get_sum_w();
  
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
    transport_sweep.sweep_mesh(phi_new,phi_old);
    
    Err_Phi err_phi;
    phi_new.normalized_difference(phi_old,err_phi);
    std::cout << "Maximum normalized difference: " << err_phi.get_worst_err() << std::endl;
    
    std::vector<double> phi_old_norm;
    std::vector<double> phi_new_norm;
    phi_new.get_phi_norm(phi_new_norm);
    phi_old.get_phi_norm(phi_old_norm);
    
    std::cout << "Phi_new norm: " << phi_new_norm[0] <<std::endl;
    std::cout << "Phi_old norm: " << phi_old_norm[0] <<std::endl;
    
    if(fabs(err_phi.get_worst_err() ) > tol)
      throw Dark_Arts_Exception(TRANSPORT , "Not getting zero change in phi");
    
    /// perform the same sweep , but check that k_i is zero
    // is_krylov = false;
    // is_get_k_i = true;
    // transport_sweep.set_sweep_type(is_krylov,is_get_k_i);
    // transport_sweep.sweep_mesh(phi_new,phi_old);
    
    
  }
  catch(const Dark_Arts_Exception& da)
  {
    val = -1;
    da.testing_message();
  }
  
  return val;
}
