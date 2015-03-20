static char help[] = "A character array that PETSc might be expecting";
#include <petscksp.h>

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
#include "Temperature_Data.h"

#include <Eigen/Dense>
#include "Diffusion_Operator.h"

#include "Dark_Arts_Exception.h"

/**
  Goal of this unit test is to check source moment formation, and reaction matrix for SLXS Lobatto scheme
  -This is a MMS problem with a spatially varying sigma_a, constant cv , zero sig_s
*/ 

int main(int argc, char** argv)
{  
  /// change to -1 if you want ctest to dump this output.  Otherwise this test passes always
  int val = 0;
  PetscErrorCode ierr;  
  PetscMPIInt size;
  
  PetscInitialize(&argc,&argv,(char*)0,help);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);
  CHKERRQ(ierr);
  if (size != 1) 
    SETERRQ(PETSC_COMM_WORLD,1,"DARK_ARTS is written to be serial only!");
      
  
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
  
  /// Initialize a Quadrule object to be able to get all of the quadrature we need
  Quadrule_New quad_fun;  
  Fem_Quadrature fem_quadrature( input_reader , quad_fun);  
  Cell_Data cell_data( input_reader );  
  Angular_Quadrature angular_quadrature( input_reader , quad_fun );    
    
  /// Create a Materials object that contains all opacity, heat capacity, and source objects
    Materials materials( input_reader, fem_quadrature , cell_data, angular_quadrature );  
    Time_Data time_data(input_reader);
    Temperature_Data t_old(fem_quadrature, input_reader, cell_data);      

    const double dt = time_data.get_dt_max();
    const double time_stage = dt + time_data.get_t_start();
    double rk_a_ii = 1.;

    const int n_p = fem_quadrature.get_number_of_interpolation_points() ;
  try{
    Diffusion_Operator diffusion_operator(input_reader, fem_quadrature, cell_data, 
      materials, angular_quadrature, 1,t_old, true,1.0E-20, 1.0E-10 , 2.);
      
    diffusion_operator.set_time_data(dt, time_stage, rk_a_ii);
    
    diffusion_operator.dump_matrix();
        
    std::vector<double> ref_norm(1,0.);
    Intensity_Moment_Data phi_old(cell_data, angular_quadrature, fem_quadrature, ref_norm);  
    Intensity_Moment_Data phi_new(phi_old);
    
    Eigen::VectorXd phi_old_vec = Eigen::VectorXd::Zero(n_p);
    Eigen::VectorXd phi_new_vec = Eigen::VectorXd::Zero(n_p);
    for(int el = 0; el < n_p ;el++)
    {
        phi_old_vec(el) = 0.1;
        phi_new_vec(el) = 0.7;
    }
    
    for(int i = 0 ; i < cell_data.get_total_number_of_cells() ; i++)
    {
      phi_old.set_cell_angle_integrated_intensity(i,0,0,phi_old_vec);
      phi_new.set_cell_angle_integrated_intensity(i,0,0,phi_new_vec);
    }
    
    // diffusion_operator.make_and_dump_rhs(phi_new , phi_old);
    
    std::cout << "Update follows\n" << std::endl;
    diffusion_operator.after_rhs_solve_system_and_dump_solution();
    
    double sig_a = 4.0;
    double sig_s = 1.0;
    double sum_sn_w = 2.;
    double c_speed = 1.;
    double cv = 0.5;
    double temp = 0.1*(1.0);
    double d_planck = 4.*pow(temp,3)/sum_sn_w;
    double pseudo_sig_t = 1./(rk_a_ii*dt*c_speed) + sig_a + sig_s;
    
    double value = sum_sn_w*rk_a_ii*dt*sig_a*d_planck;
    double nu = value/(cv + value);
    double pseudo_sig_s = sig_s + nu*sig_a;
    
    std::cout << "Pseudo sig_s: " << pseudo_sig_s << std::endl;
    std::cout << "Pseudp sig_t: " << pseudo_sig_t << std::endl;
      
      
  }
  catch(const Dark_Arts_Exception& da)
  {
    val = -1;
    da.testing_message();
  }
  ierr = PetscFinalize();
  return val;
}
