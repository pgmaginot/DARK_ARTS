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
#include "Diffusion_Matrix_Creator_Grey.h"
#include "Sweep_Matrix_Creator_Grey.h"

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
  
  /// Initialize a Quadrule object to be able to get all of the quadrature we need
  Quadrule_New quad_fun;  
  Fem_Quadrature fem_quadrature( input_reader , quad_fun);  
  Cell_Data cell_data( input_reader );  
  Angular_Quadrature angular_quadrature( input_reader , quad_fun );    
  
  
  const int n_cell = cell_data.get_total_number_of_cells();
  const int n_dir = angular_quadrature.get_number_of_dir();
  const int n_dfem_p = fem_quadrature.get_number_of_interpolation_points();
  const double sn_w = angular_quadrature.get_sum_w();
  
  /// Create a Materials object that contains all opacity, heat capacity, and source objects


  const int n_l_mom = angular_quadrature.get_number_of_leg_moments();
  Materials materials( input_reader, fem_quadrature , cell_data, angular_quadrature );  
    Time_Data time_data(input_reader);
    Temperature_Data t_old(fem_quadrature, input_reader, cell_data);  
    Temperature_Data t_star(fem_quadrature, input_reader, cell_data);  
    Intensity_Data i_old(cell_data, angular_quadrature, fem_quadrature, materials,input_reader);
    K_Temperature kt(n_cell , time_data.get_number_of_stages(), fem_quadrature);
    K_Intensity ki(n_cell ,  time_data.get_number_of_stages(),fem_quadrature, angular_quadrature);
    t_star = t_old;
      const int n_stages = time_data.get_number_of_stages();
    std::shared_ptr<V_Sweep_Matrix_Creator> sweep_matrix;
    sweep_matrix = std::make_shared<Sweep_Matrix_Creator_Grey>(fem_quadrature, materials, n_stages, angular_quadrature.get_sum_w() ,
      n_l_mom, t_old, i_old, kt, ki,t_star) ;
      
    Diffusion_Matrix_Creator_Grey diffusion_matrix(fem_quadrature, materials, angular_quadrature, t_old, n_cell, input_reader);
    
    Eigen::MatrixXd sweep_sig_t = Eigen::MatrixXd(n_dfem_p,n_dfem_p);
    Eigen::MatrixXd sweep_sig_s = Eigen::MatrixXd(n_dfem_p,n_dfem_p);
    Eigen::MatrixXd sweep_sig_a;
    
    Eigen::MatrixXd dsa_sig_s = Eigen::MatrixXd::Zero(n_dfem_p,n_dfem_p);
    Eigen::MatrixXd dsa_sig_a = Eigen::MatrixXd::Zero(n_dfem_p,n_dfem_p);
    
    const double dt = time_data.get_dt_max();
    const double time_stage = dt + time_data.get_t_start();
    std::vector<double> rk_a_of_stage_i(n_stages,0.);
    rk_a_of_stage_i[0] = 1.;

  try{
    sweep_matrix->set_time_data(dt,time_stage,rk_a_of_stage_i,0);
    diffusion_matrix.set_time_data( dt, time_stage, 1. );

    for(int cell = 0 ; cell < n_cell ; cell++)
    {      
      sweep_matrix->update_cell_dependencies(cell); 
      sweep_matrix->update_group_dependencies(0);      
      
      std::cout << "Dies here"<<std::endl;
      sweep_matrix->get_r_sig_t(sweep_sig_t); 
      sweep_matrix->get_r_sig_s(sweep_sig_s,0);
      
      sweep_sig_a = sweep_sig_t - sweep_sig_s;
      
      std::cout << "Here"<<std::endl;
      
      diffusion_matrix.set_cell_group_information( cell, 0, cell_data.get_cell_width(cell) );    
      diffusion_matrix.calculate_pseudo_r_sig_a_and_pseudo_r_sig_s(dsa_sig_a,dsa_sig_s);
      
      std::cout << "Cell: " << cell << std::endl;
      std::cout << "Transport sweep matrices. \n Sig_s:\n" << sweep_sig_s << "\n Sig_a: \n" << sweep_sig_a << std::endl;
      
      std::cout << "Diffusion matrices. \n Sig_s:\n" << dsa_sig_s << "\n Sig_a: \n" << dsa_sig_a << std::endl;
      
      
      for(int i = 0;i < n_dfem_p ; i++)
      {
        for(int j=0 ; j < n_dfem_p ; j++)
        {
          if( fabs(sweep_sig_a(i,j) - dsa_sig_a(i,j) ) > tol)
            throw Dark_Arts_Exception(MIP , "Difference in R_sig_a");
            
          if( fabs(sweep_sig_s(i,j) - dsa_sig_s(i,j) ) > tol)
            throw Dark_Arts_Exception(MIP , "Difference in R_sig_s");
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
