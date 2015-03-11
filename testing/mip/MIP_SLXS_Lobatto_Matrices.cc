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
  Materials materials( input_reader, fem_quadrature , cell_data, angular_quadrature );  
  Time_Data time_data(input_reader);
  Temperature_Data t_old(fem_quadrature, input_reader, cell_data);  

  Eigen::MatrixXd r_sig_a_reg1 = Eigen::MatrixXd::Zero(n_dfem_p,n_dfem_p);
  Eigen::MatrixXd r_sig_a_reg2 = Eigen::MatrixXd::Zero(n_dfem_p,n_dfem_p);
  
  Eigen::MatrixXd r_sig_s_reg1 = Eigen::MatrixXd::Zero(n_dfem_p,n_dfem_p);
  Eigen::MatrixXd r_sig_s_reg2 = Eigen::MatrixXd::Zero(n_dfem_p,n_dfem_p);
  
  Eigen::MatrixXd dimless_mass = Eigen::MatrixXd::Zero(n_dfem_p,n_dfem_p);
  Eigen::MatrixXd r_cv_1 = Eigen::MatrixXd::Zero(n_dfem_p,n_dfem_p);
  Eigen::MatrixXd r_cv_2 = Eigen::MatrixXd::Zero(n_dfem_p,n_dfem_p);
 
  Eigen::MatrixXd r_sig_tau_1 = Eigen::MatrixXd::Zero(n_dfem_p,n_dfem_p);
  Eigen::MatrixXd r_sig_tau_2 = Eigen::MatrixXd::Zero(n_dfem_p,n_dfem_p);
  
  Eigen::MatrixXd d_1 = Eigen::MatrixXd::Zero(n_dfem_p,n_dfem_p);
  Eigen::MatrixXd d_2 = Eigen::MatrixXd::Zero(n_dfem_p,n_dfem_p);
  
  const double sig_a_1 = 0.7;
  const double sig_a_2 = 0.9;
  const double sig_s_1 = 1.0;
  const double sig_s_2 = 0.2;
  const double cv_1 = 1.0;
  const double cv_2 = 1.2;
  const double dx_1 = 4./3.;
  const double dx_2 = 1./3.;
  
  const double w_1 = 1.;
  const double w_2 = 1.;
  
  const double dt = time_data.get_dt_min();  
  const double c = 1.0;
  const double rk_a = 1.0;
    
  const double temp_1 = 0.5;
  const double temp_2 = 0.4;
  
  dimless_mass(0,0) = w_1;
  dimless_mass(1,1) = w_2;
  
  r_sig_a_reg1 = dx_1/2.*sig_a_1*dimless_mass;
  r_sig_a_reg2 = dx_2/2.*sig_a_2*dimless_mass;  
  
  r_sig_s_reg1 = dx_1/2.*sig_s_1*dimless_mass;
  r_sig_s_reg2 = dx_2/2.*sig_s_2*dimless_mass;
  
  r_cv_1 = dx_1/2.*cv_1*dimless_mass;
  r_cv_2 = dx_2/2.*cv_2*dimless_mass;
  
  std::cout << "r_sig_a_1: \n" << r_sig_a_reg1 << std::endl;
  std::cout << "r_sig_a_2: \n" << r_sig_a_reg2 << std::endl;
  std::cout << "r_sig_s_1: \n" << r_sig_s_reg1 << std::endl;
  std::cout << "r_sig_s_2: \n" << r_sig_s_reg2 << std::endl;
  std::cout << "r_cv_1: \n" << r_cv_1 << std::endl;
  std::cout << "r_cv_2: \n" << r_cv_2 << std::endl;
  
  r_sig_tau_1 = r_sig_a_reg1+r_sig_s_reg1+(1./(c*dt*1.0))*dx_1/2.*dimless_mass;
  r_sig_tau_2 = r_sig_a_reg2+r_sig_s_reg2+(1./(c*dt*1.0))*dx_2/2.*dimless_mass;
  
  d_1(0,0) = 4.*pow(temp_1,3)/sn_w; d_1(1,1) = 4.*pow(temp_1,3)/sn_w;
  d_2(0,0) = 4.*pow(temp_2,3)/sn_w; d_2(1,1) = 4.*pow(temp_2,3)/sn_w;
 
  Eigen::MatrixXd r_pseudo_s_1 = Eigen::MatrixXd::Zero(n_dfem_p,n_dfem_p);
  Eigen::MatrixXd r_pseudo_s_2 = Eigen::MatrixXd::Zero(n_dfem_p,n_dfem_p);
  Eigen::MatrixXd r_pseudo_a_1 = Eigen::MatrixXd::Zero(n_dfem_p,n_dfem_p);
  Eigen::MatrixXd r_pseudo_a_2 = Eigen::MatrixXd::Zero(n_dfem_p,n_dfem_p);
  
  int i = 0;
  r_pseudo_s_1(i,i) = r_sig_a_reg1(i,i)*d_1(i,i)/( 1.+sn_w*dt*rk_a/r_cv_1(i,i)*r_sig_a_reg1(i,i)*d_1(i,i) )/r_cv_1(i,i)*r_sig_a_reg1(i,i) ;
  i = 1;
  r_pseudo_s_1(i,i) = r_sig_a_reg1(i,i)*d_1(i,i)/( 1.+sn_w*dt*rk_a/r_cv_1(i,i)*r_sig_a_reg1(i,i)*d_1(i,i) )/r_cv_1(i,i)*r_sig_a_reg1(i,i) ;
  r_pseudo_s_1 *= sn_w*dt*rk_a;
  r_pseudo_s_1 += r_sig_s_reg1;
  
  i=0;
  r_pseudo_s_2(i,i) = r_sig_a_reg2(i,i)*d_2(i,i)/( 1.+sn_w*dt*rk_a/r_cv_2(i,i)*r_sig_a_reg2(i,i)*d_2(i,i) )/r_cv_2(i,i)*r_sig_a_reg2(i,i) ;
  i = 1;
  r_pseudo_s_2(i,i) = r_sig_a_reg2(i,i)*d_2(i,i)/( 1.+sn_w*dt*rk_a/r_cv_2(i,i)*r_sig_a_reg2(i,i)*d_2(i,i) )/r_cv_2(i,i)*r_sig_a_reg2(i,i) ;
  r_pseudo_s_2 *= sn_w*dt*rk_a;
  r_pseudo_s_2 += r_sig_s_reg2;
  
  r_pseudo_a_1 = r_sig_tau_1 - r_pseudo_s_1 ;
  r_pseudo_a_2 = r_sig_tau_2 - r_pseudo_s_2 ;
  
  try{
    std::shared_ptr<V_Diffusion_Matrix_Creator> matrix_creator;
    matrix_creator = std::make_shared<Diffusion_Matrix_Creator_Grey>(fem_quadrature,materials,angular_quadrature,t_old);
    
    const double t_stage = time_data.get_t_start();
    
    Eigen::MatrixXd calc_r_sig_a_1 = Eigen::MatrixXd::Zero(n_dfem_p,n_dfem_p);
    Eigen::MatrixXd calc_r_sig_a_2 = calc_r_sig_a_1;
    Eigen::MatrixXd calc_r_sig_s_1 = calc_r_sig_a_1;
    Eigen::MatrixXd calc_r_sig_s_2 = calc_r_sig_a_1;
    
    matrix_creator->set_time_data(dt,t_stage,rk_a);
    matrix_creator->set_cell_group_information(1,0);
    matrix_creator->calculate_pseudo_r_sig_a_and_r_sig_s(calc_r_sig_a_1,calc_r_sig_s_1);
    
    matrix_creator->set_cell_group_information(4,0);
    matrix_creator->calculate_pseudo_r_sig_a_and_r_sig_s(calc_r_sig_a_2,calc_r_sig_s_2);
    
    std::cout << "Region 0\n"; 
    std::cout << "Calculated Pseudo r_sig_a: \n" << calc_r_sig_a_1 << "\n Expected: \n" << r_pseudo_a_1 << std::endl;
    std::cout << "Calculated Pseudo r_sig_s: \n" << calc_r_sig_s_1 << "\n Expected: \n" << r_pseudo_s_1 << std::endl;
    
    for(int i=0; i < n_dfem_p ; i++)
    {
      for(int j=0; j < n_dfem_p ; j++)
      {
        if( fabs(calc_r_sig_a_1(i,j) - r_pseudo_a_1(i,j) ) > tol)
          throw Dark_Arts_Exception(MIP, "Not calculating region 0 grey pseudo r_sig_a correctly");
          
        if( fabs(calc_r_sig_s_1(i,j) - r_pseudo_s_1(i,j) ) > tol)
          throw Dark_Arts_Exception(MIP, "Not calculating region 0 grey pseudo r_sig_s correctly");
          
        if( fabs(calc_r_sig_a_2(i,j) - r_pseudo_a_2(i,j) ) > tol)
          throw Dark_Arts_Exception(MIP, "Not calculating region 1 grey pseudo r_sig_a correctly");
          
        if( fabs(calc_r_sig_s_2(i,j) - r_pseudo_s_2(i,j) ) > tol)
          throw Dark_Arts_Exception(MIP, "Not calculating region 1 grey pseudo r_sig_s correctly");
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
