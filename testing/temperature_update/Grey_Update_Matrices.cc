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
#include "Output_Generator.h"
#include "Err_Temperature.h"
#include "L2_Error_Calculator.h"
#include "Temperature_Update.h"
#include "Grey_Temperature_Matrix_Creator.h"
#include <memory>
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
  const int n_cell = cell_data.get_total_number_of_cells();
  const int n_p = fem_quadrature.get_number_of_interpolation_points();  
  
  /// Create a Materials object that contains all opacity, heat capacity, and source objects
  Materials materials( input_reader, fem_quadrature , cell_data, angular_quadrature );  
  Time_Data time_data(input_reader);
  Temperature_Data t_old(fem_quadrature, input_reader, cell_data);  
  Intensity_Data i_old(cell_data, angular_quadrature, fem_quadrature, materials,input_reader);
    
  const int n_stages = time_data.get_number_of_stages();
  
  K_Temperature kt(n_cell , n_stages, fem_quadrature);
  K_Intensity ki(n_cell ,  n_stages, fem_quadrature, angular_quadrature);  
    
  Intensity_Moment_Data ard_phi(cell_data,angular_quadrature,fem_quadrature,i_old);  

  
  const double tol = 1.0E-8;
  try{
    if(n_p != 4)
      throw Dark_Arts_Exception(FEM, "Expecting Cubic self lumping lobatto for this test");
  
    if(fem_quadrature.get_integration_type() != SELF_LUMPING)
      throw Dark_Arts_Exception(FEM, "Expecting Cubic self lumping lobatto for this test");
      
    std::shared_ptr<V_Temperature_Matrix_Creator> temperature_matrix;
    temperature_matrix = std::make_shared<Grey_Temperature_Matrix_Creator>
      (fem_quadrature, cell_data, materials, angular_quadrature, n_stages, t_old, ard_phi);
    
    Temperature_Data t_star( cell_data.get_total_number_of_cells(), fem_quadrature);
         
    const int stage = 0;

    double dt = time_data.get_dt_max();
    
    double time_stage = time_data.get_t_start() + dt*time_data.get_c(stage);
    
    std::vector<double> rk_a_of_stage_i(n_stages,0.);
    
    for(int i = 0; i <= stage; i++)
      rk_a_of_stage_i[i] = time_data.get_a(stage,i);
      
    temperature_matrix->set_time_data(rk_a_of_stage_i, dt, time_stage, stage);
      
    const int cell = 1;
    std::vector<double> dfem_interp_points;
    fem_quadrature.get_dfem_interpolation_point(dfem_interp_points);
    
    std::vector<double> lobatto(n_p,0.);
    lobatto[0] = -1.;
    lobatto[1] = -0.4472135955;
    lobatto[2] = 0.4472135955;
    lobatto[3] = 1.;
    
    std::vector<double> lobatto_w(n_p,0.);
    lobatto_w[0] = 0.1666666667;
    lobatto_w[1] = 0.8333333333;
    lobatto_w[2] = 0.8333333333;
    lobatto_w[3] = 0.1666666667;
    
    for(int el = 0 ; el < n_p ; el++)
    {
      if( fabs(dfem_interp_points[el] - lobatto[el]) > tol)
        throw Dark_Arts_Exception(FEM , "Expecting 4 point lobatto quadrature for this test");
    }    
    
    t_star.mms_cheat(time_stage, cell_data, dfem_interp_points, input_reader);
    ard_phi.mms_cheat(time_stage,cell_data, dfem_interp_points , input_reader,angular_quadrature);
    
    Eigen::VectorXd t_star_vec = Eigen::VectorXd::Zero(n_p);
    Eigen::VectorXd phi_vec = Eigen::VectorXd::Zero(n_p);
    
    ard_phi.get_cell_angle_integrated_intensity(cell,0,0,phi_vec);
    t_star.get_cell_temperature(cell,t_star_vec);
    const double dx = cell_data.get_cell_width(cell);
    
    Eigen::MatrixXd r_cv = Eigen::MatrixXd::Zero(n_p,n_p);
    Eigen::MatrixXd r_sig_a = Eigen::MatrixXd::Zero(n_p,n_p);
    Eigen::MatrixXd d_matrix = Eigen::MatrixXd::Zero(n_p,n_p);
    Eigen::VectorXd s_t = Eigen::VectorXd::Zero(n_p);
    Eigen::VectorXd planck_vec = Eigen::VectorXd::Zero(n_p);
    
    for(int el = 0 ; el < n_p ; el++)
    {
      planck_vec(el) = pow(t_star_vec(el),4)/2.;
      r_cv(el,el) = dx/2.*1.2*lobatto_w[el];
      r_sig_a(el,el) = dx/2.*1.0*lobatto_w[el];
      d_matrix(el,el) = 4.*pow(t_star_vec(el),3)/2.;
    }
    s_t(0) =    1.3787471617e-002 ; 
    s_t(1) =    6.7388714824e-002 ; 
    s_t(2) =    6.6345941281e-002 ; 
    s_t(3) =    1.3329249497e-002 ; 
    
    Eigen::MatrixXd i_mat = Eigen::MatrixXd::Identity(n_p,n_p);
    Eigen::MatrixXd expected_coeff = Eigen::MatrixXd::Zero(n_p,n_p);
    Eigen::MatrixXd scratch = Eigen::MatrixXd::Zero(n_p,n_p);
    Eigen::MatrixXd scratch2 = Eigen::MatrixXd::Zero(n_p,n_p);
    
    for(int el = 0 ; el < n_p ; el++)
    {
      scratch2(el,el) = r_sig_a(el,el)*d_matrix(el,el)/r_cv(el,el);
    }
    
    scratch = r_cv.fullPivLu().solve( r_sig_a*d_matrix);
    
    std::cout << "Actual matrix:\n" << scratch2 << "\n fullPivLu solution: \n" << scratch << std::endl;
    for(int i = 0 ; i < n_p ; i++)
    {
      for(int j = 0 ; j < n_p ; j++)
      {
        if(fabs(scratch2(i,j) - scratch(i,j)) > tol)
          throw Dark_Arts_Exception(FEM,"Not getting inverse as expected");
      }
    }
    
    expected_coeff = i_mat + 2.*dt*1.*scratch;        
    Eigen::VectorXd expected_rhs = Eigen::VectorXd::Zero(n_p);
    Eigen::VectorXd scratch_vec1 = Eigen::VectorXd::Zero(n_p);
    Eigen::VectorXd scratch_vec2 = Eigen::VectorXd::Zero(n_p);
    scratch_vec1 = phi_vec - 2.*planck_vec;
    scratch_vec2 = r_sig_a*scratch_vec1;    
    scratch_vec2 += s_t;
    
    expected_rhs = r_cv.fullPivLu().solve(scratch_vec2);
    expected_rhs *= dt*1.;
    Eigen::VectorXd t_old_vec = Eigen::VectorXd::Zero(n_p);
    t_old.get_cell_temperature(cell,t_old_vec);
    expected_rhs += t_old_vec - t_star_vec;
    
    Eigen::MatrixXd calc_coefficient = Eigen::MatrixXd::Zero(n_p,n_p);
    Eigen::VectorXd calc_rhs = Eigen::VectorXd::Zero(n_p);
    temperature_matrix->calculate_update_quantities(cell, t_star_vec, kt, calc_coefficient , calc_rhs);
    
    std::cout << "Expected coefficient matrix:\n" << expected_coeff << "\n Calculatated coefficient matrix:\n" << calc_coefficient << std::endl;
    for(int i = 0 ; i < n_p ; i++)
    {
      for(int j = 0 ; j < n_p ; j++)
      {
        if(fabs(expected_coeff(i,j) - calc_coefficient(i,j)) > tol)
          throw Dark_Arts_Exception(FEM,"Not expected coefficient matrix");
      }
    }
    
    std::cout << "Expected rhs: \n" << expected_rhs << " \n Calculated_rhs:\n" << calc_rhs << std::endl;
    for(int i=0 ; i < n_p ; i++)
    {
      if(fabs(expected_rhs(i) - calc_rhs(i))>tol)
        throw Dark_Arts_Exception(TIME_MARCHER , "Difference in rhs");
    }
    
    Eigen::VectorXd delta_vec = Eigen::VectorXd::Zero(n_p);
    delta_vec = calc_coefficient.fullPivLu().solve(calc_rhs);
    
    std::cout << "Delta_vec\n" << delta_vec << std::endl;
    
    Eigen::VectorXd calculated_kt = Eigen::VectorXd::Zero(n_p);
    temperature_matrix->calculate_k_t(cell, t_star_vec, calculated_kt);
    
    Eigen::VectorXd expected_kt = Eigen::VectorXd::Zero(n_p);
    
    std::cout << "R_cv after inversion\n" << r_cv << std::endl;
    scratch_vec1 = Eigen::VectorXd::Zero(n_p);
    scratch_vec1 = r_sig_a*(phi_vec - 2.*planck_vec) + s_t;
    std::cout << "Expected rhs before hitting with r_cv_inv:\n" << scratch_vec1 << std::endl;
    expected_kt = r_cv.fullPivLu().solve( scratch_vec1);
    
    std::cout << "Expected kt:\n" << expected_kt << " \n Calculated kt:\n" << calculated_kt << std::endl;
    for(int el = 0 ; el < n_p ; el++)
    {
      if(fabs(expected_kt(el) - calculated_kt(el)) > tol)
        throw Dark_Arts_Exception(TIME_MARCHER, "difference in kt");
    }
    
    
  }
  catch(const Dark_Arts_Exception& da)
  {
    val = -1;
    da.testing_message();
  }
  
  

  
  return val;
}
