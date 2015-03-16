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
#include "Local_MIP_Assembler.h"

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
  
  /// Create a Materials object that contains all opacity, heat capacity, and source objects
  Materials materials( input_reader, fem_quadrature , cell_data, angular_quadrature );  
  Time_Data time_data(input_reader);
  Temperature_Data t_old(fem_quadrature, input_reader, cell_data);  

  
  
  try{
    /// test matrix integrations
    Local_MIP_Assembler local_assembler(fem_quadrature,input_reader);
    
    /// We're going to impose an arbitrary situation now.  Then feed the same situation to the ``working'' MIP DSA I have in MATLAB
    /// truth will be the MATLAB DSA's answer
    
    const double d_cm1_r = 1.3;
    const double d_c_l = 2.1;
    const double d_c_r = 3.3;
    const double d_cp1_l = 5;
    
    const double dx_cm1 = 1.5;
    const double dx_c = 4.1;
    const double dx_cp1 = 2.2;
    
    const double kappa_m12 = 0.3;
    const double kappa_p12 = 0.7;
    
    Eigen::MatrixXd r_sig_a = Eigen::MatrixXd::Zero(3,3);
    r_sig_a(0,0) = 9.3;
    r_sig_a(1,1) = 2.0;
    r_sig_a(2,2) = 1.05;
    
    
    
    Eigen::MatrixXd s_mat = Eigen::MatrixXd::Zero(3,3);
    
    s_mat(0,0) = 2.3;
    s_mat(1,1) = 1.3;
    s_mat(2,2) = 2.9;
    
    s_mat(0,1) = -2.7;
    s_mat(1,0) = -2.7;
    
    s_mat(0,2) = -2.6;
    s_mat(2,0) = -2.6; 
    
    s_mat(1,2) = -3.1;
    s_mat(2,1) = -3.1;
    
    Eigen::MatrixXd cell_cm1 = Eigen::MatrixXd::Zero(3,3);
    Eigen::MatrixXd cell_c = Eigen::MatrixXd::Zero(3,3);
    Eigen::MatrixXd cell_cp1 = Eigen::MatrixXd::Zero(3,3);
    
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(3);   
    /// left boundary matrices
    Eigen::MatrixXd expected_cell_c = Eigen::MatrixXd::Zero(3,3);
    Eigen::MatrixXd expected_cell_cp1 = Eigen::MatrixXd::Zero(3,3);
    
    expected_cell_c(0,0) = 13.436585366 ; 
    expected_cell_c(0,1) = -3.724390244 ; 
    expected_cell_c(0,2) = -2.746341463 ; 
    expected_cell_c(1,0) = -3.724390244 ; 
    expected_cell_c(1,1) =  3.300000000 ; 
    expected_cell_c(1,2) = -1.490243902 ; 
    expected_cell_c(2,0) = -2.746341463 ; 
    expected_cell_c(2,1) = -1.490243902 ; 
    expected_cell_c(2,2) =  2.235365854 ; 
    expected_cell_cp1(0,0) =  0.402439024 ; 
    expected_cell_cp1(0,1) =  0.000000000 ; 
    expected_cell_cp1(0,2) =  0.000000000 ; 
    expected_cell_cp1(1,0) = -1.609756098 ; 
    expected_cell_cp1(1,1) =  0.000000000 ; 
    expected_cell_cp1(1,2) =  0.000000000 ; 
    expected_cell_cp1(2,0) =  3.916407982 ; 
    expected_cell_cp1(2,1) = -4.545454545 ; 
    expected_cell_cp1(2,2) =  1.136363636 ; 
    
    /// right boundary matrices
    Eigen::MatrixXd middle_expected_cell_cm1 = Eigen::MatrixXd::Zero(3,3);
    Eigen::MatrixXd middle_expected_cell_c = Eigen::MatrixXd::Zero(3,3);
    Eigen::MatrixXd middle_expected_cell_cp1 = Eigen::MatrixXd::Zero(3,3);
    
    local_assembler.calculate_left_boundary_matrices(kappa_m12, kappa_p12 ,   
      dx_c, dx_cp1, d_c_l , d_c_r , d_cp1_l, r_sig_a, s_mat,
      cell_c, cell_cp1);
      
      
    std::cout << "\nLeft boundary\nExpected cell c: \n" << expected_cell_c <<
      "\n Calculated cell_c:\n" << cell_c << "\n Expected cell cp1:\n" << expected_cell_cp1 <<
      "\nCalculated cell cp1: \n" << cell_cp1 <<std::endl;
    for(int i=0 ; i < 3 ; i++)
    {
      for(int j=0 ; j < 3 ; j++)
      {
        if(fabs( expected_cell_c(i,j) - cell_c(i,j) ) > tol)
          throw Dark_Arts_Exception(MIP, "Not calculating left boundary c matrix correctly");
      
        if(fabs( expected_cell_cp1(i,j) - cell_cp1(i,j) ) > tol)
          throw Dark_Arts_Exception(MIP, "Not calculating left boundary cp1 matrix correctly");
      }
    }    
    cell_cp1 = Eigen::MatrixXd::Zero(3,3);
    cell_c = Eigen::MatrixXd::Zero(3,3);
    
    /// interior matrices    ****************************
    local_assembler.calculate_interior_matrices(kappa_m12, kappa_p12 ,   
      dx_cm1, dx_c, dx_cp1, d_cm1_r , d_c_l , d_c_r , d_cp1_l, r_sig_a, s_mat,
      cell_cm1,cell_c, cell_cp1);
      
    middle_expected_cell_cm1(0,0) =  0.433333333 ; 
    middle_expected_cell_cm1(0,1) = -1.733333333 ; 
    middle_expected_cell_cm1(0,2) =  1.768292683 ; 
    middle_expected_cell_cm1(1,0) =  0.000000000 ; 
    middle_expected_cell_cm1(1,1) =  0.000000000 ; 
    middle_expected_cell_cm1(1,2) = -1.024390244 ; 
    middle_expected_cell_cm1(2,0) =  0.000000000 ; 
    middle_expected_cell_cm1(2,1) =  0.000000000 ; 
    middle_expected_cell_cm1(2,2) =  0.256097561 ; 
    middle_expected_cell_c(0,0) = 10.363414634 ; 
    middle_expected_cell_c(0,1) = -1.675609756 ; 
    middle_expected_cell_c(0,2) = -3.258536585 ; 
    middle_expected_cell_c(1,0) = -1.675609756 ; 
    middle_expected_cell_c(1,1) =  3.300000000 ; 
    middle_expected_cell_c(1,2) = -1.490243902 ; 
    middle_expected_cell_c(2,0) = -3.258536585 ; 
    middle_expected_cell_c(2,1) = -1.490243902 ; 
    middle_expected_cell_c(2,2) =  2.235365854 ; 
    middle_expected_cell_cp1(0,0) =  0.402439024 ; 
    middle_expected_cell_cp1(0,1) =  0.000000000 ; 
    middle_expected_cell_cp1(0,2) =  0.000000000 ; 
    middle_expected_cell_cp1(1,0) = -1.609756098 ; 
    middle_expected_cell_cp1(1,1) =  0.000000000 ; 
    middle_expected_cell_cp1(1,2) =  0.000000000 ; 
    middle_expected_cell_cp1(2,0) =  3.916407982 ; 
    middle_expected_cell_cp1(2,1) = -4.545454545 ; 
    middle_expected_cell_cp1(2,2) =  1.136363636 ; 
        
    std::cout << "\nInterior\nExpected cell cm1: \n" << middle_expected_cell_cm1 <<
      "\n Calculated cell_cm1:\n" << cell_cm1 <<" \nExpected cell c: \n" << middle_expected_cell_c <<
      "\n Calculated cell_c:\n" << cell_c << "\n Expected cell cp1:\n" << middle_expected_cell_cp1 <<
      "\nCalculated cell cp1: \n" << cell_cp1 <<std::endl;
      
    for(int i=0 ; i < 3 ; i++)
    {
      for(int j=0 ; j < 3 ; j++)
      {        
        if(fabs( middle_expected_cell_cm1(i,j) - cell_cm1(i,j) ) > tol)
          throw Dark_Arts_Exception(MIP, "Not calculating interior cm1 matrix correctly");
          
        if(fabs( middle_expected_cell_c(i,j) - cell_c(i,j) ) > tol)
          throw Dark_Arts_Exception(MIP, "Not calculating interior c matrix correctly");
      
        if(fabs( middle_expected_cell_cp1(i,j) - cell_cp1(i,j) ) > tol)
          throw Dark_Arts_Exception(MIP, "Not calculating interior cp1 matrix correctly");
      }
    }    
    cell_cp1 = Eigen::MatrixXd::Zero(3,3);
    cell_c = Eigen::MatrixXd::Zero(3,3);
    cell_cm1 = Eigen::MatrixXd::Zero(3,3);
    
    /// Right boundary matrices    ******************************
    Eigen::MatrixXd right_expected_cell_cm1 = Eigen::MatrixXd::Zero(3,3);
    Eigen::MatrixXd right_expected_cell_c = Eigen::MatrixXd::Zero(3,3);
    local_assembler.calculate_right_boundary_matrices(kappa_m12, kappa_p12 ,   
      dx_cm1, dx_c, d_cm1_r , d_c_l , d_c_r , r_sig_a, s_mat,
      cell_cm1,cell_c);
      
      right_expected_cell_c(0,0) = 10.363414634 ; 
      right_expected_cell_c(0,1) = -1.675609756 ; 
      right_expected_cell_c(0,2) = -3.258536585 ; 
      right_expected_cell_c(1,0) = -1.675609756 ; 
      right_expected_cell_c(1,1) =  3.300000000 ; 
      right_expected_cell_c(1,2) = -1.490243902 ; 
      right_expected_cell_c(2,0) = -3.258536585 ; 
      right_expected_cell_c(2,1) = -1.490243902 ; 
      right_expected_cell_c(2,2) =  2.235365854 ; 
      right_expected_cell_cm1(0,0) =  0.433333333 ; 
      right_expected_cell_cm1(0,1) = -1.733333333 ; 
      right_expected_cell_cm1(0,2) =  1.768292683 ; 
      right_expected_cell_cm1(1,0) =  0.000000000 ; 
      right_expected_cell_cm1(1,1) =  0.000000000 ; 
      right_expected_cell_cm1(1,2) = -1.024390244 ; 
      right_expected_cell_cm1(2,0) =  0.000000000 ; 
      right_expected_cell_cm1(2,1) =  0.000000000 ; 
      right_expected_cell_cm1(2,2) =  0.256097561 ; 
      
    std::cout << "\nRight boundary\nExpected cell cm1: \n" << right_expected_cell_cm1 <<
      "\n Calculated cell_cm1:\n" << cell_cm1 <<"\n Expected cell c: \n" << right_expected_cell_c <<
      "\n Calculated cell_c:\n" << cell_c << std::endl;
      
    for(int i=0 ; i < 3 ; i++)
    {
      for(int j=0 ; j < 3 ; j++)
      {
        if(fabs( right_expected_cell_c(i,j) - cell_c(i,j) ) > tol)
          throw Dark_Arts_Exception(MIP, "Not calculating right boundary c matrix correctly");
      
        if(fabs( right_expected_cell_cm1(i,j) - cell_cm1(i,j) ) > tol)
          throw Dark_Arts_Exception(MIP, "Not calculating right boundary cm1 matrix correctly");
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
