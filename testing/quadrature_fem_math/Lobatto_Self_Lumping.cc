#include <iostream>
#include <vector>
#include <math.h>
#include <string>
#include <memory>

#include "Cell_Data.h"
#include "Input_Reader.h"
#include "Angular_Quadrature.h"
#include "Quadrule_New.h"
#include "Matrix_Construction_Self_Lumping.h"
#include "Materials.h"
#include "Dark_Arts_Exception.h"

int main(int argc, char** argv)
{
  int val = 0;
  
  /// get all the objects necessary to create an Angular_Quadrature object
  Input_Reader input_reader;  
  try{
    input_reader.read_xml(argv[1]);  
  }
  catch(const Dark_Arts_Exception& da_exception)
  {
    da_exception.testing_message();
    val = -1;
  }
  
  
  
  try{
    Quadrule_New quad_fun;   
    
    Fem_Quadrature fem_quadrature(input_reader, quad_fun);
    Angular_Quadrature angular_quadrature( input_reader , quad_fun );  
    Cell_Data cell_data(input_reader);
    Materials materials( input_reader, fem_quadrature , cell_data, angular_quadrature ); 
    
    const int n_dfem_p = fem_quadrature.get_number_of_interpolation_points();
      
    Eigen::MatrixXd matrix(n_dfem_p, n_dfem_p);
    Eigen::VectorXd vector(n_dfem_p);
        
    std::shared_ptr<V_Matrix_Construction> matrix_construction;
    matrix_construction = std::shared_ptr<V_Matrix_Construction> (new Matrix_Construction_Self_Lumping(fem_quadrature, materials) );     

    std::vector<std::vector<double>> expected_dimensionless_mass( n_dfem_p , std::vector<double> (n_dfem_p , 0.) );   
    std::vector<std::vector<double>> expected_L_mu_positive( n_dfem_p , std::vector<double> (n_dfem_p , 0.) );    
    std::vector<std::vector<double>> expected_L_mu_negative( n_dfem_p , std::vector<double> (n_dfem_p , 0.) ); 
    std::vector<double> expected_f_mu_positive(n_dfem_p,0.);
    std::vector<double> expected_f_mu_negative(n_dfem_p,0.);
    
    const double tol = 1.0E-6;
      
    expected_dimensionless_mass[0][0] = 0.16666667 ; 
    expected_dimensionless_mass[0][1] = 0.00000000 ; 
    expected_dimensionless_mass[0][2] = 0.00000000 ; 
    expected_dimensionless_mass[0][3] = 0.00000000 ; 
    expected_dimensionless_mass[1][0] = 0.00000000 ; 
    expected_dimensionless_mass[1][1] = 0.83333333 ; 
    expected_dimensionless_mass[1][2] = 0.00000000 ; 
    expected_dimensionless_mass[1][3] = 0.00000000 ; 
    expected_dimensionless_mass[2][0] = 0.00000000 ; 
    expected_dimensionless_mass[2][1] = 0.00000000 ; 
    expected_dimensionless_mass[2][2] = 0.83333333 ; 
    expected_dimensionless_mass[2][3] = 0.00000000 ; 
    expected_dimensionless_mass[3][0] = 0.00000000 ; 
    expected_dimensionless_mass[3][1] = 0.00000000 ; 
    expected_dimensionless_mass[3][2] = 0.00000000 ; 
    expected_dimensionless_mass[3][3] = 0.16666667 ; 
    expected_L_mu_positive[0][0] = 0.50000000 ; 
    expected_L_mu_positive[0][1] = 0.67418083 ; 
    expected_L_mu_positive[0][2] = -0.25751416 ; 
    expected_L_mu_positive[0][3] = 0.08333333 ; 
    expected_L_mu_positive[1][0] = -0.67418083 ; 
    expected_L_mu_positive[1][1] = 0.00000000 ; 
    expected_L_mu_positive[1][2] = 0.93169499 ; 
    expected_L_mu_positive[1][3] = -0.25751416 ; 
    expected_L_mu_positive[2][0] = 0.25751416 ; 
    expected_L_mu_positive[2][1] = -0.93169499 ; 
    expected_L_mu_positive[2][2] = -0.00000000 ; 
    expected_L_mu_positive[2][3] = 0.67418083 ; 
    expected_L_mu_positive[3][0] = -0.08333333 ; 
    expected_L_mu_positive[3][1] = 0.25751416 ; 
    expected_L_mu_positive[3][2] = -0.67418083 ; 
    expected_L_mu_positive[3][3] = 0.50000000 ; 
    expected_L_mu_negative[0][0] = -0.50000000 ; 
    expected_L_mu_negative[0][1] = 0.67418083 ; 
    expected_L_mu_negative[0][2] = -0.25751416 ; 
    expected_L_mu_negative[0][3] = 0.08333333 ; 
    expected_L_mu_negative[1][0] = -0.67418083 ; 
    expected_L_mu_negative[1][1] = 0.00000000 ; 
    expected_L_mu_negative[1][2] = 0.93169499 ; 
    expected_L_mu_negative[1][3] = -0.25751416 ; 
    expected_L_mu_negative[2][0] = 0.25751416 ; 
    expected_L_mu_negative[2][1] = -0.93169499 ; 
    expected_L_mu_negative[2][2] = -0.00000000 ; 
    expected_L_mu_negative[2][3] = 0.67418083 ; 
    expected_L_mu_negative[3][0] = -0.08333333 ; 
    expected_L_mu_negative[3][1] = 0.25751416 ; 
    expected_L_mu_negative[3][2] = -0.67418083 ; 
    expected_L_mu_negative[3][3] = -0.50000000 ; 
    expected_f_mu_positive[0] = 1.00000000 ; 
    expected_f_mu_positive[1] = 0.00000000 ; 
    expected_f_mu_positive[2] = -0.00000000 ; 
    expected_f_mu_positive[3] = 0.00000000 ; 
    expected_f_mu_negative[0] = 0.00000000 ; 
    expected_f_mu_negative[1] = -0.00000000 ; 
    expected_f_mu_negative[2] = 0.00000000 ; 
    expected_f_mu_negative[3] = -1.00000000 ; 
        
    matrix_construction->construct_dimensionless_mass_matrix(matrix);
    for(int i= 0 ; i < n_dfem_p ; i++)
    {
      for(int j=0; j < n_dfem_p ; j++)
      {
        if(fabs( (expected_dimensionless_mass[i][j] - matrix(i,j) )) > tol )
        {
          std::stringstream err;
          err << "Discrepancy calculating dimensionless mass matrix in row: " << i << " column: " << j <<std::endl;
          throw Dark_Arts_Exception(FEM, err.str() );
        }          
      }
    }
    
    matrix_construction->construct_pos_gradient_matrix(matrix);
    for(int i= 0 ; i < n_dfem_p ; i++)
    {
      for(int j=0; j < n_dfem_p ; j++)
      {
        if(fabs( (expected_L_mu_positive[i][j] - matrix(i,j) )) > tol )
        {
          std::stringstream err;
          err << "Discrepancy calculating L for mu > 0 in row: " << i << " column: " << j <<std::endl;
          throw Dark_Arts_Exception(FEM, err.str() );
        }          
      }
    }
    
    matrix_construction->construct_neg_gradient_matrix(matrix);
    for(int i= 0 ; i < n_dfem_p ; i++)
    {
      for(int j=0; j < n_dfem_p ; j++)
      {
        if(fabs( (expected_L_mu_negative[i][j] - matrix(i,j) )) > tol )
        {
          std::stringstream err;
          err << "Discrepancy calculating L for mu < 0 in row: " << i << " column: " << j <<std::endl;
          throw Dark_Arts_Exception(FEM, err.str() );
        }          
      }
    }
    
    matrix_construction->construct_neg_upwind_vector(vector);
    for(int i= 0 ; i < n_dfem_p ; i++)
    {
      if(fabs( (expected_f_mu_negative[i] - vector(i) )) > tol )
      {
        std::stringstream err;
        err << "Discrepancy calculating f for mu < 0 in row: " << i << std::endl;
        throw Dark_Arts_Exception(FEM, err.str() );
      }  
    }
    
    matrix_construction->construct_pos_upwind_vector(vector);
    for(int i= 0 ; i < n_dfem_p ; i++)
    {
      if(fabs( (expected_f_mu_positive[i] - vector(i) )) > tol )
      {
        std::stringstream err;
        err << "Discrepancy calculating f for mu > 0 in row: " << i << std::endl;
        throw Dark_Arts_Exception(FEM, err.str() );
      }  
    }  
      
    
  }
  catch(const Dark_Arts_Exception da_exception)
  {
    da_exception.testing_message() ;
    val = -1;
  }  

  
 
    
  // Return 0 if tests passed, something else if failing
  return val;
}


