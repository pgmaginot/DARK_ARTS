#include <iostream>
#include <vector>
#include <math.h>
#include <string>

#include "Input_Reader.h"
#include "Cell_Data.h"
#include "Intensity_Data.h"
#include "Angular_Quadrature.h"
#include "Materials.h"
#include "Fem_Quadrature.h"
#include "Quadrule_New.h"
#include "Dark_Arts_Exception.h"

int main(int argc, char** argv)
{
  int val = 0;
  const double tol = 1.0E-6;
  
  /// expected spacing we should be able to see
  std::vector<int> expected_n_cells(2,0);
  expected_n_cells[0] = 5;
  expected_n_cells[1] = 3;
  const double sn_w = 2.;
  const double I_0 = pow(0.4, 4 )/sn_w;
  const double I_1 = pow(0.98,4 )/sn_w;
  
  Input_Reader input_reader;    
    
  try{
    input_reader.read_xml(argv[1]);
  }
  catch( const Dark_Arts_Exception& da_exception )
  {
    da_exception.testing_message();
    val = -1;
  }
  
  Cell_Data cell_data( input_reader );   
  Quadrule_New quad_fun;   
  Fem_Quadrature fem_quadrature( input_reader , quad_fun);    
  Angular_Quadrature angular_quadrature( input_reader , quad_fun );  
  Materials materials( input_reader, fem_quadrature , cell_data, angular_quadrature);  
      
  const int n_p = fem_quadrature.get_number_of_interpolation_points();
  Eigen::VectorXd local_i_vec;
  Eigen::VectorXd local_zero_vec;
  
  Intensity_Data intensity_ic(cell_data, angular_quadrature, fem_quadrature, materials,input_reader);
  Intensity_Data intensity_zero(cell_data, angular_quadrature, fem_quadrature);
  
  int cell_cnt = 0;
  const int n_dir = angular_quadrature.get_number_of_dir();
  const int n_grp = 1;
  try
  {
    for(int reg=0; reg < 2 ; reg++)
    {
      double I;
      if( reg==0)
        I = I_0;
      else
        I = I_1;
        
      for(int cell = 0; cell < expected_n_cells[reg] ; cell++)
      {
        for(int dir = 0; dir < n_dir ; dir++)
        {       
          local_i_vec = Eigen::VectorXd::Zero(n_p);
          local_zero_vec = Eigen::VectorXd::Ones(n_p);
          intensity_ic.get_cell_intensity( cell_cnt, n_grp-1 , dir, local_i_vec);
          intensity_zero.get_cell_intensity(cell_cnt,n_grp-1,dir,local_zero_vec);
          for(int el = 0 ; el < n_p ; el++)
          {
            if( fabs( local_i_vec(el) - I) > tol )
              throw Dark_Arts_Exception(VARIABLE_STORAGE , "IC intensity different than expected");
              
            if( fabs(local_zero_vec(el) ) > tol )
              throw Dark_Arts_Exception(VARIABLE_STORAGE , "Zero intensity different than expected");
          }
          
          local_zero_vec(0) = double(cell_cnt) + 10.*double(dir+1) + 0.1; 
          local_zero_vec(1) = double(cell_cnt) + 10.*double(dir+1) + 0.2; 
          local_zero_vec(2) = double(cell_cnt) + 10.*double(dir+1) + 0.3; 
          local_zero_vec(3) = double(cell_cnt) + 10.*double(dir+1) + 0.4; 
          intensity_zero.set_cell_intensity(cell_cnt,n_grp-1,dir,local_zero_vec);
        }
        cell_cnt++;
      }
    }
    
    std::cout << "Checking Set/Get\n";
    
    /// check the set get's again, independently
    for(int d = (n_dir - 1); d > -1 ; d--)
    {
      for(int cell = 0; cell < 8 ; cell++)
      {
        intensity_zero.get_cell_intensity(cell,n_grp-1,d,local_zero_vec);
        for(int el = 0; el < n_p ; el++)
        {
          if( fabs(local_zero_vec(el) - (double(cell) + 10.*double(d+1) + 0.1*double(el+1) ) ) > tol)
            throw Dark_Arts_Exception(VARIABLE_STORAGE , "Set/get problems"); 
        }      
      }
    }
  } 
  catch(const Dark_Arts_Exception& da_exception )
  {
    da_exception.testing_message();
    val = -1;
  }
  
  // Return 0 if tests passed, somethnig else if failing
  return val;
}
