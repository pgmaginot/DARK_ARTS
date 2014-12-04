#include <iostream>
#include <vector>
#include <math.h>
#include <string>

#include "Input_Reader.h"
#include "Cell_Data.h"
#include "Temperature_Data.h"
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
  const double t_0 = 0.5;
  const double t_1 = 0.27;
  
  const double t_avg = (double(expected_n_cells[0])*t_0 + double(expected_n_cells[1])*t_1)/(double (expected_n_cells[0] + expected_n_cells[1]));
  
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
    
  Temperature_Data temperature_ic( fem_quadrature, input_reader, cell_data);   

  Temperature_Data temperature_zero( cell_data.get_total_number_of_cells() , fem_quadrature);
  
  Temperature_Data temperature_copy(cell_data.get_total_number_of_cells() , fem_quadrature);
  temperature_copy = temperature_ic;
  
  const int n_p = fem_quadrature.get_number_of_interpolation_points();
  Eigen::VectorXd local_t_vec;
  
  int cell_cnt = 0;
  try
  {
    for(int reg=0; reg < 2 ; reg++)
    {
      double t;
      if( reg==0)
        t = t_0;
      else
        t = t_1;
        
      for(int cell = 0; cell < expected_n_cells[reg] ; cell++)
      {
        local_t_vec = Eigen::VectorXd::Zero(n_p);
        temperature_ic.get_cell_temperature( cell_cnt, local_t_vec);
        for(int el = 0 ; el < n_p ; el++)
        {
          if( fabs( local_t_vec(el) - t) > tol )
            throw Dark_Arts_Exception(VARIABLE_STORAGE , "IC temperature different than expected");
        }
        
        local_t_vec = Eigen::VectorXd::Zero(n_p);
        temperature_copy.get_cell_temperature( cell_cnt, local_t_vec);
        for(int el = 0 ; el < n_p ; el++)
        {
          if( fabs( local_t_vec(el) - t) > tol )
            throw Dark_Arts_Exception(VARIABLE_STORAGE , "Copied temperature different than expected");
        }
        
        local_t_vec = Eigen::VectorXd::Ones(n_p);
        temperature_zero.get_cell_temperature(cell_cnt, local_t_vec);
        for(int el = 0 ; el < n_p ; el++)
        {
          if( fabs( local_t_vec(el) ) > tol )
            throw Dark_Arts_Exception(VARIABLE_STORAGE , "Zero temperature different than expected");
        }
        
        cell_cnt++;
      }
    }
    
    if(fabs(temperature_zero.calculate_average() ) > tol)
      throw Dark_Arts_Exception(VARIABLE_STORAGE , "Incorrect average for zero temperature");
    
    if(fabs(temperature_ic.calculate_average() - t_avg) > tol)
      throw Dark_Arts_Exception(VARIABLE_STORAGE , "Incorrect average for IC temperature");
    
    if(fabs(temperature_copy.calculate_average() - t_avg) > tol)
      throw Dark_Arts_Exception(VARIABLE_STORAGE , "Incorrect average for copied temperature");
      
            /// test out the set_cell_temperature functionality
    cell_cnt = 0;
    for(int reg=0; reg < 2 ; reg++)
    {  
      for(int cell = 0; cell < expected_n_cells[reg] ; cell++)
      {
        for(int el = 0; el < n_p ; el++)
          local_t_vec(el) = double(cell_cnt) + double(el)*0.1;
          
        temperature_ic.set_cell_temperature(cell_cnt, local_t_vec);
        cell_cnt++;
      }
    }
     
    cell_cnt = 0;
    for(int cell = 0; cell < 8 ; cell++)
    {
      local_t_vec(0) = double(cell) + 0.1 ;
      local_t_vec(1) = double(cell) + 0.1 ;
      local_t_vec(2) = double(cell) + 0.1 ;
      local_t_vec(3) = double(cell) + 0.1 ;
          
      temperature_ic.get_cell_temperature(cell_cnt, local_t_vec);
      for(int el = 0; el < n_p ; el++)
      {
        double exp_val = double(cell_cnt) + 0.1*double(el);
        if(fabs( local_t_vec(el) - exp_val) > tol)
          throw Dark_Arts_Exception(VARIABLE_STORAGE , "Trouble setting and retrieving the correct/same temperature values");
      }        
      cell_cnt++;
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
