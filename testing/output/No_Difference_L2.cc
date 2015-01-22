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
#include "L2_Error_Calculator.h"

#include "Dark_Arts_Exception.h"

/**
  Goal of this unit test is to verify that L2 error estimate is zero when the analytic solution is within the DFEM trial space
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
  
  /// Create a Materials object that contains all opacity, heat capacity, and source objects
  Materials materials( input_reader, fem_quadrature , cell_data, angular_quadrature );  
  Time_Data time_data(input_reader);
  
  Temperature_Data t_old(fem_quadrature, input_reader, cell_data);  
  Intensity_Data i_old(cell_data, angular_quadrature, fem_quadrature, materials,input_reader);    
  Intensity_Moment_Data phi(cell_data,angular_quadrature,fem_quadrature,i_old);  
  
  
  
  try{ 
    L2_Error_Calculator l2_error_calculator(angular_quadrature,fem_quadrature, cell_data, input_reader);
    const double t_eval = time_data.get_t_start();
    double temperature_err = l2_error_calculator.calculate_l2_error(t_eval , t_old);
    double phi_err = l2_error_calculator.calculate_l2_error(t_eval , phi);
  
    std::cout << "Phi L2 err: " << phi_err << std::endl;
    std::cout << "Temperature L2 err: " << temperature_err << std::endl;
      
    if(fabs(phi_err > tol))
      throw Dark_Arts_Exception(SUPPORT_OBJECT , "Not calculating phi zero err");
      
    if(fabs(temperature_err > tol))
      throw Dark_Arts_Exception(SUPPORT_OBJECT , "Not calculating temperature zero err");
  
    temperature_err = l2_error_calculator.calculate_cell_avg_error(t_eval , t_old);
    phi_err = l2_error_calculator.calculate_cell_avg_error(t_eval , phi);
  
    std::cout << "Phi avg err: " << phi_err << std::endl;
    std::cout << "Temperature avg err: " << temperature_err << std::endl;
      
    if(fabs(phi_err > tol))
      throw Dark_Arts_Exception(SUPPORT_OBJECT , "Not calculating phi zero err");
      
    if(fabs(temperature_err > tol))
      throw Dark_Arts_Exception(SUPPORT_OBJECT , "Not calculating temperature zero err");

  }
  catch(const Dark_Arts_Exception& da)
  {
    val = -1;
    da.testing_message();
  }
  
  return val;
}
