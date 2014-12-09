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
  
  std::vector<double> ex_temp_space_coeff(4,0.);
  std::vector<double> ex_time_coeff(2,0.);
    
  ex_temp_space_coeff[0] = 1.;
  ex_temp_space_coeff[1] = 2.;
  ex_temp_space_coeff[2] = 0.;
  ex_temp_space_coeff[3] = 3.;
  
  ex_time_coeff[0] = 0.5;
  ex_time_coeff[1] = 1.2;

  /// verify composite manufactured  solutions for a single point in space/time
  std::shared_ptr<V_MMS_Time> time_component = std::shared_ptr<V_MMS_Time> (new MMS_Time_Poly(ex_time_coeff) );
  std::shared_ptr<V_MMS_Space> temp_space_component = std::shared_ptr<V_MMS_Space> (new MMS_Space_Cos(ex_temp_space_coeff) );
  
  // <N_cells> 5 </N_cells>
      // <Left_bound> 1.0 </Left_bound>
      // <Right_bound> 22.0 </Right_bound>
      // <Spacing> Equal   </Spacing>
      // <Material_number> 0 </Material_number>
  
  const double t_start = 1.5;
  
  const int n_p = fem_quadrature.get_number_of_interpolation_points();
  Eigen::VectorXd local_t_vec;
  double ex_t, x_dfem , dx , xL;
  std::vector<double> dfem_interp;
  fem_quadrature.get_dfem_interpolation_point(dfem_interp);
  
  try
  {
    for(int cell = 0; cell < 5 ; cell++)
    {
      local_t_vec = Eigen::VectorXd::Zero(n_p);
      ex_t = 0.;
      dx = cell_data.get_cell_width(cell);
      xL = cell_data.get_cell_left_edge(cell);
      
      temperature_ic.get_cell_temperature( cell, local_t_vec);
      std::cout << "Cell: " << cell << std::endl;
      for(int el = 0 ; el < n_p ; el++)
      {
        x_dfem = xL + dx/2.*(1. + dfem_interp[el]);
        ex_t = time_component->get_time_component(t_start)*temp_space_component->get_position_component(x_dfem);
        std::cout << "calc temp: " << local_t_vec(el) << " expected temp: " << ex_t << std::endl;
        if( fabs( local_t_vec(el) - ex_t) > tol )
          throw Dark_Arts_Exception(VARIABLE_STORAGE , "IC temperature different than expected");
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
