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
#include "Temperature_Data.h"
#include "Output_Generator.h"

#include "Dark_Arts_Exception.h"

/**
  Goal of this unit test is to dump a Intensity and Temperature Restart files
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
  
  try{
    std::string input_filename = argv[1];
    unsigned int found = input_filename.find_last_of("/");
    std::string short_input_filename = input_filename.substr(found+1);  
  
    Output_Generator output_generator(angular_quadrature, fem_quadrature, cell_data,input_reader);
    
    for(int c = 0 ; c < n_cell ; c++)
    {
      Eigen::VectorXd t_loc = Eigen::VectorXd::Zero(n_p);
      t_old.get_cell_temperature(c,t_loc);
      for(int el = 0 ; el < n_p; el++)
      {
        std::cout << t_loc(el) << "  " ;
      }
      std::cout << std::endl;
    }
    
    output_generator.write_xml(false,1,t_old);
    output_generator.write_xml(false,1,i_old);

      std::cout << "ts_1 files should now be in /build diretory\n";
  }
  catch(const Dark_Arts_Exception& da)
  {
    val = -1;
    da.testing_message();
  }
  
  

  
  return val;
}
