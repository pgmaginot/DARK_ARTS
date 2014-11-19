#include <iostream>
#include <vector>
#include <math.h>
#include <string>

#include "Input_Reader.h"
#include "Cell_Data.h"
#include "Fem_Quadrature.h"
#include "Angular_Quadrature.h"
#include "Quadrule_New.h"
#include "Materials.h"
#include "Intensity_Data.h"
#include "Temperature_Data.h"
#include "Intensity_Moment_Data.h"
#include "Output_Generator.h"
#include "Dark_Arts_Exception.h"


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
  }     
  std::cout << "Input File Read" << std::endl;
  
  Quadrule_New quad_fun;  
  Fem_Quadrature fem_quadrature( input_reader , quad_fun);  
  Cell_Data cell_data( input_reader );  
  Angular_Quadrature angular_quadrature( input_reader , quad_fun );  
  Materials materials( input_reader, fem_quadrature , cell_data, angular_quadrature.get_number_of_groups() , angular_quadrature.get_sum_w() );  
  Intensity_Data intensity_old( cell_data, angular_quadrature, fem_quadrature, materials, input_reader);
  Temperature_Data temperature_old( cell_data.get_total_number_of_cells(), fem_quadrature, input_reader);
  Intensity_Moment_Data phi_ic(cell_data,angular_quadrature, fem_quadrature, intensity_old);
  
  std::string my_str = "Does nothing";
  Output_Generator output(angular_quadrature,
    fem_quadrature, cell_data, argv[1]);
  
  try{
    output.write_xml( false, 1, temperature_old);
    output.write_xml( false, 1, phi_ic);
    output.write_xml( false, 1, intensity_old);
    
    output.write_txt( false, 1, phi_ic);
    output.write_txt( false, 1, temperature_old);
    output.write_txt( false, 1, intensity_old);
  }
  catch(const Dark_Arts_Exception& da)
  {
    val = -1;
    da.testing_message();
  }
  
  // Return 0 if tests passed, somethnig else if failing
  return val;
}
