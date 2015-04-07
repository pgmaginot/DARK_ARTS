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
#include "Time_Marcher.h"

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
  }     
  std::cout << "Input File Read" << std::endl;

  try{
    /// Initialize a Quadrule object to be able to get all of the quadrature we need
    Quadrule_New quad_fun;  
    std::cout << "Quadrule object created" << std::endl;

    /// Initialize FEM data
    /// get all interpolation points, quadrature formuals, matrix formation routines, etc.
    Fem_Quadrature fem_quadrature( input_reader , quad_fun);  
    std::cout << "Fem_Quadrature object created" << std::endl;

    /// Initalize cell data (dx, xL, xR, x_ip, material_number)
    Cell_Data cell_data( input_reader );  
    std::cout << "Cell_Data object created" << std::endl;

    /// Initialize angular quadrature data.  Will include number of: directions, groups, and legendre moments.
    /// will also include evaluations of Legendre polynomials
    Angular_Quadrature angular_quadrature( input_reader , quad_fun );  
    std::cout << "Angular quadrature object created" << std::endl;
     
    /// Create a Materials object that contains all opacity, heat capacity, and source objects
    Materials materials( input_reader, fem_quadrature , cell_data, angular_quadrature );  
    std::cout << "Materials object created" << std::endl;

    /// Load SDIRK data
    Time_Data time_data( input_reader);  
    std::cout << "Time data object created" << std::endl;

    /// Initialize intensity and angle integrated intensity of previous time step
    Intensity_Data intensity_old( cell_data, angular_quadrature, fem_quadrature, materials, input_reader);
    std::cout << "Intensity object created" << std::endl;

    /// Initialize a Temperature_Data structure
    Temperature_Data temperature_old( fem_quadrature, input_reader, cell_data);
    std::cout << "Temperature object created" << std::endl;

    std::string input_filename = argv[1];
    unsigned int found = input_filename.find_last_of("/");
    std::string short_input_filename = input_filename.substr(found+1);  
    /// Time Marcher.  This object will take care of solving the problem
    Time_Marcher time_marcher(input_reader, angular_quadrature,fem_quadrature,
    cell_data, materials, temperature_old, intensity_old, time_data, short_input_filename);    
    std::cout << "Time_Marcher object created"<< std::endl;  
  }
  catch(const Dark_Arts_Exception& da)
  {
    val = -1;
    da.testing_message();
  }
  
  

  
  return val;
}
