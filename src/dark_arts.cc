#include "Input_Reader.h"
#include "Fem_Quadrature.h"
#include "Cell_Data.h"
#include "Time_Stepper.h"
#include "Angular_Quadrature.h"
#include "Intensity_Data.h"
#include "Materials.h"

#include "Eigen/Dense"

int main(int argc, char** argv)
{
  std::cout << "argc = " << argc << '\n'; 
  for(int i = 0; i < argc; i++) 
    std::cout << "argv[" << i << "] = " << argv[i] << '\n'; 
  
  Input_Reader input_reader;
  bool input_parsed = false;
    
  input_parsed = input_reader.read_xml(argv[1]);
    
  if(!input_parsed)
  {
    std::cerr << "Error reading input" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  /// Initialize a Quadrule object to be able to get all of the quadrature we need
  Quadrule_New quad_fun;
  
  /// Initialize FEM data
  /// get all interpolation points, quadrature formuals, matrix formation routines, etc.
  Fem_Quadrature fem_quadrature( input_reader , quad_fun);
    
  /// Initalize cell data (dx, xL, xR, x_ip, material_number)
  Cell_Data cell_data( input_reader );
  
  /// Initialize time-stepping scheme (SDIRK method)
  Time_Stepper time_stepper( input_reader );
  
  /// Initialize angular quadrature data.  Will include number of: directions, groups, and legendre moments.
  /// will also include evaluations of Legendre polynomials
  Angular_Quadrature angular_quadrature( input_reader , quad_fun );
  
  /// Initialize intensity and angle integrated intensity of previous time step
  Intensity_Data intensity_old( cell_data, angular_quadrature, fem_quadrature);
  
  /// Create a Materials object that contains all opacity, heat capacity, and source objects
  Materials materials( input_reader, fem_quadrature , &cell_data);
    
  Eigen::MatrixXd m(2,2);
  
  

  
  
  return 0;
}

