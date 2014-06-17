#include "Input_Reader.h"
#include "Fem_Quadrature.h"
#include "Cell_Data.h"

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
  
  /// Initialize FEM data
  /// get all interpolation points, quadrature formuals, matrix formation routines, etc.
  Fem_Quadrature fem_quadrature( input_reader );
  
  std::cout << "fem_quadrature_object created\n";
  
  /// Initalize cell data (dx, xL, xR, x_ip, material_number)
  Cell_Data cell_data( input_reader );
  
  /// Get an array of Material objects
 // Materials_Vector mat_data( &input_reader);
  
  

  
  
  return 0;
}