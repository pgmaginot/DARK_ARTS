static char help[] = "A character array that PETSc might be expecting";
#include <petscksp.h>

#include "Input_Reader.h"
#include "Output_Generator.h"
#include "Fem_Quadrature.h"
#include "Cell_Data.h"
#include "Time_Data.h"
#include "Time_Marcher.h"
#include "Angular_Quadrature.h"
#include "Intensity_Data.h"
#include "Temperature_Data.h"
#include "Materials.h"

#include "Dark_Arts_Exception.h"


int main(int argc, char** argv)
{
  /**
    Initialize PETSc
  */
  PetscErrorCode ierr;  
  PetscMPIInt size;
  
  PetscInitialize(&argc,&argv,(char*)0,help);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);
  CHKERRQ(ierr);
  // if (size != 1) 
    // SETERRQ(PETSC_COMM_WORLD,1,"DARK_ARTS is written to be serial only!");
    
  std::cout << "argc = " << argc << '\n'; 
  for(int i = 0; i < argc; i++) 
    std::cout << "argv[" << i << "] = " << argv[i] << '\n'; 

  /// DARK_ARTS objects follow
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

  
  /// Time Marcher.  This object will take care of solving the problem
  Time_Marcher time_marcher(input_reader, angular_quadrature,fem_quadrature,
    cell_data, materials, temperature_old, intensity_old, time_data);    
    std::cout << "Time_Marcher object created"<< std::endl;  
  
  /// this is the entire time loop !
  try{
    time_marcher.solve(intensity_old, temperature_old, time_data);
  }
  catch(const Dark_Arts_Exception& da_exception)
  {
    da_exception.message() ;
  }
  
  /**
    End PETSc
  */
  ierr = PetscFinalize();
 
  return 0;
}