static char help[] = "A character array that PETSc might be expecting";
#include <petscksp.h>

#include "Input_Reader.h"
#include "Fem_Quadrature.h"
#include "Cell_Data.h"
#include "Time_Data.h"
#include "Time_Marcher.h"
#include "Angular_Quadrature.h"
#include "Intensity_Data.h"
#include "Temperature_Data.h"
#include "Materials.h"
#include "Status_Generator.h"

#include "Dark_Arts_Exception.h"
#include <chrono>


int main(int argc, char** argv)
{
  auto t1 = std::chrono::high_resolution_clock::now();
  
  int val = 0;
  /**
    Initialize PETSc
  */
  try{
    PetscErrorCode ierr;  
    PetscMPIInt size;
    
    PetscInitialize(&argc,&argv,(char*)0,help);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);
    CHKERRQ(ierr);
    if (size != 1) 
      SETERRQ(PETSC_COMM_WORLD,1,"DARK_ARTS is written to be serial only!");
      
    std::cout << "argc = " << argc << '\n'; 
    for(int i = 0; i < argc; i++) 
      std::cout << "argv[" << i << "] = " << argv[i] << '\n'; 

    /// DARK_ARTS objects follow
    Input_Reader input_reader;
    input_reader.read_xml(argv[1]);  
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

    /// create filename for status generator
    /// for mesh refinement runs: output_dir/short_base_filename_refinement.xml
    /// otherwise: output_dir/
    
    std::string stat_file_path_and_name;
    if( input_reader.is_mesh_refinement() )
    {
      std::string base_file_with_path = input_reader.get_initial_input_filename();
      unsigned int found = base_file_with_path.find_last_of("/");
      std::string base_short = base_file_with_path.substr(found+1);  
      
      stat_file_path_and_name = input_reader.get_output_directory();
      stat_file_path_and_name.append(base_short);
      std::string xml_ext = ".xml";
      stat_file_path_and_name.replace( stat_file_path_and_name.find(xml_ext) , xml_ext.length() , "_");
      
      std::string input_filename = argv[1];
      found = input_filename.find_last_of("/");
      std::string short_input_filename = input_filename.substr(found+1); 
      
      stat_file_path_and_name.append( short_input_filename );
    }
    else
    {
      std::string input_filename = argv[1];
      unsigned int found = input_filename.find_last_of("/");
      std::string short_input_filename = input_filename.substr(found+1);  
      std::string output_directory;
      input_reader.get_output_directory(output_directory);
      // output_directory.append(short_input_filename);
      stat_file_path_and_name = output_directory + short_input_filename;
    }
    
    
    /// Time Marcher.  This object will take care of solving the problem
    Time_Marcher time_marcher(input_reader, angular_quadrature,fem_quadrature,
      cell_data, materials, temperature_old, intensity_old, time_data,stat_file_path_and_name);    
      std::cout << "Time_Marcher object created"<< std::endl;  
  
  /// this is the entire time loop !

    time_marcher.solve(intensity_old, temperature_old);
    /**
      End PETSc
    */
    ierr = PetscFinalize();
  }
  catch(const Dark_Arts_Exception& da_exception)
  {
    val = -1;
    da_exception.message() ;
  }
  
  auto t2 = std::chrono::high_resolution_clock::now();
  
   std::cout << "DARK_ARTS took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()
              << " milliseconds\n";
 
  return val;
}