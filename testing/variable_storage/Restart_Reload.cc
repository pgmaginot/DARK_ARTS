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
  Goal of this unit test is to read a dumped out input file
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
  
  const int n_cell = cell_data.get_total_number_of_cells();
  std::cout << "N_cell: " << n_cell << std::endl;
  const int n_p = fem_quadrature.get_number_of_interpolation_points();  
  const int n_dir = angular_quadrature.get_number_of_dir();
  /// Create a Materials object that contains all opacity, heat capacity, and source objects
  Materials materials( input_reader, fem_quadrature , cell_data, angular_quadrature );  
  Time_Data time_data(input_reader);
  
  try{
  
    Temperature_Data t_old(fem_quadrature, input_reader, cell_data);  
    Intensity_Data i_old(cell_data, angular_quadrature, fem_quadrature, materials,input_reader);
    
    MMS_Intensity intensity_reference(input_reader, angular_quadrature);
    MMS_Temperature temperature_reference(input_reader);
    
    std::vector<double> dfem_pts;
    fem_quadrature.get_dfem_interpolation_point(dfem_pts);
    const double time = time_data.get_t_start();
    for(int c = 0 ; c < n_cell ; c++)
    {
      double dx = cell_data.get_cell_width(c);
      double xL = cell_data.get_cell_left_edge(c);
      Eigen::VectorXd t_local = Eigen::VectorXd::Zero(n_p);
      t_old.get_cell_temperature(c,t_local);
      std::cout << " Cell: " << c << " dx: " << dx << " xL: " << xL << std::endl; 
      for( int el = 0 ; el < n_p ; el++)
      {
        double x_local = xL + dx/2.*(1. + dfem_pts[el]);
        double t_ref = temperature_reference.get_mms_temperature(x_local,time);
        std::cout << "Expected temperature: " << t_ref << " Restarted temperature: " << t_local(el) << std::endl;
        if( isnan(t_local(el) ))
          throw Dark_Arts_Exception(VARIABLE_STORAGE, "NAN in temperature reload");
        
        if(fabs(t_local(el) - t_ref) > tol)
          throw Dark_Arts_Exception(VARIABLE_STORAGE, "Restart temperature is incorrect");
      }
      
      for(int d = 0 ; d < n_dir ; d++)
      {
        Eigen::VectorXd i_local = Eigen::VectorXd::Zero(n_p);
        i_old.get_cell_intensity(c,0,d,i_local);
        for( int el = 0 ; el < n_p ; el++)
        {
          double x_local = xL + dx/2.*(1. + dfem_pts[el]);
          double i_ref = intensity_reference.get_mms_intensity(x_local,time,d);
          std::cout << "Expected intensity: " << i_ref << " Restarted intensity: " << i_local(el) << std::endl;
          if( isnan(i_local(el) ))
            throw Dark_Arts_Exception(VARIABLE_STORAGE, "NAN in intensity reload");
          
          if(fabs(i_local(el) - i_ref) > tol)
            throw Dark_Arts_Exception(VARIABLE_STORAGE, "Restart intensity is incorrect");
        }
        
      }
      
    }
  }
  catch(const Dark_Arts_Exception& da)
  {
    val = -1;
    da.testing_message();
  }
  
  

  
  return val;
}
