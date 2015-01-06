#include <iostream>
#include <vector>
#include <math.h>
#include <string>

#include "Input_Reader.h"
#include "Cell_Data.h"
#include "Intensity_Data.h"
#include "Intensity_Moment_Data.h"
#include "Angular_Quadrature.h"
#include "Materials.h"
#include "Fem_Quadrature.h"
#include "Quadrule_New.h"
#include "Dark_Arts_Exception.h"

int main(int argc, char** argv)
{
  int val = 0;
  const double tol = 1.0E-6;
  
  /// expected spacing we should be able to see
  
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
  Angular_Quadrature angular_quadrature( input_reader , quad_fun );  
  Materials materials( input_reader, fem_quadrature , cell_data, angular_quadrature);  
      
  const int n_p = fem_quadrature.get_number_of_interpolation_points();
  Eigen::VectorXd local_i_vec;
  
  Intensity_Data intensity_ic(cell_data, angular_quadrature, fem_quadrature, materials,input_reader);  
  
  Intensity_Moment_Data phi_ic( cell_data,angular_quadrature, fem_quadrature, intensity_ic );
  Eigen::VectorXd ex_phi_vec = Eigen::VectorXd::Zero(n_p);
  Eigen::VectorXd calc_phi_vec = Eigen::VectorXd::Zero(n_p);
  
  const int n_dir = angular_quadrature.get_number_of_dir();
  
  double x_dfem , xL, dx , I;
  std::vector<double> dfem_pt;
  fem_quadrature.get_dfem_interpolation_point(dfem_pt);
  
  std::vector<double> ex_rad_space_coeff(5,0.);
  std::vector<double> ex_rad_angle_coeff(3,0.);
  std::vector<double> ex_time_coeff(2,0.);
    
  ex_rad_space_coeff[0] = 0.6;
  ex_rad_space_coeff[1] = 2.1;
  ex_rad_space_coeff[2] = 4.5;
  ex_rad_space_coeff[3] = 5.;
  ex_rad_space_coeff[4] = 3.2;
  
  ex_rad_angle_coeff[0] = 1.1;
  ex_rad_angle_coeff[1] = 0.1;
  ex_rad_angle_coeff[2] = 1.3;
      
  ex_time_coeff[0] = 0.5;
  ex_time_coeff[1] = 1.2;
  
  /// verify composite manufactured  solutions for a single point in space/time
  std::shared_ptr<V_MMS_Time> time_component = std::shared_ptr<V_MMS_Time> (new MMS_Time_Poly(ex_time_coeff) );
  std::shared_ptr<V_MMS_Space> rad_space_component = std::shared_ptr<V_MMS_Space> (new MMS_Space_Poly(ex_rad_space_coeff) );
  std::shared_ptr<V_MMS_Angle> rad_angle_component = std::shared_ptr<V_MMS_Angle> (new MMS_Angle_Poly(ex_rad_angle_coeff, angular_quadrature) );
  
  const double t_start = 1.5;
  
  try
  {    
    for(int cell = 0; cell < 5 ; cell++)
    {
      xL = cell_data.get_cell_left_edge(cell);
      dx = cell_data.get_cell_width(cell);
      ex_phi_vec = Eigen::VectorXd::Zero(n_p);
      for(int dir = 0; dir < n_dir ; dir++)
      {       
        local_i_vec = Eigen::VectorXd::Zero(n_p);
        intensity_ic.get_cell_intensity( cell, 0 , dir, local_i_vec);
        
        for(int el = 0 ; el < n_p ; el++)
        {          
          x_dfem = xL + dx/2.*(1. + dfem_pt[el]);
          I = time_component->get_time_component(t_start);
          I *= rad_space_component->get_position_component(x_dfem);
          I *= rad_angle_component->get_angle_component(dir);
          
          std::cout << "Calc I: " << local_i_vec(el) << " Expected I: " << I << std::endl;
          if( fabs( local_i_vec(el) - I) > tol )
            throw Dark_Arts_Exception(VARIABLE_STORAGE , "IC intensity different than expected");            
        }
        ex_phi_vec += angular_quadrature.get_w(dir) * local_i_vec;
      }
      phi_ic.get_cell_angle_integrated_intensity(cell,0,0,calc_phi_vec);
      std::cout << "Difference between phi: " << std::endl;
      local_i_vec = ex_phi_vec - calc_phi_vec;
      std::cout << local_i_vec << std::endl; 
      for(int el = 0 ; el < n_p ; el++)
      { 
        if( fabs( ex_phi_vec(el) - calc_phi_vec(el)) > tol )
          throw Dark_Arts_Exception(VARIABLE_STORAGE , "phi IC different than expected");            
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
