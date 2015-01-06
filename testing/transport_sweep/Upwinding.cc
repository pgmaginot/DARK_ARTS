#include <iostream>
#include <vector>
#include <math.h>
#include <string>

#include "Input_Reader.h"
#include "Cell_Data.h"
#include "Intensity_Data.h"
#include "Angular_Quadrature.h"
#include "Materials.h"
#include "Fem_Quadrature.h"
#include "Quadrule_New.h"
#include "Solution_Saver_Flux_Moments.h"
#include "MMS_Intensity.h"

#include "Dark_Arts_Exception.h"

int main(int argc, char** argv)
{
  int val = 0;
  const double tol = 1.0E-6;
    
  /**
    Test out V_Solution_Saver::calculate_outflow
  */  
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
  
  Solution_Saver_Flux_Moments saver(fem_quadrature,angular_quadrature);

  std::vector<double> dfem_pt;
  fem_quadrature.get_dfem_interpolation_point(dfem_pt);

  MMS_Intensity i_mms(input_reader,angular_quadrature);
  
  const double xL = cell_data.get_cell_left_edge(0);
  const double dx = cell_data.get_cell_width(0);
  
  const double t_start = 1.5;
  const int n_p = fem_quadrature.get_number_of_interpolation_points();
  Eigen::VectorXd local_i_vec = Eigen::VectorXd::Zero(n_p);
  
  try{
  
    int dir = 0;
    double x_dfem;
    for(int el = 0 ; el < n_p ; el++)
    {       
      x_dfem = xL + dx/2.*(1. + dfem_pt[el]);
      local_i_vec(el) = i_mms.get_mms_intensity(x_dfem,t_start,dir);
    }
    
    const double I_left = i_mms.get_mms_intensity(xL,t_start,dir);
    const double calc_i_left = saver.calculate_outflow(dir,local_i_vec);
      
    dir = angular_quadrature.get_number_of_dir() - 1;
    for(int el = 0 ; el < n_p ; el++)
    {       
      x_dfem = xL + dx/2.*(1. + dfem_pt[el]);
      local_i_vec(el) = i_mms.get_mms_intensity(x_dfem,t_start,dir);
    }
    
    const double I_right = i_mms.get_mms_intensity(xL+dx,t_start,dir);
    const double calc_i_right = saver.calculate_outflow(dir,local_i_vec);
    std::cout << "Expected I_left: " << I_left << " calculated i_left: " << calc_i_left << std::endl;
    std::cout << "Expected I_right: " << I_right << " calculated i_right: " << calc_i_right << std::endl;
    if(fabs(calc_i_left - I_left ) > tol)
      throw Dark_Arts_Exception(TRANSPORT, "mu < 0 upwinding is wrong");
    
    if(fabs(calc_i_right - I_right ) > tol)
      throw Dark_Arts_Exception(TRANSPORT, "mu > 0 upwinding is wrong");
    
  }
  catch(const Dark_Arts_Exception& da_exception )
  {
    da_exception.testing_message();
    val = -1;
  }
  
  // Return 0 if tests passed, somethnig else if failing
  return val;
}
