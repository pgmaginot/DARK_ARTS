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
#include "MMS_Temperature.h"
#include "MMS_Intensity.h"
#include "Transport_BC_MMS.h"

#include "Dark_Arts_Exception.h"

/**
  Goal of this unit test is to check source moment formation, and reaction matrix for SLXS Lobatto scheme
  -This is a MMS problem with a spatially varying sigma_a, constant cv , zero sig_s
*/ 

int main(int argc, char** argv)
{
  int val = 0;
  const double tol = 1.0E-6;
  
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
  
  /// Initialize a Quadrule object to be able to get all of the quadrature we need
  Quadrule_New quad_fun;  
  Fem_Quadrature fem_quadrature( input_reader , quad_fun);  
  Cell_Data cell_data( input_reader );  
  Angular_Quadrature angular_quadrature( input_reader , quad_fun );
  Time_Data time_data(input_reader);
       
  MMS_Intensity i_mms(input_reader,angular_quadrature);
  try{
    const double xL = cell_data.get_cell_left_edge(0);
    int n_cell = cell_data.get_total_number_of_cells();
    const double xR = cell_data.get_cell_left_edge(n_cell) + cell_data.get_cell_width(n_cell);
    
    std::shared_ptr<V_Transport_BC> left_bc;
    std::shared_ptr<V_Transport_BC> right_bc;
    
    left_bc = std::shared_ptr<V_Transport_BC> (new Transport_BC_MMS(angular_quadrature,input_reader,xL) );
    right_bc = std::shared_ptr<V_Transport_BC> (new Transport_BC_MMS(angular_quadrature,input_reader,xR) );
    
    std::vector<double> time_trials(4,0.);
    time_trials[0] = time_data.get_t_start();
    time_trials[3] = time_data.get_t_end();
    time_trials[1] = time_trials[0] + time_data.get_dt_min() ;
    time_trials[2] = time_trials[1] + time_data.get_dt_min() ;
    
    for(int t = 0 ; t < int(time_trials.size()) ; t++)
    {
      for(int dir=0 ; dir < angular_quadrature.get_number_of_dir() ; dir++)
      {
        double bc_left = i_mms.get_mms_intensity(xL,time_trials[t],dir);
        double bc_right = i_mms.get_mms_intensity(xR,time_trials[t],dir);
        double mu = angular_quadrature.get_mu(dir);
        if(mu > 0.)
        {
          if(fabs(bc_left - left_bc->get_boundary_condition(mu, 0, time_trials[t]) ) > tol )
            throw Dark_Arts_Exception(TRANSPORT,"Left BC not matching up");
        }
        else{
          if(fabs(bc_right - right_bc->get_boundary_condition(mu, 0, time_trials[t]) ) > tol )
            throw Dark_Arts_Exception(TRANSPORT,"right BC not matching up");
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
