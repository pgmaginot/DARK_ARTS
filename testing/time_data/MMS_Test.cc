#include <iostream>
#include <vector>
#include <math.h>
#include <string>

#include "Input_Reader.h"
#include "Time_Data.h"
#include "Dark_Arts_Exception.h"
#include "MMS_Temperature.h"

int main(int argc, char** argv)
{
  int val = 0;
  
  Input_Reader input_reader;  
  try{
    input_reader.read_xml(argv[1]);  
  }
  catch(const Dark_Arts_Exception& da_exception)
  {
    da_exception.testing_message();
    val = -1;
  }
  
  Time_Data time_data(input_reader);
  MMS_Temperature temperature( input_reader );  

  /// checking SDIRK_2_2 (-)
  const int n_stage = 2;
  
  double dt = time_data.get_dt_max();
  double time = time_data.get_t_start();
  const double time_end = time_data.get_t_end();
  
  const double x = 2.0;
  double time_stage = 0.;
  
  double temperature_eval = temperature.get_mms_temperature(x,time_stage);
  std::vector<double> k_t(n_stage,0.);
  
  const double tol = 1E-12;
      
  try{
    for( ; time < time_end ; )
    {
      if( (time+dt) > time_end)
        dt = time_end - time;
        
      for(int s = 0; s < n_stage ; s++)
      {
        time_stage = time + time_data.get_c(s)*dt;
        k_t[s] = temperature.get_mms_temperature_time_derivative(x,time_stage);
        temperature_eval += dt*time_data.get_b(s)*k_t[s];
      } 
      time += dt;
    }
    
    double t_expected = temperature.get_mms_temperature(x,time_end);
    
    std::cout << "Temperature_new: " << temperature_eval << " Expected temperature: " << t_expected << std::endl;
    std::cout << "Normalized Temperature difference: " << fabs( (t_expected - temperature_eval)) << std::endl;
    
    if(fabs( temperature_eval - t_expected) > tol)
      throw Dark_Arts_Exception(TIME_MARCHER , "Debugging SDIRK");
  
  }
  catch(const Dark_Arts_Exception& da)
  {
    da.testing_message();
    val = -1;
  }
      
  // Return 0 if tests passed, something else if failing
  return val;
}
