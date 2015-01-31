#include <iostream>
#include <vector>
#include <math.h>
#include <string>

#include "Input_Reader.h"
#include "Time_Data.h"
#include "Dark_Arts_Exception.h"

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

  /// checking implicit euler, t_start = 0, t_end = 5 , dt_min = 0.1, dt_max = 1, Ramp starter using 5 steps
  const double t_start = 0.;
  const double t_end = 5.;
  const double dt_min = 0.1;
  const double dt_max = 1.1;
  const int n_stage = 1;
  
  const double a = 1.0;
  const double c = 1.0;
  const double b = 1.0;
  
  const double tol = 1.0E-6;
  
  /// for testing fake time loop
  const int too_many_steps = 10;
  std::vector<double> expected_dt(too_many_steps,0.);
  
  expected_dt[0] = 0.1;
  expected_dt[1] = 0.3;
  expected_dt[2] = 0.5;
  expected_dt[3] = 0.7;
  expected_dt[4] = 0.9;
  expected_dt[5] = 1.1;
  expected_dt[6] = 1.1;
  expected_dt[7] = 0.3;
  expected_dt[8] = 0.;
  expected_dt[9] = 0.;
  
  try{
    if( time_data.get_number_of_stages() != n_stage)
      throw Dark_Arts_Exception(SUPPORT_OBJECT, "Not getting the correct number of stages for implicit Euler");    
  }
  catch(const Dark_Arts_Exception& da)
  {
    da.testing_message();
    val = -1;
  }
  
  try{
    if( fabs(time_data.get_t_start() - t_start) > tol)
      throw Dark_Arts_Exception(SUPPORT_OBJECT, "Not getting the correct start time");    
  }
  catch(const Dark_Arts_Exception& da)
  {
    da.testing_message();
    val = -1;
  }
  
  try{
    if( fabs(time_data.get_t_end() - t_end) > tol)
      throw Dark_Arts_Exception(SUPPORT_OBJECT, "Not getting the correct end time");    
  }
  catch(const Dark_Arts_Exception& da)
  {
    da.testing_message();
    val = -1;
  }
  
  try{
    if( fabs(time_data.get_dt_min() - dt_min) > tol)
      throw Dark_Arts_Exception(SUPPORT_OBJECT, "Not getting the correct dt_min");    
  }
  catch(const Dark_Arts_Exception& da)
  {
    da.testing_message();
    val = -1;
  }
  
  try{
    if( fabs(time_data.get_dt_max() - dt_max) > tol)
      throw Dark_Arts_Exception(SUPPORT_OBJECT, "Not getting the correct dt_max");    
  }
  catch(const Dark_Arts_Exception& da)
  {
    da.testing_message();
    val = -1;
  }
  
  try{
    if( fabs(time_data.get_a(0 , 0) - a) > tol)
      throw Dark_Arts_Exception(SUPPORT_OBJECT , "Not the correct 'a' for SDIRK Implicit Euler");
  }
  catch(const Dark_Arts_Exception& da)
  {
    da.testing_message();
    val = -1;
  }
  
  try{
    if( fabs(time_data.get_b(0 ) - b) > tol)
      throw Dark_Arts_Exception(SUPPORT_OBJECT , "Not the correct 'b' for SDIRK Implicit Euler");
  }
  catch(const Dark_Arts_Exception& da)
  {
    da.testing_message();
    val = -1;
  }
  
  try{
    if( fabs(time_data.get_c(0 ) - c) > tol)
      throw Dark_Arts_Exception(SUPPORT_OBJECT , "Not the correct 'c' for SDIRK Implicit Euler");
  }
  catch(const Dark_Arts_Exception& da)
  {
    da.testing_message();
    val = -1;
  }
  
  /// test in a fake time loop
  double time_now = t_start;
  double dt;
  try{
    for(int ts = 0; ts < too_many_steps; ts++)
    {
      dt = time_data.get_dt(ts,time_now);
      
      if( fabs(dt - expected_dt[ts]) > tol)
        throw Dark_Arts_Exception(TIME_MARCHER , "Not calculating correct ramp dt");
      
      std::cout << "Time now: " << time_now << " Expected dt: " << expected_dt[ts] << " Calculated dt: " << dt << std::endl;
      time_now += dt;
    }
  }
  catch(const Dark_Arts_Exception& da)
  {
    da.testing_message();
    val = -1;
  }
    
  // Return 0 if tests passed, something else if failing
  return val;
}
