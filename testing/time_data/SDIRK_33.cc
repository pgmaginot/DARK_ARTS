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

  /// checking SDIRK_3_3 (-)
  const int n_stage = 3;
  
  std::vector<double> a_expected(n_stage*(n_stage+1)/2,0.);
  std::vector<double> b_expected(n_stage,0.);
  std::vector<double> c_expected(n_stage,0.);
  
  
  double alpha = 0.435866521508459;
  
  const double tau_2 = (1. + alpha)/2.;
  
  b_expected[0] = -(6.*alpha*alpha - 16.*alpha + 1.)/4.;
  b_expected[1] = (6.*alpha*alpha - 20.*alpha + 5.)/4.;
  b_expected[2] = alpha;
  
  c_expected[0] = alpha;
  c_expected[1] = tau_2;
  c_expected[2] = 1.;
  
  a_expected[0] = alpha;
  
  a_expected[1] = tau_2 - alpha;
  a_expected[2] = alpha;
  
  a_expected[3] = b_expected[0];
  a_expected[4] = b_expected[1];
  a_expected[5] = alpha;
  
  const double tol = 1.0E-8;
      
  try{
    int calc_n_stages = time_data.get_number_of_stages();
    std::cout << "Dark_Arts thinks Alexander_2_2 has: " << calc_n_stages<< "stages" << std::endl;
    if( calc_n_stages != n_stage)
      throw Dark_Arts_Exception(SUPPORT_OBJECT, "Not getting the correct number of stages for implicit Euler");    

    int cnt = 0; 
    for(int stage=0 ; stage < n_stage ; stage++)
    {
      std::cout << "b["<<stage<<"] expected: " << b_expected[stage] << " received: " << time_data.get_b(stage) <<std::endl;
      if( fabs(time_data.get_b(stage) - b_expected[stage]) > tol)
        throw Dark_Arts_Exception(SUPPORT_OBJECT , "Not the correct 'b' for Alexander_2_2");

      std::cout << "c["<<stage<<"] expected: " << c_expected[stage] << " received: " << time_data.get_c(stage) <<std::endl;
      if( fabs(time_data.get_c(stage) - c_expected[stage]) > tol)
        throw Dark_Arts_Exception(SUPPORT_OBJECT , "Not the correct 'c' for Alexander_2_2");
        
      for(int j = 0 ; j <= stage ; j++)
      {
        std::cout << "Count: " << cnt << std::endl;
        std::cout << "a["<<stage << "," << j <<"] expected: " << a_expected[cnt] << " received: " << time_data.get_a(stage,j) <<std::endl;
        if( fabs(time_data.get_a(stage , j) - a_expected[cnt] ) > tol)
          throw Dark_Arts_Exception(SUPPORT_OBJECT , "Not the correct 'a' for Alexander 2_2");
          
        cnt++;
      }
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
