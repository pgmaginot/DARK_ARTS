#include <iostream>
#include <vector>
#include <math.h>
#include <string>

#include "Input_Reader.h"
#include "Err_Temperature.h"
#include "Dark_Arts_Exception.h"

int main(int argc, char** argv)
{
  int val = 0;
  const double tol = 1.0E-6;
  
  const int n_p = 3;
  Err_Temperature test_err(n_p);
  
  Eigen::VectorXd dummy_vec;
  
  
    
  try{
    if( abs(test_err.get_cell_with_worst_err() - (-1) ) > 0 )
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Initialized with wrong cell");
    
    if( abs(test_err.get_element_with_worst_err() - (-1) ) > 0 )
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Initialized with wrong element");
    
    if( fabs(test_err.get_worst_err() ) > tol )
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Initialized with wrong error");
  
    test_err.get_big_delta(dummy_vec);
    for(int el = 0; el < n_p ; el++)
    {
      if( fabs(dummy_vec(el) ) > tol)
        throw Dark_Arts_Exception(SUPPORT_OBJECT , "Wrong big delta initialized");
    }
      
    const int cell = 3;
    const int element = 2;
    const double err = 2.75;
    dummy_vec = Eigen::VectorXd::Zero(n_p);
    dummy_vec(0) = -1.2;
    dummy_vec(1) = 2.1;
    dummy_vec(2) = 0.5;
    
    test_err.set_error(cell,element,err,dummy_vec);
    
    if( abs(test_err.get_cell_with_worst_err() - cell ) > 0 )
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Saved with wrong cell");
    
    if( abs(test_err.get_element_with_worst_err() - element ) > 0 )
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Saved with wrong element");

    if( fabs(test_err.get_worst_err() - err ) > tol )
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Saved with wrong error");

    Eigen::VectorXd delta_vec;
    test_err.get_big_delta(delta_vec);
    for(int el = 0; el < n_p ; el++)
    {
      if( fabs(delta_vec(el) - dummy_vec(el) ) > tol)
        throw Dark_Arts_Exception(SUPPORT_OBJECT , "Wrong big delta returned");
    } 
    
    test_err.clear();
      
    if( abs(test_err.get_cell_with_worst_err() - (-1) ) > 0 )
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Cleared with wrong cell");
    
    if( abs(test_err.get_element_with_worst_err() - (-1) ) > 0 )
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Cleared with wrong element");
    
    if( fabs(test_err.get_worst_err() ) > tol )
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Cleared with wrong error");
  
    test_err.get_big_delta(dummy_vec);
    for(int el = 0; el < n_p ; el++)
    {
      if( fabs(dummy_vec(el) ) > tol)
        throw Dark_Arts_Exception(SUPPORT_OBJECT , "Wrong big delta cleared");
    }

    const double small_number = 3.0E-8;
    test_err.set_small_number(small_number);
    if( fabs( (test_err.get_small_number() - small_number)/small_number ) > tol)
      throw Dark_Arts_Exception(SUPPORT_OBJECT , "Not setting/getting Err_Temperature small_number_correctly");
    
  }
  catch(const Dark_Arts_Exception& da_exception )
  {
    da_exception.testing_message();
    val = -1;
  }
  
  // Return 0 if tests passed, somethnig else if failing
  return val;
}
