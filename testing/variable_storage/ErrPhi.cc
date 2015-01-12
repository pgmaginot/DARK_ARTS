#include <iostream>
#include <vector>
#include <math.h>
#include <string>

#include "Input_Reader.h"
#include "Err_Phi.h"
#include "Dark_Arts_Exception.h"

int main(int argc, char** argv)
{
  int val = 0;
  const double tol = 1.0E-6;
  
  Err_Phi test_err;
  try{
    if( abs(test_err.get_cell_with_worst_err() - (-1) ) > 0 )
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Initialized with wrong cell");
    
    if( abs(test_err.get_group_with_worst_err() - (-1) ) > 0 )
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Initialized with wrong group");

    if( abs(test_err.get_legendre_moment_with_worst_err() - (-1) ) > 0 )
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Initialized with wrong moment");

    if( fabs(test_err.get_worst_err() + 1. ) > tol )
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Initialized with wrong error");
  
    const int group = 1;
    const int cell = 3;
    const int moment = 2;
    const double err = 3.75;
    test_err.set_error(cell,group,moment,err);
    
    if( abs(test_err.get_cell_with_worst_err() - cell ) > 0 )
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Saved with wrong cell");
    
    if( abs(test_err.get_group_with_worst_err() - group ) > 0 )
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Saved with wrong group");

    if( abs(test_err.get_legendre_moment_with_worst_err() - moment ) > 0 )
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Saved with wrong moment");

    if( fabs(test_err.get_worst_err() - err ) > tol )
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Saved with wrong error");

    test_err.clear();
      
    if( abs(test_err.get_cell_with_worst_err() - (-1) ) > 0 )
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Cleared with wrong cell");
    
    if( abs(test_err.get_group_with_worst_err() - (-1) ) > 0 )
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Cleared with wrong group");

    if( abs(test_err.get_legendre_moment_with_worst_err() - (-1) ) > 0 )
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Cleared with wrong moment");

    if( fabs(test_err.get_worst_err() +1. ) > tol )
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Cleared with wrong error");
  
  }
  catch(const Dark_Arts_Exception& da_exception )
  {
    da_exception.testing_message();
    val = -1;
  }
  
  // Return 0 if tests passed, somethnig else if failing
  return val;
}
