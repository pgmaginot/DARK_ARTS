#include <iostream>
#include <vector>
#include <math.h>
#include <string>

#include "Input_Reader.h"
#include "Fem_Quadrature.h"
#include "Quadrule_New.h"
#include "Dark_Arts_Exception.h"

int main(int argc, char** argv)
{
  int val = 0;
  
  /// get all the objects necessary to create an Angular_Quadrature object
  Input_Reader input_reader;  
  try{
    input_reader.read_xml(argv[1]);  
  }
  catch(const Dark_Arts_Exception& da_exception)
  {
    da_exception.testing_message();
    val = -1;
  }
  Quadrule_New quad_fun;   

  Fem_Quadrature fem_quadrature( input_reader , quad_fun);   
 
    
  // Return 0 if tests passed, something else if failing
  return val;
}


