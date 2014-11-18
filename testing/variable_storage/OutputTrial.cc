#include <iostream>
#include <vector>
#include <math.h>
#include <string>

#include "Input_Reader.h"
#include "Cell_Data.h"
#include "Dark_Arts_Exception.h"
#include "Output_Generator.h"

int main(int argc, char** argv)
{
  int val = 0;
    
  Input_Reader input_reader;    
    
  try{
    input_reader.read_xml(argv[1]);
  }
  catch( const Dark_Arts_Exception& da_exception )
  {
    da_exception.testing_message();
    val = -1;
  }
  
  std::string my_str = "Does nothing";
  Output_Generator output(input_reader, my_str);
  output.create_output(false, 0. ,1);
  
  // Return 0 if tests passed, somethnig else if failing
  return val;
}
