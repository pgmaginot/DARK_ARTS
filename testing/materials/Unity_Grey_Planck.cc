#include <iostream>
#include <vector>
#include <math.h>
#include <string>

#include "Input_Reader.h"
#include "Planck.h"
#include "Dark_Arts_Exception.h"

/// check that a=c=1
/// check that for grey planck is acT^4 (per steradian) and derivative term
int main(int argc, char** argv)
{
  int val = 0;
  
  const double t = 3.;
  const double sn_w = 2.;
  
  Input_Reader input_reader;  
  try{
    input_reader.read_xml(argv[1]);  
  }
  catch(const Dark_Arts_Exception& da_exception)
  {
    da_exception.testing_message();
    val = -1;
  }
  
  Planck planck(1.0E-15, input_reader, sn_w);
  /** test the following functions:
      \fn double get_c(void);
      \fn double integrate_B_grey(double T);
      \fn double integrate_dBdT_grey(double T);
  */
  try{
    if( fabs(planck.get_c() - 1.) > 1.0E-4 )
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Not seeting a unity c");
  }
  catch(const Dark_Arts_Exception& da_exception)
  {
    da_exception.testing_message();
    val = -1;
  }
  
  try{
    if( fabs( planck.integrate_B_grey(t) - pow(t,4)/sn_w ) > 1.0E-4 )
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Not correctly calculating unity planck");
  }
  catch(const Dark_Arts_Exception& da_exception)
  {
    da_exception.testing_message();
    val = -1;
  }
  
  try{
    if( fabs(planck.integrate_dBdT_grey(t) - 4.*pow(t,3)/sn_w ) > 1.0E-4 )
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Not correctly calculating grey unity planck derivative");
  }
  catch(const Dark_Arts_Exception& da_exception)
  {
    da_exception.testing_message();
    val = -1;
  }
  
  // Return 0 if tests passed, something else if failing
  return val;
}
