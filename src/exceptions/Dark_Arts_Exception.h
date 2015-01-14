#ifndef Dark_Arts_Exception_h
#define Dark_Arts_Exception_h

#include "Exception_Types.h"
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <petsc.h>

class Dark_Arts_Exception
{
public:
  Dark_Arts_Exception( EXCEPTION_TYPE ex_type , const std::string& ex_message);
  Dark_Arts_Exception( EXCEPTION_TYPE ex_type , const std::stringstream& ex_message);
  Dark_Arts_Exception();
  virtual ~Dark_Arts_Exception(){}

  void message(void) const;
  void testing_message(void) const;
  
private:
  std::string error_message;
};


#endif

