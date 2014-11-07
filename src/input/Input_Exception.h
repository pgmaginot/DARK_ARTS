#ifndef Input_Exceptions_h
#define Input_Exceptions_h

#include <exception>

class Input_Exception : public std::exception
{

  const char * what () const throw()
  {
    return "Input Error_Message: ";
  }
};

#endif