#include "Dark_Arts_Exception.h"

Dark_Arts_Exception::Dark_Arts_Exception( EXCEPTION_TYPE ex_type , const std::string& ex_message )
{
  switch(ex_type)
  {
    case INPUT:
    {
      error_message = "INPUT: ";
      break;
    }
    case FEM:
    {
      error_message = "FEM: ";
      break;
    }
    case VARIABLE_STORAGE:
    {
      error_message = "VARIABLE_STORAGE: ";
      break;
    }
    case SUPPORT_OBJECT:
    {
      error_message = "SUPPORT_OBJECT: ";
      break;
    }
    case TIME_MARCHER:
    {
      error_message = "TIME_MARCHER: ";
      break;
    }
    case TRANSPORT:
    {
      error_message = "TRANSPORT: ";
      break;
    }
    case MIP:
    {
      error_message = "MIP: ";
      break;
    }
    default:
    {
      error_message = "Unknown Exception Type: " ;
      break;
    }
  }
  
  error_message.append( ex_message );      
}

Dark_Arts_Exception::Dark_Arts_Exception( EXCEPTION_TYPE ex_type , const std::stringstream& ex_message )
{
  switch(ex_type)
  {
    case INPUT:
    {
      error_message = "INPUT: ";
      break;
    }
    case FEM:
    {
      error_message = "FEM: ";
      break;
    }
    case VARIABLE_STORAGE:
    {
      error_message = "VARIABLE_STORAGE: ";
      break;
    }
    case SUPPORT_OBJECT:
    {
      error_message = "SUPPORT_OBJECT: ";
      break;
    }
    case TIME_MARCHER:
    {
      error_message = "TIME_MARCHER: ";
      break;
    }
    case TRANSPORT:
    {
      error_message = "TRANSPORT: ";
      break;
    }
    case MIP:
    {
      error_message = "MIP: ";
      break;
    }
    default:
    {
      error_message = "Unknown Exception Type: " ;
      break;
    }
  }
  
  error_message.append( ex_message.str() );      
}

Dark_Arts_Exception::Dark_Arts_Exception(  )
{
  error_message = "Unknown Exception Type: " ;  
  error_message.append("Dark_Arts_Exception thrown without any initializes" ) ;     
}

void Dark_Arts_Exception::message(void) const
{
  std::cerr << error_message << std::endl;
  PetscFinalize();
  exit(EXIT_FAILURE);
  return;
}

void Dark_Arts_Exception::testing_message(void) const
{
  std::cout << "At least one test failed! " << std::endl;
  std::cout << error_message << std::endl;
  return;
}
  


