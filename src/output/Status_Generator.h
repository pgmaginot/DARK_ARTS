#ifndef Status_Generator_h
#define Status_Generator_h

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

#include "Dark_Arts_Exception.h"

/** @file   Status_Generator.h
  *   @author pmaginot
  *   @brief Declare Status_Generator class, that documents iteration progress
  *   @class Status_Generator
 */
class Status_Generator
{
public:
  Status_Generator(std::string input_file);
    
  virtual ~Status_Generator(){ m_status_stream.close(); }
  
  void write_iteration_status(void); 
protected:
  std::string m_input;
  std::string m_xml_ext;
  std::string m_stat_file;
  /// stream that will hold iteration status data
  std::ofstream m_status_stream;
  
};


#endif