#ifndef Status_Generator_h
#define Status_Generator_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
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
    
  virtual ~Status_Generator(){ write_final_counts(); m_status_stream.close(); }
  
  void write_iteration_status(const int step, const int stage, const int therm_iter , const double dt , const int inners , const double err, const double damp) ;
  
  int get_total_thermals(void) const {return m_total_thermal_iterations;}
  int get_total_sweeps(void) const {return m_total_inner_solves;}
protected:
  void write_final_counts(void);
  
  std::string m_input;
  std::string m_xml_ext;
  std::string m_stat_file;
  /// stream that will hold iteration status data
  std::ofstream m_status_stream;
  
  int m_total_thermal_iterations;
  int m_total_inner_solves;
};


#endif