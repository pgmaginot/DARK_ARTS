#ifndef Output_Generator_h
#define Output_Generator_h

#include "tinyxml.h"
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include "Input_Reader.h"
#include "Temperature_Data.h"
#include "Intensity_Data.h"
#include "Intensity_Moment_Data.h"
#include "Dark_Arts_Exception.h"

/** @file   Output_Generator.h
  *   @author pmaginot
  *   @brief Declare Input_Reader class
  *   @class Input_Reader
 */
class Output_Generator
{
public:
  Output_Generator(const Fem_Quadrature& fem_quadrature, const Cell_Data& cell_data, std::string input_filename);
  virtual ~Output_Generator(){}
  
  //! read the supplied input file, start populating data objects
  void create_output(const bool is_final, const double time, const int time_step);  
  
protected:
  void write_xml( std::string xmlFilename , const Temperature_Data& temperature);
  void write_xml( std::string xmlFilename , const Intensity_Moment_Data& phi);
  void write_xml( std::string xmlFilename , const Intensity_Data& intensity);

  const int m_n_dfem;
  const int m_n_cells;
  
  std::string m_filename_base;
};


#endif