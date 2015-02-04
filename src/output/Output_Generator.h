#ifndef Output_Generator_h
#define Output_Generator_h

#include <iostream>
#include <fstream>
#include <sstream>
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
  Output_Generator(const Angular_Quadrature& angular_quadrature,
    const Fem_Quadrature& fem_quadrature, 
    const Cell_Data& cell_data, 
    const Input_Reader& input_reader);
    
  virtual ~Output_Generator(){}
  

  void write_xml( const bool is_final, const int time_step, const Temperature_Data& temperature);  
  void write_xml( const bool is_final, const int time_step, const Intensity_Data& intensity);
  void write_xml( const bool is_final, const int time_step, const Intensity_Moment_Data& phi);

  void write_txt(const bool is_final, const int time_step, const Temperature_Data& temperature);  
  void write_txt(const bool is_final, const int time_step, const Intensity_Data& intensity);  
  void write_txt(const bool is_final, const int time_step, const Intensity_Moment_Data& phi);  

protected:
  
  const int m_n_dfem;
  const int m_n_cells;
  const int m_n_groups;
  const int m_n_l_mom;
  const int m_n_dir;
  
  std::string m_filename;
  
  const Cell_Data& m_cell_data;
  
  void construct_filename( const int data_type , const bool is_final, const int ts, std::string& output_filename) const;

  void output_cell_data(void);
};


#endif