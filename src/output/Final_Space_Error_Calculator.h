#ifndef Final_Space_Error_Calculator_h
#define Final_Space_Error_Calculator_h

#include <fstream>
#include "Input_Reader.h"
#include "Temperature_Data.h"
#include "Intensity_Data.h"
#include "Intensity_Moment_Data.h"
#include "Dark_Arts_Exception.h"
#include "Cell_Data.h"
#include "Angular_Quadrature.h"
#include "Fem_Quadrature.h"
#include "Time_Data.h"
#include "L2_Error_Calculator.h"

/** @file   Final_Space_Error_Calculator.h
  *   @author pmaginot
  *   @brief Declare Final_Space_Error_Calculator class
  *   @class Final_Space_Error_Calculator
 */
class Final_Space_Error_Calculator
{
public:
  Final_Space_Error_Calculator(const Angular_Quadrature& angular_quadrature, const Fem_Quadrature& fem_quadrature,
  const Cell_Data& cell_data, const Input_Reader& input_reader, const Time_Data& time_data, std::string filename_base);
    
  virtual ~Final_Space_Error_Calculator(){}
  
  void record_error(const double time_final , const Temperature_Data& temperature, const Intensity_Moment_Data& phi);

protected:

  L2_Error_Calculator m_space_l2_error_calculator;
  
  double m_phi_l2_err;
  double m_temperature_l2_err;
  double m_phi_A_err;
  double m_temperature_A_err;
  
  const int m_n_cell;
  const int m_dfem_ord;
  const double m_dt_max;
  
  std::string m_phi_l2_filename;
  std::string m_phi_A_filename;
  std::string m_temp_l2_filename;
  std::string m_temp_A_filename;
  
  std::ofstream m_output_file_phi_l2;
  std::ofstream m_output_file_phi_A;
  std::ofstream m_output_file_temp_l2;
  std::ofstream m_output_file_temp_A;

};


#endif