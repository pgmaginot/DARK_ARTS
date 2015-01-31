#ifndef Space_Time_Error_Calculator_h
#define Space_Time_Error_Calculator_h

#include <fstream>
#include "Input_Reader.h"
#include "Angular_Quadrature.h"
#include "Cell_Data.h"
#include "Temperature_Data.h"
#include "Intensity_Data.h"
#include "Intensity_Moment_Data.h"
#include "Fem_Quadrature.h"
#include "Time_Data.h"
#include "L2_Error_Calculator.h"
#include "Dark_Arts_Exception.h"

/** @file   Space_Time_Error_Calculator.h
  *   @author pmaginot
  *   @brief Declare Space_Time_Error_Calculator class
  *   @class Space_Time_Error_Calculator
 */
class Space_Time_Error_Calculator
{
public:
  Space_Time_Error_Calculator(const Angular_Quadrature& angular_quadrature,
    const Fem_Quadrature& fem_quadrature, 
    const Cell_Data& cell_data, 
    const Input_Reader& input_reader, 
    const Time_Data& time_data, 
    std::string filename_base);
    
  virtual ~Space_Time_Error_Calculator(){
    print_to_file();
    
    m_output_file_phi_l2.close();
    m_output_file_phi_A.close();
    m_output_file_temp_l2.close();
    m_output_file_temp_A.close();
  }

  void record_error(const double dt, const int stage, const double time_stage, const Intensity_Moment_Data& phi, const Temperature_Data& temperature);

  
protected:  
  void print_to_file(void);
  
  void update_running_error_totals(const double dt);
  
  const int m_last_stage;
  
  L2_Error_Calculator m_space_l2_error_calculator;
  
  const Time_Data& m_time_data;
  
  double m_total_phi_err;
  double m_total_temperature_err;
  double m_total_phi_A_err;
  double m_total_temperature_A_err;
  
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
  
  std::vector<double> m_err_at_each_stage_phi;
  std::vector<double> m_err_at_each_stage_temperature;
  std::vector<double> m_err_at_each_stage_phi_A;
  std::vector<double> m_err_at_each_stage_temperature_A;
  
};


#endif