/** @file   Final_Space_Error_Calculator.cc
  *   @author pmaginot
  *   @brief Implement the Final_Space_Error_Calculator class 
 */
         
#include "Final_Space_Error_Calculator.h"
// ##########################################################
// Public functions 
// ##########################################################


Final_Space_Error_Calculator::Final_Space_Error_Calculator(const Angular_Quadrature& angular_quadrature, const Fem_Quadrature& fem_quadrature,
  const Cell_Data& cell_data, const Input_Reader& input_reader, const Time_Data& time_data, std::string filename_base):
  m_space_l2_error_calculator(angular_quadrature,  
    fem_quadrature, 
    cell_data, 
    input_reader),  
  m_phi_l2_err(0.),
  m_temperature_l2_err(0.),
  m_phi_A_err(0.),
  m_temperature_A_err(0.),
  m_n_cell(cell_data.get_total_number_of_cells() ),
  m_dfem_ord(fem_quadrature.get_number_of_interpolation_points()  - 1),
  m_dt_max(time_data.get_dt_max() ),
  m_phi_l2_filename( "_final_space_phi_L2_error.txt" ),
  m_phi_A_filename( "_final_space_phi_A_error.txt" ),
  m_temp_l2_filename( "_final_space_tempeature_L2_error.txt" ),
  m_temp_A_filename( "_final_space_tempeature_A_error.txt" ),
  m_output_file_phi_l2( filename_base+=m_phi_l2_filename , std::ofstream::app),
  m_output_file_phi_A( filename_base.replace( filename_base.find(m_phi_l2_filename) , m_phi_l2_filename.length() ,m_phi_A_filename), std::ofstream::app),
  m_output_file_temp_l2( filename_base.replace( filename_base.find(m_phi_A_filename) , m_phi_A_filename.length() ,m_temp_l2_filename), std::ofstream::app),
  m_output_file_temp_A( filename_base.replace( filename_base.find(m_temp_l2_filename) , m_temp_l2_filename.length() ,m_temp_A_filename), std::ofstream::app)
{
  std::stringstream err;  

  /// check the status of different filestreams
  if(!m_output_file_phi_l2.is_open() )
  {
    err << "m_output_file_phi_l2 is not open.  Tried to open a file named: " << m_phi_l2_filename << std::endl;
    throw Dark_Arts_Exception(INPUT , err);
  }
  
  if(!m_output_file_phi_A.is_open() )
  {
    err << "m_output_file_phi_A is not open.  Tried to open a file named: " << m_phi_A_filename << std::endl;
    throw Dark_Arts_Exception(INPUT , err);
  }
  
  if(!m_output_file_temp_l2.is_open() )
  {
    err << "m_output_file_temp_l2 is not open.  Tried to open a file named: " << m_temp_l2_filename << std::endl;
    throw Dark_Arts_Exception(INPUT , err);
  }
  
  if(!m_output_file_temp_A.is_open() )
  {
    err << "m_output_file_temp_A is not open.  Tried to open a file named: " << m_temp_A_filename << std::endl;
    throw Dark_Arts_Exception(INPUT , err);
  }
}

void Final_Space_Error_Calculator::record_error(
  const double time_final , const Temperature_Data& temperature, const Intensity_Moment_Data& phi)
{  
  m_phi_l2_err = m_space_l2_error_calculator.calculate_l2_error(time_final,phi);
  m_temperature_l2_err = m_space_l2_error_calculator.calculate_l2_error(time_final,temperature);
  
  m_phi_A_err = m_space_l2_error_calculator.calculate_cell_avg_error(time_final,phi);
  m_temperature_A_err = m_space_l2_error_calculator.calculate_cell_avg_error(time_final,temperature);
  
  m_output_file_phi_l2 << "N_cells: " << m_n_cell << " Fem_ord: " << m_dfem_ord  ;
  m_output_file_phi_l2  << " Dt_max: " <<  std::scientific << std::setprecision(15) << m_dt_max;
  m_output_file_phi_l2 << " Err: " <<  std::scientific << std::setprecision(15) << m_phi_l2_err << std::endl;
                       
  m_output_file_temp_l2 << "N_cells: " << m_n_cell << " Fem_ord: " << m_dfem_ord ; 
  m_output_file_temp_l2 << " Dt_max: " << std::scientific << std::setprecision(15) << m_dt_max;
  m_output_file_temp_l2 << " Err: " << std::scientific << std::setprecision(15) << m_temperature_l2_err << std::endl;
                       
  m_output_file_phi_A  << "N_cells: " << m_n_cell << " Fem_ord: " << m_dfem_ord  ;
  m_output_file_phi_A  << " Dt_max: " << std::scientific << std::setprecision(15) << m_dt_max;
  m_output_file_phi_A  << " Err: " << std::scientific << std::setprecision(15) << m_phi_A_err << std::endl;
                       
  m_output_file_temp_A  << "N_cells: " << m_n_cell << " Fem_ord: " << m_dfem_ord  ;
  m_output_file_temp_A  << " Dt_max: " <<  std::scientific << std::setprecision(15) << m_dt_max;
  m_output_file_temp_A  << " Err: " << std::scientific << std::setprecision(15) << m_temperature_A_err << std::endl;
  
  return;
}