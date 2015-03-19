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
  m_wg_tolerance(input_reader.get_within_group_solve_tolerance() ),
  m_bg_tolerance(input_reader.get_between_group_solve_tolerance() ), 
  m_thermal_tolerance(input_reader.get_thermal_tolerance() ),
  m_phi_l2_filename( "_final_space_phi_L2_error.txt" ),
  m_phi_A_filename( "_final_space_phi_A_error.txt" ),
  m_temp_l2_filename( "_final_space_temperature_L2_error.txt" ),
  m_temp_A_filename( "_final_space_temperature_A_error.txt" ),
  m_output_file_phi_l2( filename_base += m_phi_l2_filename , std::ofstream::app),
  m_output_file_phi_A( filename_base.replace( filename_base.find(m_phi_l2_filename) , m_phi_l2_filename.length() ,m_phi_A_filename), std::ofstream::app),
  m_output_file_temp_l2( filename_base.replace( filename_base.find(m_phi_A_filename) , m_phi_A_filename.length() ,m_temp_l2_filename), std::ofstream::app),
  m_output_file_temp_A( filename_base.replace( filename_base.find(m_temp_l2_filename) , m_temp_l2_filename.length() ,m_temp_A_filename), std::ofstream::app)
{
  TIME_SOLVER time_solver = input_reader.get_time_solver();
  
  if(time_solver == IMPLICIT_EULER)
  {
    m_time_str = "EULER";
  }
  else if(time_solver == ALEXANDER_2_2)
  {
    m_time_str = "2_2";
  }
  else if(time_solver == ALEXANDER_2_2_PLUS)
  {
    m_time_str = "2_2_PLUS";
  }
  else if(time_solver == ALEXANDER_3_3)
  {
    m_time_str = "3_3";
  }

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
  const double time_final , const int n_steps, const Temperature_Data& temperature, const Intensity_Moment_Data& phi
  , const int n_thermals, const int n_sweeps)
{  
  m_phi_l2_err = m_space_l2_error_calculator.calculate_l2_error(time_final,phi);
  m_temperature_l2_err = m_space_l2_error_calculator.calculate_l2_error(time_final,temperature);
  
  m_phi_A_err = m_space_l2_error_calculator.calculate_cell_avg_error(time_final,phi);
  m_temperature_A_err = m_space_l2_error_calculator.calculate_cell_avg_error(time_final,temperature);
  
  m_output_file_phi_l2 << "N_cells: " << std::setw(4) << m_n_cell << " Fem_ord: " << std::setw(2) << m_dfem_ord  << " Time_str: " << m_time_str;
  m_output_file_phi_l2 << " N_steps: " << std::setw(6) << n_steps ;
  m_output_file_phi_l2  << " Dt_max: " <<  std::scientific << std::setprecision(15) << m_dt_max;
  m_output_file_phi_l2  << " WG_tol: " << std::scientific << std::setprecision(5) << m_wg_tolerance;
  m_output_file_phi_l2 << " BG_tol: " << std::scientific << std::setprecision(5) << m_bg_tolerance;
  m_output_file_phi_l2 << " Thermal_tol: " << std::scientific << std::setprecision(5) << m_thermal_tolerance;
  m_output_file_phi_l2 << " Err: " <<  std::scientific << std::setprecision(15) << m_phi_l2_err 
                        << " Total_thermals: " << n_thermals << " Total_sweeps: " << n_sweeps << std::endl;
                       
  m_output_file_temp_l2 << "N_cells: " << std::setw(4) << m_n_cell << " Fem_ord: " << std::setw(2) << m_dfem_ord  << " Time_str: " << m_time_str ;
  m_output_file_temp_l2 << " N_steps: " << std::setw(6) << n_steps ;
  m_output_file_temp_l2 << " Dt_max: " << std::scientific << std::setprecision(15) << m_dt_max;
  m_output_file_temp_l2  << " WG_tol: " << std::scientific << std::setprecision(5) << m_wg_tolerance ;
  m_output_file_temp_l2 << " BG_tol: " << std::scientific << std::setprecision(5) << m_bg_tolerance;
  m_output_file_temp_l2 << " Thermal_tol: " << std::scientific << std::setprecision(5) << m_thermal_tolerance;
  m_output_file_temp_l2 << " Err: " << std::scientific << std::setprecision(15) << m_temperature_l2_err 
                        << " Total_thermals: " << n_thermals << " Total_sweeps: " << n_sweeps << std::endl;
                       
  m_output_file_phi_A << "N_cells: " << std::setw(4) << m_n_cell << " Fem_ord: " << std::setw(2) << m_dfem_ord  << " Time_str: " << m_time_str ;
  m_output_file_phi_A << " N_steps: " << std::setw(6) << n_steps ;
  m_output_file_phi_A  << " Dt_max: " << std::scientific << std::setprecision(15) << m_dt_max;
  m_output_file_phi_A  << " WG_tol: " << std::scientific << std::setprecision(5) << m_wg_tolerance;
  m_output_file_phi_A << " BG_tol: " << std::scientific << std::setprecision(5) << m_bg_tolerance ;
  m_output_file_phi_A << " Thermal_tol: " << std::scientific << std::setprecision(5) << m_thermal_tolerance;
  m_output_file_phi_A  << " Err: " << std::scientific << std::setprecision(15) << m_phi_A_err 
                        << " Total_thermals: " << n_thermals << " Total_sweeps: " << n_sweeps << std::endl;
                       
  m_output_file_temp_A << "N_cells: " << std::setw(4) << m_n_cell << " Fem_ord: " << std::setw(2) << m_dfem_ord  << " Time_str: " << m_time_str ;
  m_output_file_temp_A << " N_steps: " << std::setw(6) << n_steps ;
  m_output_file_temp_A  << " Dt_max: " <<  std::scientific << std::setprecision(15) << m_dt_max;
  m_output_file_temp_A  << " WG_tol: " << std::scientific << std::setprecision(5) << m_wg_tolerance;
  m_output_file_temp_A << " BG_tol: " << std::scientific << std::setprecision(5) << m_bg_tolerance;
  m_output_file_temp_A << " Thermal_tol: " << std::scientific << std::setprecision(5) << m_thermal_tolerance ;
  m_output_file_temp_A  << " Err: " << std::scientific << std::setprecision(15) << m_temperature_A_err
                        << " Total_thermals: " << n_thermals << " Total_sweeps: " << n_sweeps << std::endl;
  
  return;
}