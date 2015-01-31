/** @file   Space_Time_Error_Calculator.cc
  *   @author pmaginot
  *   @brief Implement the Space_Time_Error_Calculator class 
 */
         
#include "Space_Time_Error_Calculator.h"
// ##########################################################
// Public functions 
// ##########################################################


Space_Time_Error_Calculator::Space_Time_Error_Calculator(
  const Angular_Quadrature& angular_quadrature,
  const Fem_Quadrature& fem_quadrature, 
  const Cell_Data& cell_data, 
  const Input_Reader& input_reader, const Time_Data& time_data, 
  std::string filename)
  :
  m_last_stage(time_data.get_number_of_stages() - 1),  
  m_space_l2_error_calculator(angular_quadrature,  
    fem_quadrature, 
    cell_data, 
    input_reader),  
  m_time_data(time_data),
  m_total_phi_err(0.),
  m_total_temperature_err(0.),
  m_total_phi_A_err(0.),
  m_total_temperature_A_err(0.),
  m_n_cell(cell_data.get_total_number_of_cells() ),
  m_dfem_ord(fem_quadrature.get_number_of_interpolation_points()  - 1),
  m_dt_max(time_data.get_dt_max() ),
  m_phi_l2_filename( "_space_time_phi_L2_error.txt" ),
  m_phi_A_filename( "_space_time_phi_A_error.txt" ),
  m_temp_l2_filename( "_space_time_tempeature_L2_error.txt" ),
  m_temp_A_filename( "_space_time_tempeature_A_error.txt" ),
  m_output_file_phi_l2( filename+=m_phi_l2_filename , std::ofstream::app),
  m_output_file_phi_A( filename.replace( filename.find(m_phi_l2_filename) , m_phi_l2_filename.length() ,m_phi_A_filename), std::ofstream::app),
  m_output_file_temp_l2( filename.replace( filename.find(m_phi_A_filename) , m_phi_A_filename.length() ,m_temp_l2_filename), std::ofstream::app),
  m_output_file_temp_A( filename.replace( filename.find(m_temp_l2_filename) , m_temp_l2_filename.length() ,m_temp_A_filename), std::ofstream::app),
  m_err_at_each_stage_phi(m_last_stage + 1 , 0.),
  m_err_at_each_stage_temperature(m_last_stage + 1,0.),
  m_err_at_each_stage_phi_A(m_last_stage + 1 , 0.),
  m_err_at_each_stage_temperature_A(m_last_stage + 1,0.)
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

void Space_Time_Error_Calculator::record_error(const double dt, const int stage, const double time_stage, const Intensity_Moment_Data& phi, const Temperature_Data& temperature)
{
  m_err_at_each_stage_phi[stage] = m_space_l2_error_calculator.calculate_l2_error(time_stage,phi);
  m_err_at_each_stage_temperature[stage] = m_space_l2_error_calculator.calculate_l2_error(time_stage,temperature);
  
  m_err_at_each_stage_phi_A[stage] = m_space_l2_error_calculator.calculate_cell_avg_error(time_stage,phi);
  m_err_at_each_stage_temperature_A[stage] = m_space_l2_error_calculator.calculate_cell_avg_error(time_stage,temperature);
  
  std::cout << "New Time step: \n:" ;
  std::cout << "Err1: " << m_err_at_each_stage_phi[stage];
  std::cout << "\nErr2: " << m_err_at_each_stage_temperature[stage];
  std::cout << "\nErr3: " << m_err_at_each_stage_phi_A[stage];
  std::cout << "\nErr4: " << m_err_at_each_stage_temperature_A[stage];
  
  if(stage == m_last_stage)
    update_running_error_totals(dt);
    
  return;
}

void Space_Time_Error_Calculator::update_running_error_totals(const double dt)
{  
  for(int s = 0 ; s <= m_last_stage ; s++)
  {
    m_total_phi_err += m_time_data.get_b(s)*dt*m_err_at_each_stage_phi[s];
    m_total_temperature_err += m_time_data.get_b(s)*dt*m_err_at_each_stage_temperature[s];
    m_total_phi_A_err +=m_time_data.get_b(s)*dt*m_err_at_each_stage_phi_A[s];
    m_total_temperature_A_err += m_time_data.get_b(s)*dt*m_err_at_each_stage_temperature_A[s];
    
    m_err_at_each_stage_phi[s] = 0.0;
    m_err_at_each_stage_temperature[s] = 0.0;
    m_err_at_each_stage_phi_A[s]= 0.0;
    m_err_at_each_stage_temperature_A[s]= 0.0;
  }
  
  
  
  return;
}

void Space_Time_Error_Calculator::print_to_file(void)
{
  m_output_file_phi_l2 << "N_cells: " << m_n_cell << " Fem_ord: " << m_dfem_ord  ;
  m_output_file_phi_l2  << " Dt_max: " <<  std::scientific << std::setprecision(15) << m_dt_max;
  m_output_file_phi_l2 << " Err: " <<  std::scientific << std::setprecision(15) << m_total_phi_err << std::endl;
                       
  m_output_file_temp_l2 << "N_cells: " << m_n_cell << " Fem_ord: " << m_dfem_ord ; 
  m_output_file_temp_l2 << " Dt_max: " << std::scientific << std::setprecision(15) << m_dt_max;
  m_output_file_temp_l2 << " Err: " << std::scientific << std::setprecision(15) << m_total_temperature_err << std::endl;
                       
  m_output_file_phi_A  << "N_cells: " << m_n_cell << " Fem_ord: " << m_dfem_ord  ;
  m_output_file_phi_A  << " Dt_max: " << std::scientific << std::setprecision(15) << m_dt_max;
  m_output_file_phi_A  << " Err: " << std::scientific << std::setprecision(15) << m_total_phi_A_err << std::endl;
                       
  m_output_file_temp_A  << "N_cells: " << m_n_cell << " Fem_ord: " << m_dfem_ord  ;
  m_output_file_temp_A  << " Dt_max: " <<  std::scientific << std::setprecision(15) << m_dt_max;
  m_output_file_temp_A  << " Err: " << std::scientific << std::setprecision(15) << m_total_temperature_A_err << std::endl;
                        
  return;
}

