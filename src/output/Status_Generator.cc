/** @file   Input_Reader.cc
  *   @author pmaginot
  *   @brief Implement the Input_Reader class that reads XML files with TinyXml
 */
         
#include "Status_Generator.h"
// ##########################################################
// Public functions 
// ##########################################################


Status_Generator::Status_Generator(std::string input_file):
  m_input(input_file) , 
  m_xml_ext(".xml"),
  m_stat_file( input_file.replace( m_input.find(m_xml_ext) , m_xml_ext.length() , "_iteration_status.txt") ),
  m_status_stream( m_stat_file ),
  m_total_thermal_iterations(0),
  m_total_inner_solves(0)
{
  std::cout << "This is the stat file: " << m_stat_file << std::endl;
  if(!m_status_stream)
    throw Dark_Arts_Exception(SUPPORT_OBJECT, "Could not create iteration Status_Generator filestream");

}

void Status_Generator::write_iteration_status(const int step, const int stage, const int therm_iter , 
  const double dt , const int inners , const double err, const double damp) 
{
  m_total_thermal_iterations++;
  m_total_inner_solves += inners;
  
  m_status_stream <<  "Time_step: " << std::setw(6) << step  << 
                      " Stage: " << std::setw(2) << stage <<     
                      " Thermal_iter: " << std::setw(4) << therm_iter << 
                      " dt: " ;
  m_status_stream << std::scientific << std::setprecision(15) << dt ;
  m_status_stream << " number_of_inner_solves: " << inners ;
  m_status_stream << " Relative_Error: " << std::scientific << std::setprecision(15) << err ;
  m_status_stream << " Damping: " << damp << std::endl;
  return;
}  

void Status_Generator::write_final_counts(void)
{
  m_status_stream << "\n\n Total_thermal_iters: " << m_total_thermal_iterations <<
                    "\n Total_Transport_sweeps: " << m_total_inner_solves << std::endl;
  return;
}

