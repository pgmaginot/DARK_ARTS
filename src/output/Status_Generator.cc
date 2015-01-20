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
  m_status_stream( m_stat_file )
{
  if(!m_status_stream)
    throw Dark_Arts_Exception(SUPPORT_OBJECT, "Could not create iteration Status_Generator filestream");

}

void Status_Generator::write_iteration_status(void)
{
  m_status_stream << "Time_step:  Stage: dt: number_of_inner_solves \n";
  return;
}  

