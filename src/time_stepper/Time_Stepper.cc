/** @file   Time_Stepper.cc
  *   @author pmaginot
  *   @brief Implement the Time_Stepper class, SDIRK information, time step controlling
*/
#include "Time_Stepper.h"

Time_Stepper::Time_Stepper(Input_Reader&  input_reader)
{
  m_time_solver = input_reader.get_time_solver();
  
  if(m_time_solver == INVALID_TIME_SOLVER)
  {
    std::cerr << "Bad time solver in Time_Stepper constructor\n";
    exit(EXIT_FAILURE);
  }
  
  if( m_time_solver == IMPLICIT_EULER)
  {
    m_number_stages = 1;
  }
  
  m_a.resize(m_number_stages*(m_number_stages +1)/2 , 0. );
  m_b.resize(m_number_stages,0.);
  m_c.resize(m_number_stages,0.);
  
  fill_sdirk_vectors();
}


int Time_Stepper::get_number_of_stages(void) const
{
  return m_number_stages;
}

void Time_Stepper::fill_sdirk_vectors(void)
{
  if( m_time_solver == IMPLICIT_EULER)
  {
    m_a[0] = 1.;
    m_b[0] = 1.;
    m_c[0] = 1.;
  }
  
  return;
}
