/** @file   Time_Data.cc
  *   @author pmaginot
  *   @brief Implement the Time_Stepper class, SDIRK information, time step controlling
*/
#include "Time_Data.h"
Time_Data::Time_Data(const Input_Reader&  input_reader)
  :
  m_number_stages{-1},
  m_time_solver{input_reader.get_time_solver()}
{

  
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
  
  /// determine if grey or MF, set-up solvers appropriately

}


int Time_Data::get_number_of_stages(void) const
{
  return m_number_stages;
}

void Time_Data::fill_sdirk_vectors(void)
{
  if( m_time_solver == IMPLICIT_EULER)
  {
    m_a[0] = 1.;
    m_b[0] = 1.;
    m_c[0] = 1.;
  }
  
  return;
}

double Time_Data::get_a(const int stage, const int index)
{
  if( (stage >= m_number_stages) || (stage < 0) )
  {
    std::cerr << "Invalid stage request for DIRK a constant\n";
    exit(EXIT_FAILURE);
  }
  
  if( (index < 0) || (index > stage) )
  {
    std::cerr << "Index > stage  or index < 0 in request for DIRK a constant\n";
    exit(EXIT_FAILURE);
  }
  
  if(stage == 0)
  {
    return m_a[0];
  }
  else
  {
    return m_a[stage*(stage+1)/2 + index];
  }
}