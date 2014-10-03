/** @file   Time_Data.cc
  *   @author pmaginot
  *   @brief Implement the Time_Stepper class, SDIRK information, time step controlling
*/
#include "Time_Data.h"
Time_Data::Time_Data(const Input_Reader&  input_reader)
  :
  m_number_stages{-1},
  m_time_solver{input_reader.get_time_solver()},
  m_dt_min{input_reader.get_dt_min() },
  m_dt_max{input_reader.get_dt_max() },
  m_t_end{input_reader.get_t_end() },
  m_t_start{input_reader.get_t_start() }
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
  
  /// set-up time starter
  STARTING_METHOD starting_method = input_reader.get_starting_time_method();
  if(starting_method == RAMP)
  {
    m_calculate_dt = std::shared_ptr<V_DT_Calculator> (new DT_Calculator_Ramp() );
  }
  else if(starting_method == EXPONENTIAL)
  {
    m_calculate_dt = std::shared_ptr<V_DT_Calculator> (new DT_Calculator_Exponential() );
  }
  else if(starting_method == VECTOR)
  {
    m_calculate_dt = std::shared_ptr<V_DT_Calculator> (new DT_Calculator_Vector() );
  }
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

double Time_Data::get_a(const int stage, const int index) const
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

double Time_Data::get_c(const int stage) const
{ 
  return m_c[stage];
}

double Time_Data::get_dt(const int step)
{
  return m_calculate_dt->calculate_dt(step);
}

double Time_Data::get_t_start(void) const
{
  return m_t_start;
}

double Time_Data::get_t_end(void) const
{
  return m_t_end;
}

double Time_Data::get_dt_min(void) const
{
  return m_dt_min;
}

double Time_Data::get_dt_max(void) const
{
  return m_dt_max;
}