#include "DT_Calculator_Vector.h"

DT_Calculator_Vector::DT_Calculator_Vector(const Input_Reader& input_reader)
  :
  V_DT_Calculator( input_reader ),
  m_n_vector_stages( input_reader.get_number_of_vector_stages() )
{
  /// get the values we need from Input_Reader
  input_reader.get_time_start_vectors(m_dt_full_divisors, m_small_steps);
  
  /// store a cumulative sum of steps in each vector stage
  for(int i=1; i< m_n_vector_stages ; i++)
    m_small_steps[i] += m_small_steps[i-1];
}

double DT_Calculator_Vector::calculate_dt(const int step, const double dt_old, const double adapt_criteria)
{
  double dt = m_dt_max;  
  
  if( step < m_small_steps[m_n_vector_stages -1] )
  {  
    for(int i = 0 ; i < m_n_vector_stages ; i++)
    {
      if( step < m_small_steps[i] )
      {
        dt /= m_dt_full_divisors[i];
        break;
      }
    }
  }  
  
  check_dt(dt, step);
  
  return dt;
}