#include "DT_Calculator_Exponential.h"

DT_Calculator_Exponential::DT_Calculator_Exponential(const Input_Reader& input_reader)
  :
  V_DT_Calculator( input_reader ),
  m_ratio{input_reader.get_time_start_exponential_ratio() } ,
  m_step_max{ int(ceil( log(m_dt_max/m_dt_min) / log(m_ratio) )) }
{

}

double DT_Calculator_Exponential::calculate_dt(const int step, const double dt_old)
{
  double dt;
  if(step ==0)
  {
    dt = m_dt_min;
  }
  else 
  {
    dt = dt_old*m_ratio;
  }
  
  if( dt > m_dt_max)
    dt = m_dt_max;

  return dt;
}