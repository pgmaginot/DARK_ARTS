#include "DT_Calculator_Exponential.h"

DT_Calculator_Exponential::DT_Calculator_Exponential(const Input_Reader& input_reader)
  :
  V_DT_Calculator( input_reader ),
  m_ratio{input_reader.get_time_start_exponential_ratio() } ,
  m_step_max{ int(ceil( log(m_dt_max/m_dt_min) / log(m_ratio) )) }
{

}

double DT_Calculator_Exponential::calculate_dt(const int step)
{
  double dt = 0.;
  if( step < m_step_max)
  {
    dt = m_dt_min * pow(m_ratio,step);
  }
  else
  {
    dt = m_dt_max;
  }

  return dt;
}