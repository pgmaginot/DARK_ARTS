#include "V_DT_Calculator.h"

V_DT_Calculator::V_DT_Calculator(const Input_Reader& input_reader)
:
  m_dt_min( input_reader.get_dt_min() ),
  m_dt_max( input_reader.get_dt_max() )
{

}

void V_DT_Calculator::check_dt(double& dt , const int step)
{
  if(dt < m_dt_min)
  {
    std::stringstream err;
    err << "Calculated dt= : " << dt << "for step: " << step << std::endl;
    err << "Minimum dt= " << m_dt_min;
    throw Dark_Arts_Exception(SUPPORT_OBJECT, err.str() ) ;
  }
  if(dt > m_dt_max)
  {
    dt = m_dt_max;
  }
  
  return;
}