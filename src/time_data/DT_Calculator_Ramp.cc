#include "DT_Calculator_Ramp.h"

DT_Calculator_Ramp::DT_Calculator_Ramp(const Input_Reader& input_reader)
:
  V_DT_Calculator( input_reader ),
  m_n_ramp_steps{ input_reader.get_number_of_ramp_steps() },
  m_slope{ (m_dt_max - m_dt_min)/( double(m_n_ramp_steps) ) }
{

}

double DT_Calculator_Ramp::calculate_dt(const int step)
{
  double dt = 0.;
  /// linearly ramp time step from dt_min (step==0) to dt_max step >= 
  if( step < m_n_ramp_steps )
  {
    dt = m_dt_min + m_slope* ( double(step) );
  }
  else
  {
    dt = m_dt_max;
  }
  
  return dt;
}
