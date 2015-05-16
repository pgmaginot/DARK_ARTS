#include "DT_Calculator_Temperature_Change.h"

DT_Calculator_Temperature_Change::DT_Calculator_Temperature_Change(const Input_Reader& input_reader)
  :
  V_DT_Calculator( input_reader ),
  m_goal_delta_temperature( input_reader.get_t_change_adaptive_goal() )
{

}

double DT_Calculator_Temperature_Change::calculate_dt(const int step, const double dt_old, const double delta_temperature_last)
{
  double dt = m_goal_delta_temperature/delta_temperature_last*dt_old;
  
  check_dt(dt, step);
  
  return dt;
}