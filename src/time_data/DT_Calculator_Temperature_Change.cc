#include "DT_Calculator_Temperature_Change.h"

DT_Calculator_Temperature_Change::DT_Calculator_Temperature_Change(const Input_Reader& input_reader)
  :
  V_DT_Calculator( input_reader ),
  m_goal_delta_temperature( input_reader.get_t_change_adaptive_goal() )
{

}

double DT_Calculator_Temperature_Change::calculate_dt(const int step, const double dt_old, const double delta_temperature_last)
{
  double dt = 0.;
  if(step==0)
  {
    dt = m_dt_min*1000.;
  }
  else 
  {
    // std::cout << "delta t in dt calculator: " << std::scientific << std::setprecision(5) <<delta_temperature_last << std::endl;
    dt = std::min( m_goal_delta_temperature/delta_temperature_last*dt_old , dt_old*2.) ;
  }
  check_dt(dt, step);
  // std::cout << "dt: " << dt << std::endl;
  return dt;
}