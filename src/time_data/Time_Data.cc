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
  m_t_start{input_reader.get_t_start() },
  m_current_dump_goal(0),
  m_n_extra_dumps(input_reader.get_n_data_dumps() ) , 
  m_times_to_dump(input_reader.get_dump_times_vector() )
{  
  if(m_time_solver == INVALID_TIME_SOLVER)
    throw Dark_Arts_Exception( SUPPORT_OBJECT , "Bad time solver in Time_Stepper constructor");
  
  if( m_time_solver == IMPLICIT_EULER)
  {
    m_number_stages = 1;
  }
  else if( m_time_solver == ALEXANDER_2_2)
  {
    m_number_stages = 2;
  }
  else if( m_time_solver == ALEXANDER_2_2_PLUS)
  {
    m_number_stages = 2;
  }
  else if( m_time_solver == ALEXANDER_3_3)
  {
    m_number_stages = 3;
  }
  
  m_a.resize(m_number_stages*(m_number_stages +1)/2 , 0. );
  m_b.resize(m_number_stages,0.);
  m_c.resize(m_number_stages,0.);
  
  fill_sdirk_vectors();
  
  /// set-up time starter
  STARTING_METHOD starting_method = input_reader.get_starting_time_method();
  if(starting_method == RAMP)
  {
    m_calculate_dt = std::make_shared<DT_Calculator_Ramp>( input_reader) ;
  }
  else if(starting_method == EXPONENTIAL)
  {
    m_calculate_dt = std::make_shared<DT_Calculator_Exponential>( input_reader );
  }
  else if(starting_method == VECTOR)
  {
    m_calculate_dt = std::make_shared<DT_Calculator_Vector>( input_reader ) ;
  }
  else if(starting_method == ADAPTIVE)
  {    
    ADAPTIVE_TIME_STEP_CONTROL meth = input_reader.get_adaptive_time_method();
    if(meth == CHANGE_IN_T)
    {
      m_calculate_dt = std::make_shared<DT_Calculator_Temperature_Change>( input_reader ) ;
    }
    else if(meth == CHANGE_IN_T_VOLUMETRIC)
    {
      m_calculate_dt = std::make_shared<DT_Calculator_Temperature_Change>( input_reader ) ;
     
    }
  }
  
  if(!m_calculate_dt)
    throw Dark_Arts_Exception( SUPPORT_OBJECT, "DT_Calculator not initilaized in Time_Data");
  
  /// put t_end at the end of dump times vector
  m_times_to_dump.push_back(m_t_end);
  
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
  if( m_time_solver == ALEXANDER_2_2)
  {
    double alpha = 1. - (sqrt(2.)/2.);
    
    m_a[0] = alpha;
    m_a[1] = 1.- alpha;   m_a[2] = alpha;
    
    m_b[0] = 1. - alpha;
    m_b[1] = alpha;
    
    m_c[0] = alpha;
    m_c[1] = 1.;
  }
  if( m_time_solver == ALEXANDER_2_2_PLUS)
  {
  
    double alpha = 1. + (sqrt(2.)/2.);
    
    m_a[0] = alpha;
    m_a[1] = 1.- alpha;   m_a[2] = alpha;
    
    m_b[0] = 1. - alpha;
    m_b[1] = alpha;
    
    m_c[0] = alpha;
    m_c[1] = 1.;
  }
  if( m_time_solver == ALEXANDER_3_3)
  {
    /// taken form SIAM J Num. Anal. V 14, n. 6 (1977)
    const double alpha = 0.4358665215084590;
    const double tau2 = (1.+alpha)/2.;
    const double b1 = -1.*(6.*alpha*alpha - 16.*alpha + 1.)/4.;
    const double b2 = (6.*alpha*alpha - 20.*alpha + 5.)/4.;
    m_a[0] = alpha;
    m_a[1] = tau2 - alpha;  m_a[2] = alpha;
    m_a[3] = b1;            m_a[4] = b2;  m_a[5] = alpha;
    
    m_b[0] = b1;            m_b[1] = b2;  m_b[2] = alpha;
    
    m_c[0] = alpha;
    m_c[1] = tau2;
    m_c[2] = 1.;
  }
  
  return;
}

double Time_Data::get_a(const int stage, const int index) const
{
  double a = 0.;
  if( (stage >= m_number_stages) || (stage < 0) )
    throw Dark_Arts_Exception( SUPPORT_OBJECT , "Invalid stage request for DIRK a constant");
  
  if( (index < 0) || (index > stage) )
    throw Dark_Arts_Exception( SUPPORT_OBJECT , "Index > stage  or index < 0 in request for DIRK a constant");
  
  if(stage == 0)
  {
    a = m_a[0];
  }
  else
  {
    int offset = stage*(stage+1)/2;    
    a =  m_a[offset + index];
  }
    
  return a;
}

double Time_Data::get_b(const int stage) const
{
  return m_b[stage];
}

double Time_Data::get_c(const int stage) const
{ 
  return m_c[stage];
}

double Time_Data::get_dt(const int step, const double time_now, const double dt_old, const double adapt_criteria)
{

  // std::cout << "Time_data adaptive criterion: " << adapt_criteria << std::endl;
  double dt = m_calculate_dt->calculate_dt(step,dt_old,adapt_criteria);
  if( (time_now + dt + m_dt_min) > m_times_to_dump[m_current_dump_goal] )
  {
    /// taking full time step suggested by starting method will end time past desired t_end
    dt = m_times_to_dump[m_current_dump_goal] - time_now;
    m_current_dump_goal++;
  }
  if( dt < 0. )
    throw Dark_Arts_Exception(TIME_MARCHER , "calculating negative dt");
    
  if(dt > m_dt_max)
    dt = m_dt_max;
    
  return dt;
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
