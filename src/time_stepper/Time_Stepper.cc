/** @file   Time_Stepper.cc
  *   @author pmaginot
  *   @brief Implement the Time_Stepper class, SDIRK information, time step controlling
*/
#include "Time_Stepper.h"

Time_Stepper::Time_Stepper(const Input_Reader&  input_reader, const Angular_Quadrature& angular_quadrature,
    const Fem_Quadrature& fem_quadrature, Cell_Data* const cell_data, Materials* const materials)
  :
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
  if( angular_quadrature.get_number_of_groups() > 1){
    m_intensity_update = std::shared_ptr<V_Intensity_Update> (new Intensity_Update_MF(input_reader, fem_quadrature, cell_data, materials, angular_quadrature, m_number_stages) );
    m_temperature_update = std::shared_ptr<V_Temperature_Update> (new Temperature_Update_MF(fem_quadrature, cell_data, materials, angular_quadrature, m_number_stages) );
  }
  else{
    m_intensity_update = std::shared_ptr<V_Intensity_Update> (new Intensity_Update_Grey(input_reader,fem_quadrature, cell_data, materials, angular_quadrature, m_number_stages ) );
    m_temperature_update = std::shared_ptr<V_Temperature_Update> (new Temperature_Update_Grey( fem_quadrature, cell_data, materials, angular_quadrature, m_number_stages ) );
  }
}


int Time_Stepper::get_number_of_stages(void) const
{
  return m_number_stages;
}

void Time_Stepper::fill_sdirk_vectors(void)
{
  if( m_time_solver == IMPLICIT_EULER)
  {
    m_a[0] = 1.;
    m_b[0] = 1.;
    m_c[0] = 1.;
  }
  
  return;
}

double Time_Stepper::get_a(const int stage, const int index)
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