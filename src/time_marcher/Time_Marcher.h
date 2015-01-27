
#ifndef Time_Marcher_h
#define Time_Marcher_h

#include "Time_Data.h"
#include "Intensity_Update_Grey.h"
#include "Intensity_Update_MF.h"

#include "Temperature_Update_Grey.h"
#include "Temperature_Update_MF.h"
#include "Status_Generator.h"
#include "Output_Generator.h"
#include "Final_Space_Error_Calculator.h"
#include "Space_Time_Error_Calculator.h"

class Time_Marcher
{

public:
  Time_Marcher(const Input_Reader&  input_reader, const Angular_Quadrature& angular_quadrature,
    const Fem_Quadrature& fem_quadrature, const Cell_Data& cell_data, Materials& materials,  
    Temperature_Data& t_old, Intensity_Data& i_old,
    const Time_Data& time_data, std::string input_filename);
    
  virtual ~Time_Marcher(){}
  
  void solve(Intensity_Data& i_old, Temperature_Data& t_old, Time_Data& time_data);
private:
  const int m_n_stages;

  const Time_Data& m_time_data;
  
  const double m_thermal_tolerance;
  
  K_Intensity m_k_i;
  K_Temperature m_k_t;
  Temperature_Data m_t_star;
  Intensity_Moment_Data m_ard_phi;
  
  /// damping variable, and temperature update err object used to track thermal iteration progress
  double m_damping;
  int m_iters_before_damping;
  const double m_damping_decrease_factor;
  const int m_iteration_increase_factor;
  const int m_checkpoint_frequency;
  const int m_max_damps;
  
  Err_Temperature m_err_temperature;
  Status_Generator m_status_generator;
  Output_Generator m_output_generator;
  
  const bool m_calculate_space_time_error;
  const bool m_calculate_final_space_error;
  
  std::shared_ptr<Space_Time_Error_Calculator> m_space_time_error_calculator;
  std::shared_ptr<Final_Space_Error_Calculator> m_final_space_error_calculator;
  
  std::shared_ptr<V_Intensity_Update> m_intensity_update;
  std::shared_ptr<V_Temperature_Update> m_temperature_update;
};

#endif