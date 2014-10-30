
#ifndef Time_Marcher_h
#define Time_Marcher_h

#include "Time_Data.h"
#include "Intensity_Update_Grey.h"
#include "Intensity_Update_MF.h"

#include "Temperature_Update_Grey.h"
#include "Temperature_Update_MF.h"

class Time_Marcher
{

public:
  Time_Marcher(const Input_Reader&  input_reader, const Angular_Quadrature& angular_quadrature,
    const Fem_Quadrature& fem_quadrature, const Cell_Data& cell_data, Materials& materials,  
    Temperature_Data& t_old, Intensity_Data& i_old,
    const Time_Data& time_data);
    
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
  Err_Temperature m_err_temperature;
  
  std::shared_ptr<V_Intensity_Update> m_intensity_update;
  std::shared_ptr<V_Temperature_Update> m_temperature_update;
};

#endif