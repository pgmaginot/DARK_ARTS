
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
  Time_Marcher(const Input_Reader&  input_reader, Angular_Quadrature& angular_quadrature,
    const Fem_Quadrature& fem_quadrature, Cell_Data* const cell_data, Materials* const materials,  
    Temperature_Data& t_old, Intensity_Data& i_old,
    const Time_Data& time_data);
    
  ~Time_Marcher(){}
private:
  const int m_n_stages;

  K_Intensity m_k_i;
  K_Temperature m_k_t;
  Temperature_Data m_t_star;
  
  std::shared_ptr<V_Intensity_Update> m_intensity_update;
  std::shared_ptr<V_Temperature_Update> m_temperature_update;
};

#endif