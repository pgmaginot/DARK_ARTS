#ifndef MMS_Temperature_h
#define MMS_Temperature_h

#include "Input_Reader.h"
#include "MMS_Time_Poly.h"
#include "MMS_Time_Cos.h"
#include "MMS_Space_Cos.h"
#include "MMS_Space_Poly.h"

class MMS_Temperature
{
public:
  /// Will set n_grp, n_el, n_dir, n_leg, m_i will be zero
  MMS_Temperature(const Input_Reader& input_reader);    
  virtual ~MMS_Temperature(){}  
  double get_mms_temperature(const double position , const double time);
  
private:    
  std::shared_ptr<V_MMS_Time> m_time_dep;
  std::shared_ptr<V_MMS_Space> m_temp_space_dep;  
};

#endif
