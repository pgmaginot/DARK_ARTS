#ifndef MMS_Intensity_h
#define MMS_Intensity_h

#include "Angular_Quadrature.h"
#include "Input_Reader.h"
#include "MMS_Time_Poly.h"
#include "MMS_Time_Cos.h"
#include "MMS_Angle_Poly.h"
#include "MMS_Angle_Isotropic.h"
#include "MMS_Space_Cos.h"
#include "MMS_Space_Poly.h"

class MMS_Intensity
{
public:
  /// Will set n_grp, n_el, n_dir, n_leg, m_i will be zero
  MMS_Intensity(const Input_Reader& input_reader, const Angular_Quadrature& angular_quadrature);    
  virtual ~MMS_Intensity(){}  
  double get_mms_intensity(const double position , const double time , const int dir);
  double get_mms_intensity(const double position , const double time , const double mu);

  double get_mms_phi(const double position , const double time );
  double get_mms_intensity_time_derivative(const double position , const double time , const int dir);
  double get_mms_intensity_space_derivative(const double position , const double time , const int dir);
private:  
  
  const Angular_Quadrature& m_ang_quad;
  std::shared_ptr<V_MMS_Time> m_time_dep;
  std::shared_ptr<V_MMS_Space> m_rad_space_dep;
  std::shared_ptr<V_MMS_Angle> m_angle_dep;
  double m_angle_integration;
  
};

#endif
