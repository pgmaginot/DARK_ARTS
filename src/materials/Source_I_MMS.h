#ifndef Source_I_MMS_h
#define Source_I_MMS_h

#include <vector>
#include <stdlib.h>
#include "VSource_I.h"
#include "Input_Reader.h"
#include "Angular_Quadrature.h"
#include "Planck.h"
#include "VAbsorption_Opacity.h"
#include "VScattering_Opacity.h"
#include "MMS_Time_Poly.h"
#include "MMS_Time_Cos.h"
#include "MMS_Angle_Poly.h"
#include "MMS_Angle_Isotropic.h"
#include "MMS_Space_Cos.h"
#include "MMS_Space_Poly.h"

class Source_I_MMS : public VSource_I
{
public:
  Source_I_MMS(const Input_Reader& input_reader, const Angular_Quadrature& angular_quadrature,
  std::vector<std::shared_ptr<VAbsorption_Opacity>>& abs_op , 
  std::vector<std::shared_ptr<VScattering_Opacity>>& scat_op, 
  const int mat_num , const Planck& planck);
  
  virtual ~Source_I_MMS(){}

  double get_intensity_source(const double position, 
    const int group, const int dir, const double time) override;
private:
  const Planck& m_planck;
  const Angular_Quadrature& m_ang_quad;
  const double m_c;
  double m_angle_integral=0.;
  std::shared_ptr<VScattering_Opacity> m_scat_op;
  std::shared_ptr<VAbsorption_Opacity> m_abs_op;
  
  std::shared_ptr<V_MMS_Time> m_time_dep;
  std::shared_ptr<V_MMS_Space> m_rad_space_dep;
  std::shared_ptr<V_MMS_Space> m_temp_space_dep;
  std::shared_ptr<V_MMS_Angle> m_angle_dep;
  
};

#endif