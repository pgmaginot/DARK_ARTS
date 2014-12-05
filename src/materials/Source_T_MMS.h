#ifndef Source_T_MMS_h
#define Source_T_MMS_h

#include "VSource_T.h"
#include "Input_Reader.h"
#include "MMS_Space_Poly.h"
#include "MMS_Space_Cos.h"
#include "MMS_Time_Poly.h"
#include "MMS_Time_Cos.h"
#include "MMS_Angle_Poly.h"
#include "MMS_Angle_Isotropic.h"
#include "VCv.h"
#include "VAbsorption_Opacity.h"
#include "Angular_Quadrature.h"
#include "Planck.h"

class Source_T_MMS : public VSource_T
{

public:
  Source_T_MMS(const Input_Reader& input_reader, 
    const Angular_Quadrature& angular_quadrature,
    std::vector<std::shared_ptr<VAbsorption_Opacity> >& abs_op, 
    std::vector<std::shared_ptr<VCv> >& cv, 
    const int mat_num, const Planck& planck) ;
  virtual ~Source_T_MMS(){}

  double get_temperature_source(const double position, const double time) override;
private:
  const Planck& m_planck;
  const double m_sn_w;
  double m_angle_integration;
  /// intesnity = time_dep*ang_dep*space_dep
  /// phi = \sum_{d} {w_d intensity_d} = \sum_d {time_dep*space_dep*(ang_dep_d) }
  /// phi = space_dep*time_dep* ( \sum_{d} w_d ang_dep_d )
  std::shared_ptr<V_MMS_Space> m_rad_space_dep;
  std::shared_ptr<V_MMS_Space> m_temp_space_dep;
  std::shared_ptr<V_MMS_Time> m_time_dep;  
  std::shared_ptr<VAbsorption_Opacity> m_abs_op;
  std::shared_ptr<VCv> m_cv;
};

#endif