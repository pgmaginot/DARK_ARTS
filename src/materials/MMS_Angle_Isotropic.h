#ifndef MMS_Angle_Isotropic_h
#define MMS_Angle_Isotropic_h

#include <vector>
#include <stdlib.h>
#include "V_MMS_Angle.h"

class MMS_Angle_Isotropic : public V_MMS_Angle
{

public:
  MMS_Angle_Isotropic(const double sum_w);
  virtual ~MMS_Angle_Isotropic(){}

  double get_angle_component(const int dir) override;
  double get_angle_component(const double mu) override;
private:  
  const double m_sum_w_inv;
};

#endif