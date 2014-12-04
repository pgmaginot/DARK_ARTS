#ifndef MMS_Angle_Poly_h
#define MMS_Angle_Poly_h

#include <vector>
#include <stdlib.h>
#include "V_MMS_Angle.h"
#include "Angular_Quadrature.h"

class MMS_Angle_Poly : public V_MMS_Angle
{

public:
  MMS_Angle_Poly(const std::vector<double>& poly_coeff, const Angular_Quadrature& angular_quadrature);
  virtual ~MMS_Angle_Poly(){}

  double get_angle_component(const int dir) override;
private:  
  const std::vector<double> m_poly_coeff;
  const int m_max_poly_degree_p1;  
  std::vector<double> m_angles;
  double m_mu;
  double m_val;
  double m_pow;
};

#endif