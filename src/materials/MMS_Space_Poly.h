#ifndef MMS_Space_Poly_h
#define MMS_Space_Poly_h

#include <vector>
#include <stdlib.h>
#include "V_MMS_Space.h"

class MMS_Space_Poly : public V_MMS_Space
{

public:
  MMS_Space_Poly(const std::vector<double>& coeff);
  virtual ~MMS_Space_Poly(){}

  double get_position_component(const double position) override;
  double get_position_derivative(const double position) override;
private:  
  const std::vector<double> m_poly_coeff;
  const int m_max_poly_degree_p1;
  double m_val;
  double m_pow;
};

#endif