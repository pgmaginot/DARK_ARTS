#ifndef MMS_Time_Poly_h
#define MMS_Time_Poly_h

#include <vector>
#include <stdlib.h>
#include "V_MMS_Time.h"

class MMS_Time_Poly : public V_MMS_Time
{

public:
  MMS_Time_Poly(const std::vector<double>& coeff);
  virtual ~MMS_Time_Poly(){}

  double get_time_component(const double time) override;
  double get_time_derivative(const double time) override;
private:  
  const std::vector<double> m_poly_coeff;
  const int m_max_poly_degree_p1;
  double m_val;
  double m_pow;
};

#endif