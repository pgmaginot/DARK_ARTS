#ifndef V_XS_Treatment_h
#define V_XS_Treatment_h

#include <vector>
#include <stdlib.h>

class V_XS_Treatment
{
public:
  V_XS_Treatment();
  virtual ~V_XS_Treatment();
  
  virtual void calculate_xs_at_integration_points(const std::vector<double>& xs_evaluations, 
    std::vector<double>& xs_at_dfem_integration_points) = 0;
};

#endif