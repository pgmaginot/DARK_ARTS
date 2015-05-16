#ifndef V_Adaptive_Check_h
#define V_Adaptive_Check_h

#include "Input_Reader.h"
#include "Dark_Arts_Exception.h"

/** @file   V_Adaptive_Check.h
  *   @author pmaginot
  *   @brief base class to allow for check whether time step must be cut
 */

class V_Adaptive_Check
{
public:
  V_Adaptive_Check(){}
    
  virtual ~V_Adaptive_Check(){}

  virtual bool adaptive_check(const int stage, const double dt, double& adapt_criteria) = 0;  
protected:
  
};

#endif