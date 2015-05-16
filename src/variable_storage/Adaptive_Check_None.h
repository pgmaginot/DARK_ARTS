#ifndef Adaptive_Check_None_h
#define Adaptive_Check_None_h


/** @file   Adaptive_Check_None.h
  *   @author pmaginot
  *   @brief To use when there is no adaptive time check for any variable type
 */
#include "V_Adaptive_Check.h"
class Adaptive_Check_None : public V_Adaptive_Check
{
public:
  Adaptive_Check_None( ): V_Adaptive_Check() {}
    
  virtual ~Adaptive_Check_None(){}

  bool adaptive_check(const int stage, const double dt) override { return false;}
  
private:
};

#endif