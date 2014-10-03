#ifndef V_DT_Calculator_h
#define V_DT_Calculator_h


/** @file   V_DT_Calculator.h
  *   @author pmaginot
  *   @brief provide a interfaces for calculating \f$ \Delta t \f$ during time marching.  Right now only incorporates the need to
  * start with small time steps.  Could be changed/expanded to allow for adaptivity....
 */

class V_DT_Calculator
{
public:
  V_DT_Calculator(){}
    
  virtual ~V_DT_Calculator(){}

  virtual double calculate_dt(const int step) = 0;
};

#endif