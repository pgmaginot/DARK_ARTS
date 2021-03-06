#ifndef V_DT_Calculator_h
#define V_DT_Calculator_h

#include "Input_Reader.h"
#include "Dark_Arts_Exception.h"


/** @file   V_DT_Calculator.h
  *   @author pmaginot
  *   @brief provide a interfaces for calculating \f$ \Delta t \f$ during time marching.  Right now only incorporates the need to
  * start with small time steps.  Could be changed/expanded to allow for adaptivity....
 */

class V_DT_Calculator
{
public:
  V_DT_Calculator(const Input_Reader& input_reader);
    
  virtual ~V_DT_Calculator(){}

  virtual double calculate_dt(const int step, const double dt_old, const double adapt_criteria) = 0;
  
protected:
  const double m_dt_min;
  const double m_dt_max;  
  
  void check_dt(double& dt , const int step);
};

#endif