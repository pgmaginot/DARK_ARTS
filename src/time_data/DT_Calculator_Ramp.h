#ifndef DT_Calculator_Ramp_h
#define DT_Calculator_Ramp_h


/** @file   DT_Calculator_Ramp.h
  *   @author pmaginot
  *   @brief provide a interfaces for calculating \f$ \Delta t \f$ during time marching.  Right now only incorporates the need to
  * start with small time steps.  Could be changed/expanded to allow for adaptivity....
 */
#include "V_DT_Calculator.h"
class DT_Calculator_Ramp : public V_DT_Calculator
{
public:
  DT_Calculator_Ramp(const Input_Reader& input_reader);
    
  virtual ~DT_Calculator_Ramp(){}

  double calculate_dt(const int step, const double dt_old) override;
  
protected:
  const int m_n_ramp_steps;
  const double m_slope;
};

#endif