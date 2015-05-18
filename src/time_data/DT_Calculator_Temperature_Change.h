#ifndef DT_Calculator_Temperature_Change_h
#define DT_Calculator_Temperature_Change_h


/** @file   DT_Calculator_Vector.h
  *   @author pmaginot
  *   @brief provide a interfaces for calculating \f$ \Delta t \f$ during time marching.  Right now only incorporates the need to
  * start with small time steps.  Could be changed/expanded to allow for adaptivity....
 */
#include "V_DT_Calculator.h"
#include <iomanip>

class DT_Calculator_Temperature_Change : public V_DT_Calculator
{
public:
  DT_Calculator_Temperature_Change( const Input_Reader& input_reader);
    
  virtual ~DT_Calculator_Temperature_Change(){}

  double calculate_dt(const int step, const double dt_old, const double adapt_criteria) override;
  
private:
  const double m_goal_delta_temperature;
};

#endif