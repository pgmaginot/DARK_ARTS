#ifndef DT_Calculator_Vector_h
#define DT_Calculator_Vector_h


/** @file   DT_Calculator_Vector.h
  *   @author pmaginot
  *   @brief provide a interfaces for calculating \f$ \Delta t \f$ during time marching.  Right now only incorporates the need to
  * start with small time steps.  Could be changed/expanded to allow for adaptivity....
 */
#include "V_DT_Calculator.h"
class DT_Calculator_Vector : public V_DT_Calculator
{
public:
  DT_Calculator_Vector( const Input_Reader& input_reader);
    
  virtual ~DT_Calculator_Vector(){}

  double calculate_dt(const int step) override;
  
private:
  const int m_n_vector_stages;
  std::vector<double> m_dt_full_divisors;
  std::vector<int> m_small_steps;  
};

#endif