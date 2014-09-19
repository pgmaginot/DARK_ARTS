#ifndef Time_Stepper_h
#define Time_Stepper_h

#include "Inputs_Allowed.h"
#include "Input_Reader.h"
#include "Angular_Quadrature.h"
#include "Fem_Quadrature.h"
#include "Cell_Data.h"
#include "Materials.h"

#include "V_Temperature_Update.h"
#include "Temperature_Update_Grey.h"
#include "Temperature_Update_MF.h"

#include "V_Intensity_Update.h"
#include "Intensity_Update_Grey.h"
#include "Intensity_Update_MF.h"

#include <vector>
#include <memory>
#include <stdlib.h>
#include <iostream>

class Time_Stepper
{
public:
  /// Only able to initialize if given input object and quadrature object
  Time_Stepper(const Input_Reader&  input_reader, const Angular_Quadrature& angular_quadrature,
    const Fem_Quadrature& fem_quadrature, Cell_Data* const cell_data, Materials* const materials);
    
  ~Time_Stepper(){}
  
  /* ***************************************************
  *
  *   Public Functions
  *
  *************************************************** */
  
  int get_number_of_stages(void) const;
  
  double get_a(const int stage, const int index);
  
protected:

  /**
   SDIRK time integration basics
   \f{eqnarray}{
        y^{n+1} &=& y^{n} + \Delta t \sum_{s=1}^{N_{stages}}{ b_s k_s} \\
        k_s &=& f(t^n + c_s \Delta t , y_s) \\
        y_s &=& y^n + \Delta t \sum_{j=1}^{s-1}{a_{s,j} k_j } + \Delta t a_{s,s} k_s 
   \f}
  */
  
  /// vector of length n_stages
  std::vector<double> m_b;
  /// vector of length n_stages
  std::vector<double> m_c;
  /// lower triangular matrix (represented as a vector, of size(n_stages*(n_stages+1)/2)
  std::vector<double> m_a;   
  
  int m_number_stages = -1;
  
  const TIME_SOLVER m_time_solver = INVALID_TIME_SOLVER;
  
  void fill_sdirk_vectors(void);
  
/* ***************************************************
  *
  *   Protected Functions
  *
  *************************************************** */

  void fill_sdirk_vectors(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c);
  
  std::shared_ptr<V_Intensity_Update> m_intensity_update;
  std::shared_ptr<V_Temperature_Update> m_temperature_update;
};

#endif