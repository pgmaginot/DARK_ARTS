#ifndef Time_Data_h
#define Time_Data_h

#include "Input_Reader.h"
#include "Inputs_Allowed.h"
// #include "Angular_Quadrature.h"
// #include "Fem_Quadrature.h"
// #include "Cell_Data.h"
// #include "Materials.h"

// #include "V_Temperature_Update.h"
// #include "Temperature_Update_Grey.h"
// #include "Temperature_Update_MF.h"

// // #include "V_Intensity_Update.h"
// #include "Intensity_Update_Grey.h"
// #include "Intensity_Update_MF.h"

// #include <vector>
// #include <memory>
// #include <stdlib.h>
// #include <iostream>

class Time_Data
{
public:
  /* ***************************************************
  *
  *   Public Functions
  *
  *************************************************** */
  Time_Data(const Input_Reader&  input_reader);
  ~Time_Data(){}
  
  int get_number_of_stages(void) const;
  
  double get_a(const int stage, const int index);
  
protected:  
  int m_number_stages;
  
  const TIME_SOLVER m_time_solver = INVALID_TIME_SOLVER;
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
  

  
  void fill_sdirk_vectors(void);
  
/* ***************************************************
  *
  *   Protected Functions
  *
  *************************************************** */

  void fill_sdirk_vectors(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c);
  
};

#endif