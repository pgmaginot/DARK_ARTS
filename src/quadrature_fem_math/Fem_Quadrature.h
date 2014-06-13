#ifndef Fem_Quadrature_h
#define Fem_Quadrature_h

#include "Inputs_Allowed.h"
#include "Input_Reader.h"

#include <vector>
#include <stdlib.h>
#include <iostream>
#include <sstream>

/** @file   Fem_Quadrature.h
  *   @author pmaginot
  *   @brief Implement the Fem_Quadrature class that owns all things quadrature related
 */

class Fem_Quadrature
{
public:
  /// Only able to initialize if given an Input_Reader object
  /// constructor defined in 
  Fem_Quadrature(Input_Reader&  input_reader);
  ~Fem_Quadrature(){}
  
  
  
protected:
  const int m_xs_extra_points = 10;
  
  std::vector<double> m_dfem_interpolation_points;
  std::vector<double> m_dfem_interpolation_weights;
  std::vector<double> m_integration_points;
  std::vector<double> m_integration_weights;
  
  std::vector<double> m_xs_eval_points;
  std::vector<double> m_xs_eval_weights;
};

#endif