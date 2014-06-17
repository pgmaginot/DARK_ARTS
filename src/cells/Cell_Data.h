#ifndef Cell_Data_h
#define Cell_Data_h

#include "Inputs_Allowed.h"
#include "Input_Reader.h"
#include "Fem_Quadrature.h"

#include <vector>
#include <stdlib.h>
#include <iostream>
#include <sstream>

class Fem_Quadrature
{
public:
  /// Only able to initialize if given an Input_Reader object
  Fem_Quadrature(const Input_Reader&  input_reader);
  ~Fem_Quadrature();
  
protected:
  std::vector<double> m_x_l;
  std::vector<double> m_x_r;
  std::vector<double> m_dx;
  std::vector<int> m_material_num;
    
};

#endif