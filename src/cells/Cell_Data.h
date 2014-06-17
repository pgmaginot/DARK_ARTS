#ifndef Cell_Data_h
#define Cell_Data_h

#include "Inputs_Allowed.h"
#include "Input_Reader.h"

#include <vector>
#include <stdlib.h>
#include <iostream>
#include <sstream>

class Cell_Data
{
public:
  /// Only able to initialize if given input object and quadrature object
  Cell_Data(Input_Reader&  input_reader);
  ~Cell_Data(){}
  
protected:
  std::vector<double> m_x_l;
  std::vector<double> m_x_r;
  std::vector<int> m_material_num;    
  
  int m_total_cells=0;
};

#endif