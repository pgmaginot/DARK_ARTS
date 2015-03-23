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
  virtual ~Cell_Data(){}
  
  /// Public accesor functions
  double get_cell_width(int cell_num) const;
  double get_cell_left_edge(int cell_num) const;
  int get_cell_material_number(int cell_num) const;
  
  int get_total_number_of_cells(void) const;  
protected:
  std::vector<double> m_x_l;
  std::vector<double> m_dx;
  std::vector<int> m_material_num;    
  
  int m_total_cells=0;
  
/* ***************************************************
  *
  *   Protected Functions
  *
  *************************************************** */
  
  void determine_cell_properties(const int n_reg, const std::vector<int>& cell_reg,
    const Input_Reader& input_reader);
    
  // void read_mesh_and_refine(const Input_Reader& input_reader);
};

#endif