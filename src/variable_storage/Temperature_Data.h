#ifndef Temperature_Data_h
#define Temperature_Data_h

#include "Fem_Quadrature.h"
#include "Cell_Data.h"

#include <vector>
#include <stdlib.h>

class Temperature_Data
{
public:
  /// Will set n_grp, n_el, n_dir, n_leg as static values
  Temperature_Data(const Cell_Data& cell_data, const Fem_Quadrature& fem_quad);
  ~Temperature_Data(){}
  
  /// Public accessor functions
  double get_temperature(const int el, const int cell) const;

  /// Public functions to save values
  void set_temperature(const int el, const int cell, const double val);
  
protected:
  std::vector<double> m_t;
  
  /// total number of cells in the problem
  int m_cells;
  
  /// number of DFEM unknowns in each cell
  int m_el_per_cell;
  
  
  /// total length of the intensity data
  int m_t_length;
  
/* ***************************************************
  *
  *   Protected Functions
  *
  *************************************************** */
  
  bool temperature_range_check(const int el, const int cell) const;
  
  int temperature_data_locator(const int el, const int cell) const;
   
  bool temperature_bounds_check(const int loc) const;

};

#endif