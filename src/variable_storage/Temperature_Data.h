#ifndef Temperature_Data_h
#define Temperature_Data_h

#include "Fem_Quadrature.h"
#include "Cell_Data.h"

#include "Eigen/Dense"

#include <vector>
#include <stdlib.h>

class Temperature_Data
{
public:
  Temperature_Data(const Cell_Data& cell_data, const Fem_Quadrature& fem_quad, const int n_groups);
  ~Temperature_Data(){}
  
  /// Single element set/get functions
  double get_temperature(const int el, const int cell) const;

  void set_temperature(const int el, const int cell, const double val);

  /// single cell set/get functions
  void get_cell_temperature(const int cell, Eigen::VectorXd& vec) const;
  
  void set_cell_temperature(const int cell, const Eigen::VectorXd& vec);
  
  void get_cell_ard(const int cell, Eigen::VectorXd& vec) const;
  
  void set_cell_ard(const int cell, const Eigen::VectorXd& vec);
  
protected:
    /// total number of cells in the problem
  const int m_cells;
  
  /// number of DFEM unknowns in each cell
  const int m_el_per_cell;
  
  /// total length of the intensity data
  const int m_t_length;
  
  /// number of frequency groups  
  const int m_n_groups;
  
  /// total length of ard.  0 if a grey problem!
  const int m_ard_length;
  
  /// vector to hold temperature unknowns
  std::vector<double> m_t;
  
  /// vector to hold the absorption rate density (if a MF problem only)
  std::vector<double> m_ard;

  
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