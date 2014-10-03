#ifndef K_Temperature_h
#define K_Temperature_h

#include "Cell_Data.h"
#include "Fem_Quadrature.h"


#include <vector>
#include <stdlib.h>

#include "Eigen/Dense"

class K_Temperature
{
public:
  /// Will set n_grp, n_el, n_dir, n_leg as static values
  K_Temperature(const int n_cells, const int n_stages, 
    const Fem_Quadrature& fem_quadrature);
  ~K_Temperature(){}
  
  /// Public accessor functions
  void get_kt(const int cell, const int stage, Eigen::VectorXd& kt) const;

  /// Public functions to save values
  void set_kt(const int cell, const int stage, Eigen::VectorXd& kt);
  
protected:
  
  /// total number of cells in the problem
  const int m_cells;
    
  /// number of DFEM unknowns in each cell
  const int m_el_per_cell;  
  
  /// number of DIRK stages
  const int m_n_stages; 
  
  /// total length of the intensity data
  const int m_k_length;  
  
  std::vector<double> m_kt;
  
/* ***************************************************
  *
  *   Protected Functions
  *
  *************************************************** */
  
  bool kt_range_check(const int cell, const int stage) const;
  
  int kt_data_locator(const int cell, const int stage) const;
   
  bool kt_bounds_check(const int loc) const;

};

#endif