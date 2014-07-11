#ifndef Intensity_Data_h
#define Intensity_Data_h

#include "Fem_Quadrature.h"
#include "Angular_Quadrature.h"
#include "Cell_Data.h"

#include <vector>
#include <stdlib.h>

class Intensity_Data
{
public:
  /// Will set n_grp, n_el, n_dir, n_leg as static values
  Intensity_Data(const Cell_Data& cell_data, const Angular_Quadrature& ang_quad,
    const Fem_Quadrature& fem_quad);
  ~Intensity_Data(){}
  
  /// Public accessor functions
  double get_intensity(const int el, const int cell,
    const int group, const int dir) const;
    
  double get_angle_integrated_intensity(const int el, const int cell,
    const int group, const int l_mom) const;
  
  /// Public functions to save values
  void set_intensity(const int el, const int cell,
    const int group, const int dir, const double val);
    
  void set_angle_integrated_intensity(const int el, const int cell,
    const int group, const int l_mom, const double val);
  
protected:
  std::vector<double> m_i;
  std::vector<double> m_phi;
  
  /// total number of cells in the problem
  int m_cells;
  /// total number of groups in the problem
  int m_groups;
  /// total number of directions in the problem
  int m_dir;
  
  /// total length of the intensity data
  int m_leg;
  
  /// number of DFEM unknowns in each cell
  int m_el_per_cell;
  
  
  /// total length of the intensity data
  int m_i_length;
  /// total length of the angle integrated intensity data
  int m_phi_length;
  
  /* ***************************************************
  *
  *   Protected Functions
  *
  *************************************************** */
  
  bool intensity_range_check(const int el, const int cell, const int grp, const int dir) const;
  
  bool angle_integrated_range_check(const int el, const int cell, const int grp, const int l_mom) const;
  
  int intensity_data_locator(const int el, const int cell, const int group, const int dir) const;
  
  int angle_integrated_data_locator(const int el, const int cell, const int group, const int leg_mom) const;
  
  bool intensity_bounds_check(const int loc) const;
  
  bool angle_integrated_bounds_check(const int loc) const;

};

#endif
