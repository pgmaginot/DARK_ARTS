#ifndef Intensity_Data_h
#define Intensity_Data_h

#include "Fem_Quadrature.h"
#include "Angular_Quadrature.h"
#include "Cell_Data.h"
#include "Eigen/Dense"

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
    
  void get_cell_intensity(const int cell, const int group, const int dir, 
    std::vector<double>& loc_i_vec) const;
    
  /// return an Eigen compatible vector
  void get_cell_intensity(const int cell, const int group, const int dir, 
    Eigen::VectorXd& loc_i_vec) const;
    
  double get_angle_integrated_intensity(const int el, const int cell,
    const int group, const int l_mom) const;
    
  void get_cell_angle_integrated_intensity(const int cell, const int group, const int l_mom, 
    std::vector<double>& loc_phi_vec) const;
  
  /// return an Eigen compatible vector
  void get_cell_angle_integrated_intensity(const int cell, const int group, const int l_mom, 
    Eigen::VectorXd&  loc_phi_vec) const;
  
  /// Public functions to save values
  void set_intensity(const int el, const int cell,
    const int group, const int dir, const double val);
    
  void set_cell_intensity(const int cell,
    const int group, const int dir, const std::vector<double>& val);
    
  void set_cell_intensity(const int cell,
    const int group, const int dir, const Eigen::VectorXd& val);
    
  void set_angle_integrated_intensity(const int el, const int cell,
    const int group, const int l_mom, const double val);
    
  void set_cell_angle_integrated_intensity(const int cell,
    const int group, const int l_mom, const std::vector<double>& val);
    
  void set_cell_angle_integrated_intensity(const int cell,
    const int group, const int l_mom, const Eigen::VectorXd& val);
  
private:
  std::vector<double> m_i;
  std::vector<double> m_phi;
  
  /// total number of cells in the problem
  const int m_cells;
  /// total number of groups in the problem
  const int m_groups;
  /// total number of directions in the problem
  const int m_dir;
    /// number of legendre moments to store of the full intensity
  const int m_leg;
    /// number of DFEM unknowns in each cell
  const int m_el_per_cell;
  
  /// quantities necessary for faster indexing of intensity
  const int m_dir_div_2;
  
  /// offset from negative mu to positive mu data ordering
  const int m_offset;
  const int m_el_times_dir_div_2;
  const int m_el_times_dir_div_2_times_grp;
  
  /// quantities necessary for faster indexing of angle integrated intensity
  const int m_el_times_l_mom;
  const int m_el_times_l_mom_times_group;
  
  /// total length of the intensity data
  const int m_i_length;
  /// total length of the angle integrated intensity data
  const int m_phi_length;
  
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
