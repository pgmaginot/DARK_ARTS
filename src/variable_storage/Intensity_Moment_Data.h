#ifndef Intensity_Moment_Data_h
#define Intensity_Moment_Data_h

#include "Fem_Quadrature.h"
#include "Angular_Quadrature.h"
#include "Cell_Data.h"
#include "Eigen/Dense"

#include <vector>
#include <stdlib.h>

class Intensity_Moment_Data
{
public:
  /// Will set n_grp, n_el, n_dir, n_leg as static values
  Intensity_Moment_Data(const Cell_Data& cell_data, const Angular_Quadrature& ang_quad,
    const Fem_Quadrature& fem_quad);
  ~Intensity_Moment_Data(){}
  
  /// Public accessor functions    
  double get_angle_integrated_intensity(const int el, const int cell,
    const int group, const int l_mom) const;
    
  void get_cell_angle_integrated_intensity(const int cell, const int group, const int l_mom, 
    std::vector<double>& loc_phi_vec) const;
  
  /// return an Eigen compatible vector
  void get_cell_angle_integrated_intensity(const int cell, const int group, const int l_mom, 
    Eigen::VectorXd&  loc_phi_vec) const;
    
  /// Public saver functions
    
  void set_angle_integrated_intensity(const int el, const int cell,
    const int group, const int l_mom, const double val);
    
  void set_cell_angle_integrated_intensity(const int cell,
    const int group, const int l_mom, const std::vector<double>& val);
    
  void set_cell_angle_integrated_intensity(const int cell,
    const int group, const int l_mom, const Eigen::VectorXd& val);
  
private:
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
    
  /// quantities necessary for faster indexing of angle integrated intensity
  const int m_el_times_l_mom;
  const int m_el_times_l_mom_times_group;
  
  /// total length of the angle integrated intensity data
  const int m_phi_length;
   
  
  bool angle_integrated_range_check(const int el, const int cell, const int grp, const int l_mom) const;
    
  int angle_integrated_data_locator(const int el, const int cell, const int group, const int leg_mom) const;
    
  bool angle_integrated_bounds_check(const int loc) const;

};

#endif
