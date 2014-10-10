#ifndef Intensity_Moment_Data_h
#define Intensity_Moment_Data_h

#include "Fem_Quadrature.h"
#include "Angular_Quadrature.h"
#include "Cell_Data.h"
#include "Err_Phi.h"
#include "Eigen/Dense"

#include <vector>
#include <stdlib.h>

class Intensity_Moment_Data
{
public:
  /// Will set n_grp, n_el, n_dir, n_leg as static values
  Intensity_Moment_Data(const Cell_Data& cell_data, const Angular_Quadrature& ang_quad,
    const Fem_Quadrature& fem_quad);
  
  Intensity_Moment_Data(const int n_cells, const int n_grp, 
  const int n_leg_mom, const int n_el_cell);
  
  ~Intensity_Moment_Data(){}
  
  /// Copy constructor
  Intensity_Moment_Data(const Intensity_Moment_Data& intensity_moment);
  
  /// assignment operator
  Intensity_Moment_Data& operator= (const Intensity_Moment_Data& intensity_moment);
  
  /// Public accessor functions    
  double get_angle_integrated_intensity(const int el, const int cell,
    const int group, const int l_mom) const;
    
  void get_cell_angle_integrated_intensity(const int cell, const int group, const int l_mom, 
    std::vector<double>& loc_phi_vec) const;
  
  /// return an Eigen compatible vector
  void get_cell_angle_integrated_intensity(const int cell, const int group, const int l_mom, 
    Eigen::VectorXd&  loc_phi_vec) const;
    
  void get_all_moments(std::vector<Eigen::VectorXd>& local_phi, 
    const int cell, const int grp) const;
    
  /// Public saver functions
    
  void set_angle_integrated_intensity(const int el, const int cell,
    const int group, const int l_mom, const double val);
    
  void set_cell_angle_integrated_intensity(const int cell,
    const int group, const int l_mom, const std::vector<double>& val);
    
  void set_cell_angle_integrated_intensity(const int cell,
    const int group, const int l_mom, const Eigen::VectorXd& val);
    
  void clear_angle_integrated_intensity(void);
  
  void normalized_difference(Intensity_Moment_Data& phi_compare, Err_Phi& err_phi) const;
  
  void add_contribution(const int cell, const int grp, const int l_mom, Eigen::VectorXd& contrib);
  
private:
  std::vector<double> m_phi;
  
  /// total number of cells in the problem
  const int m_cells;
  /// total number of groups in the problem
  const int m_groups;
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
