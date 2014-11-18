#ifndef Intensity_Moment_Data_h
#define Intensity_Moment_Data_h

#include "Fem_Quadrature.h"
#include "Angular_Quadrature.h"
#include "Cell_Data.h"
#include "Err_Phi.h"
#include "Eigen/Dense"
#include "Intensity_Data.h"

#include <vector>
#include <stdlib.h>

class Intensity_Moment_Data
{
public:
  /// used only by the phi_ard that is intiialized in Time_Marcher object
  Intensity_Moment_Data(const Cell_Data& cell_data, const Angular_Quadrature& ang_quad,
    const Fem_Quadrature& fem_quad, const Intensity_Data& i_old);

  /// used by subsequent phi .  Needed so that error comparisons can be made with some knowledge of what physically
  /// sized quantities are (aka what order of magnitude is really zero since we have unit issues)
  Intensity_Moment_Data(const Cell_Data& cell_data, const Angular_Quadrature& ang_quad,
    const Fem_Quadrature& fem_quad, const std::vector<double>& reference_phi_norm);  
  
  /// Copy constructor
  Intensity_Moment_Data(const Intensity_Moment_Data& intensity_moment);
  
  virtual ~Intensity_Moment_Data(){}
  
  /// assignment operator used during iteration processes
  Intensity_Moment_Data& operator= (const Intensity_Moment_Data& intensity_moment);
    
  /// return an Eigen compatible vector
  void get_cell_angle_integrated_intensity(const int cell, const int group, const int l_mom, 
    Eigen::VectorXd&  loc_phi_vec) const;
    
  void get_all_moments(std::vector<Eigen::VectorXd>& local_phi, 
    const int cell, const int grp) const;
    
  /// Public saver     
  void set_cell_angle_integrated_intensity(const int cell,
    const int group, const int l_mom, const Eigen::VectorXd& val);
    
  void clear_angle_integrated_intensity(void);
  
  void normalized_difference(Intensity_Moment_Data& phi_compare, Err_Phi& err_phi) const;
  
  void add_contribution(const int cell, const int grp, const int l_mom, Eigen::VectorXd& contrib);
  
  void get_phi_norm(std::vector<double>& norm_vec) const;  

private:  
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
   
  std::vector<double> m_norm_for_err;
  
  const double m_small_ratio;
  
  std::vector<double> m_phi;
  
  double get_angle_integrated_intensity(const int el, const int cell,
    const int group, const int l_mom) const;
  
  bool angle_integrated_range_check(const int el, const int cell, const int grp, const int l_mom) const;
    
  int angle_integrated_data_locator(const int el, const int cell, const int group, const int leg_mom) const;
    
  bool angle_integrated_bounds_check(const int loc) const;
  
  void calculate_reference_phi_norms(const Intensity_Data& i_old, 
    const Angular_Quadrature& ang_quad, 
    const Fem_Quadrature& fem_quad);
};

#endif
