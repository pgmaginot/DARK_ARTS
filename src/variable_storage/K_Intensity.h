#ifndef K_Intensity_h
#define K_Intensity_h

#include "Cell_Data.h"
#include "Fem_Quadrature.h"
#include "Angular_Quadrature.h"
#include "Intensity_Data.h"
#include "Time_Data.h"

#include <vector>
#include <stdlib.h>

#include "Eigen/Dense"

class K_Intensity
{
public:
  /// Will set n_grp, n_el, n_dir, n_leg as static values
  K_Intensity(const int n_cells, const int n_stages, 
    const Fem_Quadrature& fem_quadrature, const Angular_Quadrature& angular_quadrature);
  ~K_Intensity(){}
  
  /// Public accessor functions
  void get_ki(const int cell, const int grp, const int dir, const int stage, Eigen::VectorXd& ki) const;

  /// Public functions to save values
  void set_ki(const int cell, const int grp, const int dir, const int stage, Eigen::VectorXd& ki);
  
  void advance_intensity(Intensity_Data& i_old, const double dt, const Time_Data* time_data);
private:
  
  /// total number of cells in the problem
  const int m_cells;
    
  /// number of DFEM unknowns in each cell
  const int m_el_per_cell;  
  
  /// number of DIRK stages
  const int m_n_stages; 
  
  const int m_n_dir;
  
  const int m_n_grp;
  
  /// total length of the intensity data
  const int m_k_length;  
  
  /**
    Constants for faster k_I indexing
  */
  const int m_el_stage;
  const int m_el_stage_dir_div_2;
  const int m_el_stage_dir_div_2_grp;
  const int m_offset;
  const int m_dir_div_2;
  
  Eigen::VectorXd m_vec_sum;
  Eigen::VectorXd m_vec_retrieve;
  
  std::vector<double> m_rk_b;
  
  std::vector<double> m_ki;
  
/* ***************************************************
  *
  *   Protected Functions
  *
  *************************************************** */
  
  bool ki_range_check(const int cell, const int grp, const int dir, const int stage) const;
  
  int ki_data_locator(const int cell, const int grp, const int dir, const int stage) const;
   
  bool ki_bounds_check(const int loc) const;

};

#endif