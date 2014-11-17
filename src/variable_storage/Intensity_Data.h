#ifndef Intensity_Data_h
#define Intensity_Data_h

#include "Fem_Quadrature.h"
#include "Angular_Quadrature.h"
#include "Cell_Data.h"
#include "Eigen/Dense"
#include "Materials.h"
#include "Input_Reader.h"

class Intensity_Data
{
public:
  /// Will set n_grp, n_el, n_dir, n_leg, m_i will be zero
  Intensity_Data(const Cell_Data& cell_data, 
    const Angular_Quadrature& ang_quad,
    const Fem_Quadrature& fem_quad);
  
  Intensity_Data(const Cell_Data& cell_data, 
    const Angular_Quadrature& ang_quad,
    const Fem_Quadrature& fem_quad, 
    Materials& materials,
    const Input_Reader& input_reader);
    
  virtual ~Intensity_Data(){}
  
  double get_intensity(const int el, const int cell, const int group, const int dir) const;
  
 
  
  /// Public accessor functions
  void get_cell_intensity(const int cell, const int group, const int dir, 
    Eigen::VectorXd& loc_i_vec) const;
  
  /// Public functions to save values    
  void set_cell_intensity(const int cell,
    const int group, const int dir, const Eigen::VectorXd& val);
     
private:  
  /// total number of cells in the problem
  const int m_cells;
  /// total number of groups in the problem
  const int m_groups;
  /// total number of directions in the problem
  const int m_dir;
    /// number of DFEM unknowns in each cell
  const int m_el_per_cell;
  
  /// quantities necessary for faster indexing of intensity
  const int m_dir_div_2;
  
  /// offset from negative mu to positive mu data ordering
  const int m_offset;
  const int m_el_times_dir_div_2;
  const int m_el_times_dir_div_2_times_grp;
  
  /// total length of the intensity data
  const int m_i_length;
  
  /// actual storage structure to store intensities
  std::vector<double> m_i;
  
  /* ***************************************************
  *
  *   Protected Functions
  *
  *************************************************** */
  void set_cell_intensity(const int cell,
    const int group, const int dir, const double val);  
    
  bool intensity_range_check(const int el, const int cell, const int grp, const int dir) const;
  
  int intensity_data_locator(const int el, const int cell, const int group, const int dir) const;
  
  bool intensity_bounds_check(const int loc) const;

};

#endif
