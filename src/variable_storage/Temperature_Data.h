#ifndef Temperature_Data_h
#define Temperature_Data_h

#include "Fem_Quadrature.h"
#include "Input_Reader.h"
#include "Eigen/Dense"


class Temperature_Data
{
public:
  Temperature_Data(const int n_cells, 
    const Fem_Quadrature& fem_quad);
  
  /// Initiali condition constructor
  Temperature_Data(const int n_cells, 
    const Fem_Quadrature& fem_quad, 
    const Input_Reader& input_reader);
  
  ~Temperature_Data(){}
  
  /// Single element set/get functions
  double get_temperature(const int el, const int cell) const;

  void set_temperature(const int el, const int cell, const double val);

  /// single cell set/get functions
  void get_cell_temperature(const int cell, Eigen::VectorXd& vec) const;
  
  void set_cell_temperature(const int cell, const Eigen::VectorXd& vec);
  
  /// assignment operator
  Temperature_Data& operator= (const Temperature_Data& t_data);
  
  double calculate_average(void);
  
protected:
    /// total number of cells in the problem
  const int m_cells;
  
  /// number of DFEM unknowns in each cell
  const int m_el_per_cell;
  
  /// total length of the intensity data
  const int m_t_length;
    
  /// vector to hold temperature unknowns
  std::vector<double> m_t;
  
  std::vector<double> m_dfem_w;
  
/* ***************************************************
  *
  *   Protected Functions
  *
  *************************************************** */
  void set_cell_temperature(const int cell, const double val);
  
  bool temperature_range_check(const int el, const int cell) const;
  
  int temperature_data_locator(const int el, const int cell) const;
   
  bool temperature_bounds_check(const int loc) const;

};

#endif