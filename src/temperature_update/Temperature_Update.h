#ifndef Temperature_Update_h
#define Temperature_Update_h

#include "Grey_Temperature_Matrix_Creator.h"
#include "MF_Temperature_Matrix_Creator.h"
#include "Temperature_Data.h"
#include "Intensity_Moment_Data.h"
#include "K_Temperature.h"
#include "Fem_Quadrature.h"
#include "Cell_Data.h"
#include "Materials.h"
#include "Angular_Quadrature.h"
#include "Eigen/Dense"
#include "Err_Temperature.h"

#include <iomanip>
#include <iostream>

#include <memory>

/** @file   Temperature_Update.h
  *   @author pmaginot
  *   @brief Declare the Temperautre_Update class that will update a Temperature_Object given an Intensity_Object
 */

class Temperature_Update
{
public:
  Temperature_Update(const Fem_Quadrature& fem_quadrature, const Cell_Data& cell_data, Materials& material, 
    const Angular_Quadrature& angular_quadrature, const int n_stages, const Temperature_Data& t_old,
  const Intensity_Moment_Data& phi);
    
  virtual ~Temperature_Update(){}

  void update_temperature(Temperature_Data& t_star, const K_Temperature& k_t, const double damping, Err_Temperature& err_temperature);
    
  void set_time_data( const double dt, const double time_stage, const std::vector<double>& rk_a_of_stage_i, const int stage );

  void calculate_k_t(const Temperature_Data& t_star, K_Temperature& k_t);
protected:  
  /// number of DFEM points per cell
  const int m_np;
    
    /// number of cells
  const int m_n_cells;
  
  bool m_time_data_set;
  
  int m_stage;
  
  Materials& m_materials;
  
  Eigen::VectorXd m_t_star_vec;
  Eigen::VectorXd m_delta_vec;
  Eigen::VectorXd m_rhs_vec;
  Eigen::VectorXd m_k_vec;
  Eigen::MatrixXd m_coefficient_matrix;
  
  std::shared_ptr<V_Temperature_Matrix_Creator> m_matrix_creator;
    
  bool check_eigen_variables_finite(void) const;
  
  
};

#endif