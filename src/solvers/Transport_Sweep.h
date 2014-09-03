#ifndef Transport_Sweep_h
#define Transport_Sweep_h

// #include "V_Sweep_Matrix_Creator.h"

#include "Temperature_Data.h"
#include "Intensity_Data.h"

#include "K_Temperature.h"
#include "K_Intensity.h"

#include "Fem_Quadrature.h"

#include "Cell_Data.h"

#include "Materials.h"

#include "Angular_Quadrature.h"
#include "Time_Stepper.h"

#include "Eigen/Dense"

#include <memory>

/** @file   Transport_Sweep.h
  *   @author pmaginot
  *   @brief Declare the Transport sweep operator
 */

class Transport_Sweep
{
public:
  Transport_Sweep(const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* material, 
    const Angular_Quadrature& angular_quadrature, const Time_Stepper& time_stepper);
    
    
  ~Transport_Sweep(){}

  void sweep_mesh(const bool is_krylov, Intensity_Data& intensity_new, const Intensity_Data& intensity_old,
    Temperature_Data& t_new, const Temperature_Data& t_star, const Temperature_Data& t_n,
    const K_Temperature& k_t, const K_Intensity& k_i, const int stage, 
    const std::vector<double>& rk_a, const double time, const double dt);
  
private:  
  /// number of cells
  const int m_n_cells;
  
  /// number of groups
  const int m_n_groups;
  
  /// number of directions
  const int m_n_dir;
  
  /// number of Legendre moments
  const int m_n_l_mom;
  
  /// number of DFEM points per cell
  const int m_np;    
    
  /// scratch matrix holder to be passed to matrix creator
  Eigen::MatrixXd m_matrix_scratch;
  
  /// rhs vector
  Eigen::VectorXd m_rhs_vec;
  
  /// lhs matrix
  Eigen::MatrixXd m_lhs_mat;
  
  /// new, local intensity solution
  Eigen::VectorXd m_local_soln;
  
  /// Creator of the linear boltzmann matrices (and source moment vector) that describe the transport sweep
  // std::shared_ptr<V_Sweep_Matrix_Creator> m_sweep_matrix_creator;
};

#endif