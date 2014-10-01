#ifndef Transport_Sweep_h
#define Transport_Sweep_h

#include "Fem_Quadrature.h"
#include "Cell_Data.h"
#include "Materials.h"
#include "Angular_Quadrature.h"

#include "Temperature_Data.h"
#include "Intensity_Data.h"
#include "Intensity_Moment_Data.h"
#include "K_Temperature.h"
#include "K_Intensity.h"

#include "V_Sweep_Matrix_Creator.h"
#include "Sweep_Matrix_Creator_Grey.h"
#include "Sweep_Matrix_Creator_MF.h"

#include "Psi_In.h"

#include <memory>

/** @file   Transport_Sweep.h
  *   @author pmaginot
  *   @brief Declare the Transport sweep operator
 */

class Transport_Sweep
{
public:
  Transport_Sweep(const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* material, 
    Angular_Quadrature& angular_quadrature, const int n_stages);
    
  ~Transport_Sweep(){}

  void sweep_mesh(const Intensity_Moment_Data& phi_old, Intensity_Moment_Data& phi_new, const bool is_krylov);
  
  void set_ard_phi_ptr(Intensity_Moment_Data* ard_phi_ptr);
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
  
  const double m_sum_w;
  
  Angular_Quadrature* const m_ang_quad;
    
  /// scratch matrix holder to be passed to matrix creator
  Eigen::MatrixXd m_matrix_scratch;
  
  /// rhs vector
  Eigen::VectorXd m_rhs_vec;
  
  /// lhs matrix
  Eigen::MatrixXd m_lhs_mat;
  
  /// new, local intensity solution
  Eigen::VectorXd m_local_soln;
  
  double m_time;
  
  Psi_In m_psi_in;
    
  
  /// Creator of the linear boltzmann matrices (and source moment vector) that describe the transport sweep
  std::shared_ptr<V_Sweep_Matrix_Creator> m_sweep_matrix_creator;
  
  void get_boundary_conditions(Psi_In& psi_in,const bool is_krylov);
  
  
};

#endif