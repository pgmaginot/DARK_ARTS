#ifndef Transport_Sweep_h
#define Transport_Sweep_h

#include "Angular_Quadrature.h"

/// this file is included from Sweep_Matrix_Creator_Grey and Sweep_Matrix_Creator_MF
#include "Sweep_Matrix_Creator_Grey.h"
#include "Sweep_Matrix_Creator_MF.h"

#include "Psi_In.h"

#include "Sweep_Fixed_Source_None.h"
#include "Sweep_Fixed_Source_Linearization.h"

#include "Solution_Saver_K_I.h"
#include "Solution_Saver_Flux_Moments.h"

/** @file   Transport_Sweep.h
  *   @author pmaginot
  *   @brief Declare the Transport sweep operator
 */

class Transport_Sweep
{
public:
  Transport_Sweep(const Fem_Quadrature& fem_quadrature, 
    const Cell_Data& cell_data, 
    Materials& material, 
    const Angular_Quadrature& angular_quadrature, 
    const int n_stages, 
    const Temperature_Data& t_old, 
    const Intensity_Data& i_old,
    const K_Temperature& kt, 
    K_Intensity& ki,
    const Temperature_Data& t_star);
    
  ~Transport_Sweep(){}

  /** set all of the first time data that the sweep will need
    * this includes: time stepping data (current stage number, rk_a, dt), 

    *
  */
  void set_time_data( const double dt, const double time_stage, const std::vector<double>& rk_a_of_stage_i, const int stage );
  
  /// sweep the mesh, calculating a phi_new
  void sweep_mesh(const Intensity_Moment_Data& phi_old, Intensity_Moment_Data& phi_new, const bool is_krylov, const bool is_k_i_sweep);
  
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
  
  const Angular_Quadrature& m_ang_quad;
    
  /// scratch matrix holder to be passed to matrix creator
  Eigen::MatrixXd m_matrix_scratch;
  
  /// scratch vector holder to be used as necessary
  Eigen::VectorXd m_vector_scratch;
  
  /// to hold all the flux moments in a given cell/group
  std::vector<Eigen::VectorXd> m_local_phi;
  
  /// rhs vector
  Eigen::VectorXd m_rhs_vec;
  
  /// lhs mamtrix
  Eigen::MatrixXd m_lhs_mat;
  
  /// new, local intensity solution
  Eigen::VectorXd m_local_soln;
  
  double m_time;
  
  Psi_In m_psi_in;    
  
  /// Creator of the linear boltzmann matrices (and source moment vector) that describe the transport sweep
  std::shared_ptr<V_Sweep_Matrix_Creator> m_sweep_matrix_creator;
  
  /// pointer used in call during @fn sweep_mesh()
  std::shared_ptr<V_Sweep_Fixed_Source> m_sweep_source;  
  /// use for parts of krylov solves, returns a zero source, smart ptr to a Sweep_Fixed_Source_None object that is created in Transport_Sweep constructor
  std::shared_ptr<V_Sweep_Fixed_Source> m_no_source;
  /// ptr that has access to the appropriate Grey/MF Sweep_Matrix_Creator shared_ptr that allows access to get_s_i
  std::shared_ptr<V_Sweep_Fixed_Source> m_fixed_source;
  
  
  /// pointer used in call during sweep  mesh to save (or not) the local intensity solution or the angle integrated intensity
  std::shared_ptr<V_Solution_Saver> m_sweep_saver;
  /// used to save \f$ k_I \f$ .  Does not change, overrwrite, or save angle integrated intensity.  Calculates k_I ONLY
  std::shared_ptr<V_Solution_Saver> m_k_i_saver;
  /// used during most normal sweeps to save angle integrated moments of the local solution
  std::shared_ptr<V_Solution_Saver> m_angle_integrated_saver;  
  
  void get_boundary_conditions(const bool is_krylov);
};

#endif