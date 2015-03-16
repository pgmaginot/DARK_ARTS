#ifndef Diffusion_Operator_h
#define Diffusion_Operator_h

#include "Inputs_Allowed.h"
#include "Fem_Quadrature.h"
#include "Cell_Data.h"
#include "Materials.h"
#include "Angular_Quadrature.h"
#include "Temperature_Data.h"


#include "Eigen/Dense"

#include "Intensity_Moment_Data.h"

#include "Diffusion_Matrix_Creator_Grey.h"

#include "Grey_Diffusion_Ordering.h"
#include "MG_WG_Diffusion_Ordering.h"
#include "MG_LMFGA_Group_Collapse_Diffusion_Ordering.h"
#include "MG_LMFGA_No_Collapse_Cell_Outer_Diffusion_Ordering.h"
#include "MG_LMFGA_No_Collapse_Group_Outer_Diffusion_Ordering.h"

#include "Local_MIP_Assembler.h"
#include "MIP_Kappa_Calculator.h"

#include "MIP_Left_Boundary_Reflective.h"
#include "MIP_Left_Boundary_Incident_Flux.h"
#include "MIP_Right_Boundary_Incident_Flux.h"

#include <petscksp.h>

/** @file   Diffusion_Operator.h
  *   @author pmaginot
  *   @brief a class that creates a MIP diffusion matrix and inverts it, updating an intensity_moment_data object 
 */
class Diffusion_Operator
{
public:
  Diffusion_Operator(const Input_Reader& input_reader, const Fem_Quadrature& fem_quadrature, 
    const Cell_Data& cell_data, Materials& materials, const Angular_Quadrature& angular_quadrature,
    const int n_groups, const Temperature_Data& t_eval, const bool is_wg_solve);
    
  virtual ~Diffusion_Operator();

  bool check_all_eigen_variables_for_finite(void);
  
  void build_matrix(const int mip_loop_number);
  
  void build_rhs(const int mip_loop_number, const Intensity_Moment_Data& phi_new, const Intensity_Moment_Data& phi_old);
  
  void solve_system();
    
protected:
  const double m_sn_w;
  /// number of unknowns per DFEM 
  const int m_np;
  /// number of spatial cells
  const int m_n_cell;
  /// number of loops (groups) per solve (primarily for MG within group scattering DSA, all other times in 1)
  const int m_n_mip_loops;
  /// n_cell most times or n_cell x n_groups if doing LMFGA without group collapse
  const int m_n_mip_blocks;
  /// n_el x m_n_mip_blocks
  const int m_n_mip_sys_size;
  
  Eigen::MatrixXd m_r_sig_s;
  Eigen::MatrixXd m_r_sig_a;
  Eigen::MatrixXd m_s_matrix;
  
  Eigen::MatrixXd m_cell_cm1;
  Eigen::MatrixXd m_cell_c;
  Eigen::MatrixXd m_cell_cp1;
  
  double m_d_r_cm1 ,  m_d_l_c , m_d_r_c , m_d_l_cp1;
  double m_kappa_cm12 , m_kappa_cp12;
  double m_dx_cm1 , m_dx_c , m_dx_cp1;
  
  const Cell_Data& m_cell_data;
    
  double m_sdirk_a_stage;
  int m_dt;
  double m_time_stage;
  
  bool m_matrix_initial_build;
  
  MIP_Kappa_Calculator m_kappa_calculator;
    
  std::shared_ptr<V_Diffusion_Matrix_Creator> m_diffusion_matrix_creator;  
  
  std::shared_ptr<V_Diffusion_Ordering> m_diffusion_ordering;
  
  Local_MIP_Assembler m_local_assembler;  
  /// PETSc variables 
  PetscErrorCode m_petsc_err;  
  Mat m_mip_global;
  KSP m_krylov_solver;
  Vec m_mip_solution;
  Vec m_mip_rhs;
  
  PetscReal *m_pointer_to_eigen_data;
  int *m_vec_index_array;
  
};

#endif