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
    
  virtual ~Diffusion_Operator(){}

  bool check_all_eigen_variables_for_finite(void);
    
protected:
  const double m_sn_w;
  const int m_np;
  
  Eigen::MatrixXd m_r_sig_s;
  Eigen::MatrixXd m_r_sig_a;
  
  Eigen::MatrixXd m_no_mu_pos_l_matrix;
  Eigen::MatrixXd m_no_mu_neg_l_matrix;
  
  double m_sdirk_a_stage;
  int m_dt;
  double m_time_stage;
  
  std::shared_ptr<V_Diffusion_Matrix_Creator> m_diffusion_matrix_creator;  
};

#endif