#ifndef Diffusion_Matrix_Creator_Grey_h
#define Diffusion_Matrix_Creator_Grey_h

#include "V_Diffusion_Matrix_Creator.h"

/** @file   Diffusion_Matrix_Creator_Grey.h
  *   @author pmaginot
  *   @brief Provide a class that constructs grey TRT MIP matrices with planck lienarization
 */
 

class Diffusion_Matrix_Creator_Grey: public V_Diffusion_Matrix_Creator
{
public:
  /// Only able to initialize if given an Input_Reader object
  /// constructor defined in Fem_Quadrature.cc
  Diffusion_Matrix_Creator_Grey(const Fem_Quadrature& fem_quadrature, Materials& materials,
    const Angular_Quadrature& angular_quadrature, const Temperature_Data& t_star, const int n_cells, const Input_Reader& input_reader);
    
  virtual ~Diffusion_Matrix_Creator_Grey(){}
   
  /** MF needs to update M, r_cv, and spectrium \f$ \sum_{g=0}^G{\mathbf{R}_{\sigma_{a,g}} \mathbf{D}_g}   \f$
    Grey needs to update M, r_cv only
  */
  void set_time_data( const double dt, const double time_stage, const double sdirk_a_of_stage ) override;
    
  void calculate_pseudo_r_sig_a_and_r_sig_s(Eigen::MatrixXd& r_sig_a, Eigen::MatrixXd& r_sig_s) override;
    
protected:
  double m_rk_a_ii;
  double m_dt;
  const double m_c_speed;
  const double m_sum_sn_w;
  double m_tau;
  
  Eigen::MatrixXd m_mass;
  Eigen::MatrixXd m_dx_mass;
  Eigen::MatrixXd m_r_cv;  
  Eigen::MatrixXd m_d_matrix;
  Eigen::MatrixXd m_coefficient;
  const Eigen::MatrixXd m_identity_matrix;
  Eigen::MatrixXd m_temporary_mat;
  
  /// to evaluate \widetilde{D}, we need sig_a and sig_s  Store sig_s in d_vector (that V_Diffusion_Matrix_Creator posseses).
  /// store sigma_a in this vector
  std::vector<double> m_sig_a_for_d_coeff;
  
  void evaluate_diffusion_coefficents(double& d_r_cm1, double& d_l_c , double& d_r_c , double& d_l_cp1) override;
};

#endif