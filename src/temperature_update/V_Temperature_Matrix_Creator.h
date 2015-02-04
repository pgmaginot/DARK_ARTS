#ifndef V_Temperature_Matrix_Creator_h
#define V_Temperature_Matrix_Creator_h


#include "Temperature_Data.h"
#include "Intensity_Moment_Data.h"
#include "Matrix_Construction_Exact.h"
#include "Matrix_Construction_Self_Lumping.h"
#include "Matrix_Construction_Trad_Lumping.h"
#include "Eigen/Dense"
#include "K_Temperature.h"
#include "Cell_Data.h"
#include "Angular_Quadrature.h"
#include "Fem_Quadrature.h"

/** @file   V_Temperature_Matrix_Creator.h
  *   @author pmaginot
 */

class V_Temperature_Matrix_Creator
{
public:
  V_Temperature_Matrix_Creator(const Fem_Quadrature& fem_quadrature, const Cell_Data& cell_data, Materials& material, 
    const Angular_Quadrature& angular_quadrature, const int n_stages, const Temperature_Data& t_old, const Intensity_Moment_Data& phi);
    
  virtual ~V_Temperature_Matrix_Creator(){}
  
  virtual void calculate_update_quantities(const int cell, const Eigen::VectorXd& t_star, const K_Temperature& k_t,
    Eigen::MatrixXd& coefficient , Eigen::VectorXd& rhs) = 0;
    
  virtual void calculate_k_t(const int cell, const Eigen::VectorXd& t_star, Eigen::VectorXd& k_t) = 0;
    
  void set_time_data(const std::vector<double>& rk_a, const double dt, const double time_stage, const int stage_num);    
protected:
  const int m_np;
  std::vector<double> m_rk_a;
  double m_dt;
  double m_time;
  int m_stage;  
  Materials& m_materials;
  const double m_sn_w;
  const Temperature_Data& m_t_old;
  const Intensity_Moment_Data& m_phi;
  const Eigen::MatrixXd m_identity_matrix;
  Eigen::MatrixXd m_d_matrix;
  Eigen::MatrixXd m_r_sig_a;  
  Eigen::MatrixXd m_r_cv;
  
  Eigen::VectorXd m_t_old_vec;
  Eigen::VectorXd m_driving_source_vec;
  Eigen::VectorXd m_phi_vec;
  Eigen::VectorXd m_planck_vec;
  Eigen::VectorXd m_k_vec;  
  
  bool m_time_data_set;
  
  std::shared_ptr<V_Matrix_Construction> m_mtrx_builder;
    
};

#endif