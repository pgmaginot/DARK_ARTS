#ifndef V_Temperature_Update_h
#define V_Temperature_Update_h

#include "Matrix_Construction_Exact.h"
#include "Matrix_Construction_Self_Lumping.h"
#include "Matrix_Construction_Trad_Lumping.h"

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

class V_Temperature_Update
{
public:
  V_Temperature_Update(const Fem_Quadrature& fem_quadrature, const Cell_Data& cell_data, Materials& material, 
    const Angular_Quadrature& angular_quadrature, const int n_stages);
    
    
  virtual ~V_Temperature_Update(){}

  virtual void update_temperature(const Intensity_Moment_Data& phi, Temperature_Data& t_star, 
    const Temperature_Data& t_n, const K_Temperature& k_t, const double damping, Err_Temperature& err_temperature) = 0;
    
  void set_time_data( const double dt, const double time_stage, const std::vector<double>& rk_a_of_stage_i, const int stage );

  virtual void calculate_k_t(const Temperature_Data& t_star, K_Temperature& k_t, const Intensity_Moment_Data& ard_phi) = 0;
protected:
  /// need to access material objects
  Materials& m_material;
  
  const double m_sn_w;
  
  /// number of DFEM points per cell
  const int m_np;
    
    /// number of cells
  const int m_n_cells;
  
  const int m_n_source_quad_pts;

  const Eigen::MatrixXd m_i_matrix;
  
  bool m_time_data_set;
  
  Eigen::VectorXd m_t_star_vec;
  Eigen::VectorXd m_t_old_vec;
  Eigen::VectorXd m_driving_source_vec;
  Eigen::VectorXd m_phi_vec;
  Eigen::VectorXd m_planck_vec;
  Eigen::VectorXd m_k_vec;
  
  Eigen::MatrixXd m_r_sig_a;
  Eigen::MatrixXd m_r_cv;
  Eigen::MatrixXd m_d_matrix;
    
  /// object that will build all local matrices for the temperature update
  std::shared_ptr<V_Matrix_Construction> m_mtrx_builder;
    
  /// the size of this vector is equal to the number of DFEM integration (quadrature points)
  std::vector<double> m_temp_mat_vec;  
  
  /// time stepping scheme data
  std::vector<double> m_rk_a;
  double m_dt = 0.;
  double m_time = 0.;
  int m_stage  = 0;
  
  bool check_eigen_variables_finite(void) const;

private:  

};

#endif