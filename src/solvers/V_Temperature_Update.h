#ifndef V_Temperature_Update_h
#define V_Temperature_Update_h

#include "V_Matrix_Construction.h"
#include "Matrix_Construction_Exact.h"
#include "Matrix_Construction_Self_Lumping.h"
#include "Matrix_Construction_Trad_Lumping.h"

#include "Temperature_Data.h"
#include "Intensity_Data.h"
#include "K_Temperature.h"

#include "Fem_Quadrature.h"

#include "Cell_Data.h"

#include "Materials.h"

#include "Angular_Quadrature.h"

#include "Eigen/Dense"

#include <memory>

/** @file   Temperature_Update.h
  *   @author pmaginot
  *   @brief Declare the Temperautre_Update class that will update a Temperature_Object given an Intensity_Object
 */

class V_Temperature_Update
{
public:
  V_Temperature_Update(const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* material, 
    const Angular_Quadrature& angular_quadrature);
    
    
  ~V_Temperature_Update(){}

  virtual void update_temperature(const Intensity_Data& intensity, 
    Temperature_Data& t_new, const Temperature_Data& t_star, const Temperature_Data& t_n,
    const K_Temperature& k_t, const int stage, const double time, const double dt) = 0;

  /// need to access material objects, save a ptr to avoid passing it all the time
  Materials* m_material;
  
 
  /// object that will build all local matrices for the temperature update
  std::shared_ptr<V_Matrix_Construction> m_mtrx_builder;
  
    /// \f$ \mathbf{R}_{\sigma_a} \f$
  Eigen::MatrixXd m_r_sig_a = Eigen::MatrixXd::Zero(m_np,m_np);
  /// \f$ \mathbf{R}_{C_v} \f$
  Eigen::MatrixXd m_r_cv = Eigen::MatrixXd::Zero(m_np,m_np);
  /// \f$ \mathbf{I} \f$
  Eigen::MatrixXd m_i_matrix = Eigen::MatrixXd::Zero(m_np,m_np);
  /// \f$ \mathbf{D} \f$
  Eigen::MatrixXd m_d_matrix = Eigen::MatrixXd::Zero(m_np,m_np);
  /// \f$ \mathbf{I} + 4\pi \Delta t a_{ii} \mathbf{R}_{C_v}^{-1} \mathbf{R}_{\sigma_a} \mathbf{D}  \f$
  Eigen::MatrixXd m_coeff_matrix = Eigen::MatrixXd::Zero(m_np,m_np);
  
  /// \f$ \vec{\widehat{B}} \f$
  Eigen::VectorXd m_planck = Eigen::VectorXd::Zero(m_np);
  /// \f$ \vec{T}_n \f$
  Eigen::VectorXd m_t_old = Eigen::VectorXd::Zero(m_np);
  /// \f$ \vec{T}^* \f$
  Eigen::VectorXd m_t_star = Eigen::VectorXd::Zero(m_np);
  /// \f$ \vec{S}_T \f$
  Eigen::VectorXd m_driving_source = Eigen::VectorXd::Zero(m_np);
  /// \f$ \vec{T}_i
  Eigen::VectorXd m_t_new = Eigen::VectorXd::Zero(m_np);
  
  std::vector<double> m_temp_mat_vec;
private:  
  /// need to access cell data, save a ptr to skip passing it all the time
  Cell_Data* m_cell_data;
  
  /// number of DFEM points per cell
  const int m_np;
  
  /// lumping type
  MATRIX_INTEGRATION m_matrix_type;
  
  /// number of cells
  const int m_n_cells;
  
  /// Sum of \f$ S_n \f$ weights
  const double m_sn_w;
  

  

/* ****************************************************
*
*     Protected Functions
*
  **************************************************** */
  
  
};

#endif