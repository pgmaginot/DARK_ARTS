#ifndef V_Temperature_Update_h
#define V_Temperature_Update_h

#include "V_Matrix_Construction.h"
#include "Matrix_Construction_Exact.h"
#include "Matrix_Construction_Self_Lumping.h"
#include "Matrix_Construction_Trad_Lumping.h"

#include "Temperature_Data.h"
#include "Intensity_Data.h"

#include "Fem_Quadrature.h"

#include "Cell_Data.h"

#include "Materials.h"

#include "Eigen/Dense"

#include <memory>

/** @file   Temperature_Update.h
  *   @author pmaginot
  *   @brief Declare the Temperautre_Update class that will update a Temperature_Object given an Intensity_Object
 */

class V_Temperature_Update
{
public:
  V_Temperature_Update(const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* material);
    
    
  ~V_Temperature_Update(){}

  void update_temperature(const Intensity_Data& intensity, Temperature_Data& t_star, Temperature_Data& t_n);    
private:  
    /// need to access material objects, save a ptr to avoid passing it all the time
  Materials* m_material;
  
  /// need to access cell data, save a ptr to skip passing it all the time
  Cell_Data* m_cell_data;
  
  const int m_np;
  
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
  Eigen::VectorXd m_planck;
  /// \f$ \vec{T}_n \f$
  Eigen::VectorXd m_t_old;
  /// \f$ \vec{T}^* \f$
  Eigen::VectorXd m_t_star;
  /// \f$ \vec{S}_T \f$
  Eigen::VectorXd m_driving_source;
  /// \f$ \vec{T}_i
  Eigen::VectorXd m_t_new;
  
  

/* ****************************************************
*
*     Protected Functions
*
  **************************************************** */
  
  
};

#endif