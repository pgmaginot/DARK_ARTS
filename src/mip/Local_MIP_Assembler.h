#ifndef Local_MIP_Assembler_h
#define Local_MIP_Assembler_h

#include <memory>
#include "Inputs_Allowed.h"
#include "Fem_Quadrature.h"
#include "MIP_Left_Boundary_Incident_Flux.h"
#include "MIP_Left_Boundary_Reflective.h"
#include "MIP_Right_Boundary_Incident_Flux.h"
#include "Eigen/Dense"

/** @file   Local_MIP_Assembler.h
  *   @author pmaginot
  *   @brief Assemble the local MIP matrices edge integrals of jump, average, cell integration operators
 */
 

class Local_MIP_Assembler
{
public:
  /// Only able to initialize if given an Input_Reader object
  /// constructor defined in Fem_Quadrature.cc
  Local_MIP_Assembler(const Fem_Quadrature& fem_quadrature, const Input_Reader& input_reader);
    
  virtual ~Local_MIP_Assembler(){}
    
  void calculate_left_boundary_matrices(const double kappa_12, const double kappa_32 ,   
    const double dx_c, const double dx_cp1, 
    const double d_c_l , const double d_c_r , const double d_cp1_l,
    const Eigen::MatrixXd& r_sig_a, const Eigen::MatrixXd& s_matrix,
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& cell_c, 
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& cell_cp1);
    
  void additive_update_left_boundary_rhs(Eigen::VectorXd& rhs);
    
  void calculate_interior_matrices(const double kappa_cm12, const double kappa_cp12 , 
    const double dx_cm1, const double dx_c, const double dx_cp1, 
    const double d_cm12_r, const double d_c_l , const double d_c_r , const double d_cp1_l,
    const Eigen::MatrixXd& r_sig_a, const Eigen::MatrixXd& s_matrix,
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& cell_cm1, 
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& cell_c, 
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& cell_cp1);
    
  void calculate_right_boundary_matrices(const double kappa_nm12, const double kappa_np12 , 
    const double dx_cm1, const double dx_c, 
    const double d_cm1_r, const double d_c_l , const double d_c_r ,
    const Eigen::MatrixXd& r_sig_a, const Eigen::MatrixXd& s_matrix,
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& cell_cm1, 
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& cell_c);
    
  void additive_update_right_boundary_rhs(Eigen::VectorXd& rhs);
    
protected:
  /** ****************************************************************
    * Variables that are initialzed in the constructor initialization list
    ****************************************************************
   */   
  const Eigen::RowVectorXd m_L;
  const Eigen::RowVectorXd m_Ls;
  const Eigen::RowVectorXd m_R;
  const Eigen::RowVectorXd m_Rs;
  
  const Eigen::MatrixXd m_Lt_R;
  const Eigen::MatrixXd m_Lst_R;
  const Eigen::MatrixXd m_Lt_Rs;
  
  const Eigen::MatrixXd m_Lt_L;
  const Eigen::MatrixXd m_Rt_R;
  const Eigen::MatrixXd m_Lst_L;
  const Eigen::MatrixXd m_Rst_R;
  const Eigen::MatrixXd m_Lt_Ls;
  const Eigen::MatrixXd m_Rt_Rs;
  
  const Eigen::MatrixXd m_Rt_L;
  const Eigen::MatrixXd m_Rst_L;
  const Eigen::MatrixXd m_Rt_Ls;
  
  std::shared_ptr<V_MIP_Left_Boundary> m_left_boundary;
  std::shared_ptr<V_MIP_Right_Boundary> m_right_boundary;
    
};

#endif