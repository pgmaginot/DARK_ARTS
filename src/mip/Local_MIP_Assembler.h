#ifndef Local_MIP_Assembler_h
#define Local_MIP_Assembler_h

#include "Inputs_Allowed.h"
#include "Fem_Quadrature.h"
#include "Eigen/Dense"
// #include "Diffusion_Matrix_Creator_Grey.h"
/** @file   Local_MIP_Assembler.h
  *   @author pmaginot
  *   @brief Assemble the local MIP matrices edge integrals of jump, average, cell integration operators
 */
 

class Local_MIP_Assembler
{
public:
  /// Only able to initialize if given an Input_Reader object
  /// constructor defined in Fem_Quadrature.cc
  Local_MIP_Assembler(const Fem_Quadrature& fem_quadrature);
    
  virtual ~Local_MIP_Assembler(){}
    
  void calculate_left_boundary_matrices(const double kappa_cm12, const double kappa_cp12 , 
    const double d_cm12_r, const double d_c_l , const double d_c_r , const double d_cp1_l,
    const Eigen::MatrixXd& r_sig_a, const Eigen::MatrixXd& s_matrix,
    Eigen::MatrixXd& cell_c, Eigen::MatrixXd& cell_cp1, Eigen::VectorXd& rhs);
    
  void calculate_interior_matrices(const double kappa_cm12, const double kappa_cp12 , 
    const double d_cm12_r, const double d_c_l , const double d_c_r , const double d_cp1_l,
    const Eigen::MatrixXd& r_sig_a, const Eigen::MatrixXd& s_matrix,
    Eigen::MatrixXd& cell_cm1, Eigen::MatrixXd& cell_c, Eigen::MatrixXd& cell_cp1);
    
  void calculate_right_boundary_matrices(const double kappa_cm12, const double kappa_cp12 , 
    const double d_cm12_r, const double d_c_l , const double d_c_r , const double d_cp1_l,
    const Eigen::MatrixXd& r_sig_a, const Eigen::MatrixXd& s_matrix,
    Eigen::MatrixXd& cell_cm1, Eigen::MatrixXd& cell_c);
    
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
  
  Eigen::MatrixXd cell_cm1;
  Eigen::MatrixXd cell_c;
  Eigen::MatrixXd cell_cp1;
  
};

#endif