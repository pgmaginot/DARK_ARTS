#ifndef MIP_Left_Boundary_Reflective_h
#define MIP_Left_Boundary_Reflective_h

#include "V_MIP_Left_Boundary.h"

/** @file   MIP_Left_Boundary_Reflective.h
  *   @author pmaginot
  *   @brief Provide a base class that implements MIP diffusion left boundary conditions
 */
 

class MIP_Left_Boundary_Reflective : public V_MIP_Left_Boundary
{
public:
  /// Only able to initialize if given an Input_Reader object
  /// constructor defined in Fem_Quadrature.cc
  MIP_Left_Boundary_Reflective(const Fem_Quadrature& fem_quadrature);
    
  virtual ~MIP_Left_Boundary_Reflective(){}
    
  void add_left_boundary_contributions(const double kappa_12, const double kappa_32 ,   
    const double d_1_l , const double d_1_r , const double d_2_l,
    const double dx_1, const double dx_2, 
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& cell_c, 
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& cell_cp1) override;
    
  void add_left_boundary_rhs_contributions(Eigen::VectorXd& rhs) override;
    
protected:
  /** ****************************************************************
    * Variables that are initialzed in the constructor initialization list
    ****************************************************************
   */   
  const Eigen::MatrixXd m_Rt_L;
  const Eigen::MatrixXd m_Rt_R;
  const Eigen::MatrixXd m_Rst_L;
  const Eigen::MatrixXd m_Rst_R;
  const Eigen::MatrixXd m_Rt_Rs;
  const Eigen::MatrixXd m_Rt_Ls;
};

#endif