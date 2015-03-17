#ifndef V_MIP_Right_Boundary_h
#define V_MIP_Right_Boundary_h

#include "Fem_Quadrature.h"
#include <Eigen/Dense>
#include "Dark_Arts_Exception.h"

/** @file   V_MIP_Left_Boundary.h
  *   @author pmaginot
  *   @brief Provide a base class that implements MIP diffusion left boundary conditions
 */
 

class V_MIP_Right_Boundary
{
public:
  /// Only able to initialize if given an Input_Reader object
  /// constructor defined in Fem_Quadrature.cc
  V_MIP_Right_Boundary(const Fem_Quadrature& fem_quadrature);
    
  virtual ~V_MIP_Right_Boundary(){}
    
  virtual void add_right_boundary_contributions(const double kappa_nm12, const double kappa_np12 ,   
    const double d_cm1_r, const double d_c_l , const double d_c_r , 
    const double dx_cm1, const double dx_c, 
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& cell_cm1, 
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& cell_c) = 0;
    
  virtual void add_right_boundary_rhs_contributions(Eigen::VectorXd& rhs) = 0;
    
  void set_time(const double time);
  
  void set_incoming_current(const double current);
protected:
  /** ****************************************************************
    * Variables that are initialzed in the constructor initialization list
    ****************************************************************
   */   
  double m_time;
  double m_incoming_current;
   
  const Eigen::RowVectorXd m_L;
  const Eigen::RowVectorXd m_R;
  const Eigen::RowVectorXd m_Ls;
  const Eigen::RowVectorXd m_Rs;

};

#endif