#ifndef V_MIP_Left_Boundary_h
#define V_MIP_Left_Boundary_h

#include <Eigen/Dense>
#include "Fem_Quadrature.h"
#include "Dark_Arts_Exception.h"

/** @file   V_MIP_Left_Boundary.h
  *   @author pmaginot
  *   @brief Provide a base class that implements MIP diffusion left boundary conditions
 */
 

class V_MIP_Left_Boundary
{
public:
  /// Only able to initialize if given an Input_Reader object
  /// constructor defined in Fem_Quadrature.cc
  V_MIP_Left_Boundary(const Fem_Quadrature& fem_quadrature);
    
  virtual ~V_MIP_Left_Boundary(){}
    
  virtual void add_left_boundary_contributions(const double kappa_12, const double kappa_32 ,   
    const double d_1_l , const double d_1_r , const double d_2_l,
    const double dx_1, const double dx_2, 
    Eigen::MatrixXd& cell_c, Eigen::MatrixXd& cell_cp1, Eigen::VectorXd& rhs) = 0;
  
  void set_time(const double time);
  void set_incoming_current(const double incoming_current);
protected:
  double m_time;
  double m_incoming_current;
  /** ****************************************************************
    * Variables that are initialzed in the constructor initialization list
    ****************************************************************
   */   
  const Eigen::RowVectorXd m_L;
  const Eigen::RowVectorXd m_R;
  const Eigen::RowVectorXd m_Rs;
  const Eigen::RowVectorXd m_Ls;

};

#endif