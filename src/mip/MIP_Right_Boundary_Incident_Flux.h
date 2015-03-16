#ifndef MIP_Right_Boundary_Incident_Flux_h
#define MIP_Right_Boundary_Incident_Flux_h

#include "Input_Reader.h"
#include <Eigen/Dense>
#include "Angular_Quadrature.h"
#include "Dark_Arts_Exception.h"
#include "V_MIP_Right_Boundary.h"

/** @file   MIP_Right_Boundary_Incident_Flux.h
  *   @author pmaginot
  *   @brief Provide a base class that implements MIP diffusion left boundary conditions
 */
 

class MIP_Right_Boundary_Incident_Flux : public V_MIP_Right_Boundary
{
public:
  /// Only able to initialize if given an Input_Reader object
  /// constructor defined in Fem_Quadrature.cc
  MIP_Right_Boundary_Incident_Flux(const Fem_Quadrature& fem_quadrature);
    
  virtual ~MIP_Right_Boundary_Incident_Flux(){}
    
  void add_right_boundary_contributions(const double kappa_nm12, const double kappa_np12 ,   
    const double d_cm1_r, const double d_c_l , const double d_c_r , 
    const double dx_cm1, const double dx_c, 
    Eigen::MatrixXd& cell_cm1, Eigen::MatrixXd& cell_c) override;
    
    void add_right_boundary_rhs_contributions(Eigen::VectorXd& rhs) override;
    
protected:
  /** ****************************************************************
    * Variables that are initialzed in the constructor initialization list
    ****************************************************************
   */   
  const Eigen::MatrixXd m_Lt_R;
  const Eigen::MatrixXd m_Lt_L;
  
  const Eigen::MatrixXd m_Lst_R;
  const Eigen::MatrixXd m_Lst_L;
  
  const Eigen::MatrixXd m_Lt_Rs;
  const Eigen::MatrixXd m_Lt_Ls;
  
  const Eigen::MatrixXd m_Rt_R;
  const Eigen::MatrixXd m_Rst_R;
  const Eigen::MatrixXd m_Rt_Rs;
};

#endif