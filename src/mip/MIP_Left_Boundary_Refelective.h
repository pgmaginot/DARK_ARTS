#ifndef MIP_Left_Boundary_Reflective_h
#define MIP_Left_Boundary_Reflective_h

#include "V_MIP_Left_Boundary.h"
#include "Input_Reader.h"
#include <Eigen/Dense>
#include "Angular_Quadrature.h"
#include "Dark_Arts_Exception.h"

/** @file   V_MIP_Left_Boundary.h
  *   @author pmaginot
  *   @brief Provide a base class that implements MIP diffusion left boundary conditions
 */
 

class MIP_Left_Boundary_Reflective : public V_MIP_Left_Boundary
{
public:
  /// Only able to initialize if given an Input_Reader object
  /// constructor defined in Fem_Quadrature.cc
  MIP_Left_Boundary_Reflective(const Angular_Quadrature& angular_quadrature);
    
  virtual ~MIP_Left_Boundary_Reflective(){}
    
  void add_left_boundary_contributions(Eigen::MatrixXd& cell_c, Eigen::MatrixXd& cell_cp1, Eigen::VectorXd& rhs) override;
    
protected:
  /** ****************************************************************
    * Variables that are initialzed in the constructor initialization list
    ****************************************************************
   */   

};

#endif