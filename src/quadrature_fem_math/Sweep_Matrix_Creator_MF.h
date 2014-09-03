#ifndef Sweep_Matrix_Creator_MF_h
#define Sweep_Matrix_Creator_MF_h

#include "V_Sweep_Matrix_Creator.h"

/** @file   Sweep_Matrix_Creator_MF.h
  *   @author pmaginot
  *   @brief Provide a concrete class that matrices and sources for MF radiative transfer sweeps
 */

class Sweep_Matrix_Creator_MF : private V_Sweep_Matrix_Creator
{
public:
  /// Only able to initialize if given an Input_Reader object
  /// constructor defined in Fem_Quadrature.cc
  Sweep_Matrix_Creator_MF(const Fem_Quadrature& fem_quadrature);
  virtual ~Sweep_Matrix_Creator_MF(){}
  
  void construct_r_sig_t(const int cell, const int grp, Eigen::MatrixXd& r_sig_t) override;
  
  void construct_r_sig_s(const int cell, const int grp, const int l_mom, Eigen::MatrixXd& r_sig_s) override;
  
  void construct_s_i(const int cell,const int grp, const int l_mom, Eigen::VectorXd& s_i) override;
   
  
private:

};

#endif