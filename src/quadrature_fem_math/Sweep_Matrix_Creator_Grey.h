#ifndef Sweep_Matrix_Creator_Grey_h
#define Sweep_Matrix_Creator_Grey_h

#include "V_Sweep_Matrix_Creator.h"

/** @file   Sweep_Matrix_Creator_Grey.h
  *   @author pmaginot
  *   @brief Provide a base class that constructs mass, reaction, and gradient matrices, as well as upwind vectors
  *     Concrete cases for  SELF_LUMPING , TRAD_LUMPING , and EXACT integration techniques
 */

class Sweep_Matrix_Creator_Grey : private V_Sweep_Matrix_Creator
{
public:
  /// Only able to initialize if given an Input_Reader object
  /// constructor defined in Fem_Quadrature.cc
  Sweep_Matrix_Creator_Grey(const Fem_Quadrature& fem_quadrature);
  virtual ~Sweep_Matrix_Creator_Grey(){}
  
  void construct_r_sig_t(const int cell, const int grp, Eigen::MatrixXd& r_sig_t) override;
  
  void construct_r_sig_s(const int cell, const int grp, const int l_mom, Eigen::MatrixXd& r_sig_s) override;
  
  void construct_s_i(const int cell,const int grp, const int l_mom, Eigen::VectorXd& s_i) override;
   
  
private:

};

#endif