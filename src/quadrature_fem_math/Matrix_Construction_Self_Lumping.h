#ifndef Matrix_Construction_Self_Lumping_h
#define Matrix_Construction_Self_Lumping_h

#include "V_Matrix_Construction.h"

/** @file   Matrix_Construction_Self_Lumping.h
  *   @author pmaginot
  *   @brief Implement the SELF_LUMPING matrix integration techniques
 */

class Matrix_Construction_Self_Lumping: public V_Matrix_Construction
{
public:
  /// Only able to initialize if given an Input_Reader object
  /// constructor defined in Fem_Quadrature.cc
  Matrix_Construction_Self_Lumping(const Fem_Quadrature& fem_quadrature,  Materials& materials);
  virtual ~Matrix_Construction_Self_Lumping(){}
  
  void construct_dimensionless_mass_matrix(Eigen::MatrixXd& mass) override;
  
protected:  
  
  void construct_reaction_matrix(Eigen::MatrixXd& rx_mat, std::vector<double>& xs) override;
};

#endif