#ifndef Matrix_Construction_Trad_Lumping_h
#define Matrix_Construction_Trad_Lumping_h

#include "V_Matrix_Construction.h"

/** @file   Matrix_Construction_Trad_Lumping.h
  *   @author pmaginot
  *   @brief Implement the Exact matrix integration techniques (exact to the accuracy of the quadrature at lest)
 */

class Matrix_Construction_Trad_Lumping : public V_Matrix_Construction
{
public:
  Matrix_Construction_Trad_Lumping(const Fem_Quadrature& fem_quadrature, Materials& materials);
  virtual ~Matrix_Construction_Trad_Lumping(){}
    
  void construct_dimensionless_mass_matrix(Eigen::MatrixXd& mass) override;

protected:  
  
  void construct_reaction_matrix(Eigen::MatrixXd& rx_mat, std::vector<double>& xs) override;
  
};

#endif