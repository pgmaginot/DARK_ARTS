#ifndef Sweep_Fixed_Source_Linearization_h
#define Sweep_Fixed_Source_Linearization_h

#include "V_Sweep_Fixed_Source.h"
#include "V_Sweep_Matrix_Creator.h"
#include <memory>

/** @file   Sweep_Fixed_Source_Linearization.h
  *   @author pmaginot
  *   @brief allow for access to either the Grey of MF sweep_matrix_creator_object
 */

class Sweep_Fixed_Source_Linearization : public V_Sweep_Fixed_Source
{
public:
  Sweep_Fixed_Source_Linearization(const Fem_Quadrature& fem_quadrature, std::shared_ptr<V_Sweep_Matrix_Creator> sweep_matrix_ptr);
    
  virtual ~Sweep_Fixed_Source_Linearization(){}

  void get_source(Eigen::VectorXd& source_vec) override;
protected:
  std::shared_ptr<V_Sweep_Matrix_Creator> m_sweep_matrix_ptr;
};

#endif