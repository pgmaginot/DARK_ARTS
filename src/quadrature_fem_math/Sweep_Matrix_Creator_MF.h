#ifndef Sweep_Matrix_Creator_MF_h
#define Sweep_Matrix_Creator_MF_h

#include "V_Sweep_Matrix_Creator.h"

/** @file   Sweep_Matrix_Creator_MF.h
  *   @author pmaginot
  *   @brief Provide a concrete class that matrices and sources for MF radiative transfer sweeps
 */

class Sweep_Matrix_Creator_MF : public V_Sweep_Matrix_Creator
{
public:
  /// Only able to initialize if given an Input_Reader object
  /// constructor defined in Fem_Quadrature.cc
  Sweep_Matrix_Creator_MF(const Fem_Quadrature& fem_quadrature, Materials* const materials,
    Cell_Data* const cell_data, const int n_stages);
  virtual ~Sweep_Matrix_Creator_MF(){}
  

   
  
private:
  void construct_r_sig_t(void) override;
  
  void construct_r_sig_s(void) override;
  
  void construct_s_i(void) override;
};

#endif