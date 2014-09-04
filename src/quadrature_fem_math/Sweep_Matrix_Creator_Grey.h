#ifndef Sweep_Matrix_Creator_Grey_h
#define Sweep_Matrix_Creator_Grey_h

#include "V_Sweep_Matrix_Creator.h"

/** @file   Sweep_Matrix_Creator_Grey.h
  *   @author pmaginot
  *   @brief Provide a base class that constructs mass, reaction, and gradient matrices, as well as upwind vectors
  *     Concrete cases for  SELF_LUMPING , TRAD_LUMPING , and EXACT integration techniques
 */

class Sweep_Matrix_Creator_Grey : public V_Sweep_Matrix_Creator
{
public:
  /// Only able to initialize if given an Input_Reader object
  /// constructor defined in Fem_Quadrature.cc
  Sweep_Matrix_Creator_Grey(const Fem_Quadrature& fem_quadrature, Materials* const materials,
    Cell_Data* const cell_data, const int n_stages);
  virtual ~Sweep_Matrix_Creator_Grey(){}
  

   
  
private:
  void construct_r_sig_t(void);
  
  void construct_r_sig_s(void);
  
  void construct_s_i(void);
};

#endif