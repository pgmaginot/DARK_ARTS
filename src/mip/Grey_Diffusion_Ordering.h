#ifndef Grey_Diffusion_Ordering_h
#define Grey_Diffusion_Ordering_h

#include "V_Diffusion_Ordering.h"

/** @file   Grey_Diffusion_Ordering.h
  *   @author pmaginot
  *   @brief Ordering for grey diffusion operator
 */
 

class Grey_Diffusion_Ordering : public V_Diffusion_Ordering
{
public:
  /// Only able to initialize if given an Input_Reader object
  /// constructor defined in Fem_Quadrature.cc
  Grey_Diffusion_Ordering(const Cell_Data& cell_data, const Angular_Quadrature& angular_quadrature);
    
  virtual ~Grey_Diffusion_Ordering(){}
    
  void get_cell_and_group(const int block_i , int& cell , int& group) override;
  
protected:

};

#endif