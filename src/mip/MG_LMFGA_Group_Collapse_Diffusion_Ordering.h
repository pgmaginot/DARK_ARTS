#ifndef MG_LMFGA_Group_Collapse_Diffusion_Ordering_h
#define MG_LMFGA_Group_Collapse_Diffusion_Ordering_h

#include "V_Diffusion_Ordering.h"

/** @file   MG_LMFGA_Group_Collapse_Diffusion_Ordering.h
  *   @author pmaginot
  *   @brief Ordering for grey diffusion operator
 */
 

class MG_LMFGA_Group_Collapse_Diffusion_Ordering : public V_Diffusion_Ordering
{
public:
  /// Only able to initialize if given an Input_Reader object
  /// constructor defined in Fem_Quadrature.cc
  MG_LMFGA_Group_Collapse_Diffusion_Ordering(const Cell_Data& cell_data, const Angular_Quadrature& angular_quadrature);
    
  virtual ~MG_LMFGA_Group_Collapse_Diffusion_Ordering(){}
    
  virtual void get_cell_and_group(const int block_i , const int mip_loop_number,  int& cell , int& group) override;
  
protected:

};

#endif