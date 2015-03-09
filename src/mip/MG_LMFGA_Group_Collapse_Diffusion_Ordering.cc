/** @file   MG_LMFGA_Group_Collapse_Diffusion_Ordering.cc
  *   @author pmaginot
  *   @brief ordering of MIP matrix for grey radiative transfer
 */

#include "MG_LMFGA_Group_Collapse_Diffusion_Ordering.h"

MG_LMFGA_Group_Collapse_Diffusion_Ordering::MG_LMFGA_Group_Collapse_Diffusion_Ordering(const Cell_Data& cell_data, const Angular_Quadrature& angular_quadrature)
:
  V_Diffusion_Ordering(cell_data,angular_quadrature,1)
{

}

void MG_LMFGA_Group_Collapse_Diffusion_Ordering::get_cell_and_group(const int block_i , const int mip_loop_number,  int& cell , int& group)
{
  cell = block_i;
  group = 0;
  
  if( ~(mip_loop_number == 0) )
    throw Dark_Arts_Exception(TIME_MARCHER , "Trying to use non zero mip_loop_number for MG_LMFGA_Group_Collapse ");
  
  check_bounds(cell,group);    
}

