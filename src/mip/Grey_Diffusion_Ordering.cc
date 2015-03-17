/** @file   Grey_Diffusion_Ordering.cc
  *   @author pmaginot
  *   @brief ordering of MIP matrix for grey radiative transfer
 */

#include "Grey_Diffusion_Ordering.h"

Grey_Diffusion_Ordering::Grey_Diffusion_Ordering(const Cell_Data& cell_data, const Angular_Quadrature& angular_quadrature)
:
  V_Diffusion_Ordering(cell_data,angular_quadrature)
{

}

void Grey_Diffusion_Ordering::get_cell_and_group(const int block_i , int& cell , int& group)
{
  cell = block_i;
  group = 0;
  
  check_bounds(cell,group);    
}

