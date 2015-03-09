/** @file   MG_WG_Diffusion_Ordering.cc
  *   @author pmaginot
  *   @brief ordering of MIP matrix for grey radiative transfer
 */

#include "MG_WG_Diffusion_Ordering.h"

MG_WG_Diffusion_Ordering::MG_WG_Diffusion_Ordering(const Cell_Data& cell_data, const Angular_Quadrature& angular_quadrature)
:
  V_Diffusion_Ordering(cell_data,angular_quadrature,angular_quadrature.get_number_of_groups() )
{

}

void MG_WG_Diffusion_Ordering::get_cell_and_group(const int block_i , const int mip_loop_number,  int& cell , int& group)
{
  cell = block_i;
  group = mip_loop_number;
  
  check_bounds(cell,group);    
}

