/** @file   V_Diffusion_Matrix_Creator.cc
  *   @author pmaginot
  *   @brief Implement the (partially) virtual Diffusion_Matrix_Cretator class, 
  *     The goal of this clas is to construct MIP diffusion operator matrices in a generic format
  * concrete instantiations of this class will evaluate the following :
  *   \f$ \mathbf{R}_{\widetilde{\Sigma}_a} \f$, \f$ \mathbf{R}_{\widetilde{\Sigma}_s} \f$, 
  *  and $\widetilde{D}$ at cell edges and quadratureintegration points
 */

#include "V_Diffusion_Ordering.h"

V_Diffusion_Ordering::V_Diffusion_Ordering(const Cell_Data& cell_data, const Angular_Quadrature& angular_quadrature)
:
  m_n_groups(angular_quadrature.get_number_of_groups() ),
  m_n_cells(cell_data.get_total_number_of_cells() )
{

}

void V_Diffusion_Ordering::check_bounds(const int cell_num , const int group_num)
{
  if( (cell_num >= m_n_cells) || (cell_num < 0) )
    throw Dark_Arts_Exception( TIME_MARCHER, "Diffusion ordering resulting in non-logical cell number");
    
  if( (group_num >= m_n_groups) || (group_num < 0) )
    throw Dark_Arts_Exception( TIME_MARCHER, "Diffusion ordering resulting in non-logical group number");
    
}

