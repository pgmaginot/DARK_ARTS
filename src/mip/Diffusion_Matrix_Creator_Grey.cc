/** @file   Diffusion_Matrix_Creator_Grey.cc
  *   @author pmaginot
  *   @brief Implement the Grey_Diffusion_Matrix_Cretator class, this accounts for the full SDRIK linearization
  *   this class will evaluate the following :
  *   \f$ \mathbf{R}_{\widetilde{\Sigma}_a} \f$, 
  *   \f$ \mathbf{R}_{\widetilde{\Sigma}_s} \f$, 
  *  and $\widetilde{D}$ at cell edges and quadrature integration points
 */

#include "Diffusion_Matrix_Creator_Grey.h"

Diffusion_Matrix_Creator_Grey::Diffusion_Matrix_Creator_Grey(const Fem_Quadrature& fem_quadrature, Materials& materials,
  const Temperature_Data& t_star)
:
  V_Diffusion_Matrix_Creator(fem_quadrature,materials,t_star)
{    

}

void Diffusion_Matrix_Creator_Grey::set_time_data( const double dt, const double time_stage, const double sdirk_a_of_stage )
{
  return;
}
  
void Diffusion_Matrix_Creator_Grey::calculate_pseudo_r_sig_a(void)
{
  return;
}
  

void Diffusion_Matrix_Creator_Grey::calculate_pseudo_r_sig_s(void)
{
  return;
}
  
void Diffusion_Matrix_Creator_Grey::evaluate_all_pseudo_d_coefficients(void)
{
  return;
}
    
