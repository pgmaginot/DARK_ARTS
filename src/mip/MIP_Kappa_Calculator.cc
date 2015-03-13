/** @file   Local_MIP_Assembler.cc
  *   @author pmaginot
  *   @brief Implement the (partially) virtual Diffusion_Matrix_Cretator class, 
  *     The goal of this clas is to construct MIP diffusion operator matrices in a generic format
  * concrete instantiations of this class will evaluate the following :
  *   \f$ \mathbf{R}_{\widetilde{\Sigma}_a} \f$, \f$ \mathbf{R}_{\widetilde{\Sigma}_s} \f$, 
  *  and $\widetilde{D}$ at cell edges and quadratureintegration points
 */

#include "MIP_Kappa_Calculator.h"

MIP_Kappa_Calculator::MIP_Kappa_Calculator(const int p_ord , const double z_mip)
:
  m_ord( double( p_ord ) ),
  m_z_mip( z_mip )
{

}

double MIP_Kappa_Calculator::calculate_interior_edge_kappa(const double dx_l, const double dx_r , const double d_l , const double d_r) const
{
  return std::max(0.25, m_z_mip/2.*(m_ord*(m_ord+1.))*(d_l/dx_l + d_r/dx_r) );
}
  
double MIP_Kappa_Calculator::calculate_boundary_kappa(const double dx , const double d) const
{
  return std::max(0.25, m_z_mip*(m_ord*(m_ord+1.))*(d/dx) );
}
  