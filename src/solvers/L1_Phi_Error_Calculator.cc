#include "L1_Phi_Error_Calculator.h"

L1_Phi_Error_Calculator::L1_Phi_Error_Calculator( const Fem_Quadrature& fem_quadrature, 
    const Cell_Data& cell_data,
    const Angular_Quadrature& angular_quadrature)
    :
  V_Phi_Error_Calculator(fem_quadrature, cell_data , angular_quadrature)
{

}


double L1_Phi_Error_Calculator::calculate_phi_error_norm(const Intensity_Moment_Data& phi_new , const Intensity_Moment_Data& phi_old, const int iter)
{
  double err = 0.;
  Eigen::VectorXd local_phi_old = Eigen::VectorXd::Zero(m_n_el);
  Eigen::VectorXd local_phi_new = Eigen::VectorXd::Zero(m_n_el);
  
  double phi_normalizer = 0.;
  for(int cell = 0; cell < m_n_cell; cell++)
  {
    for(int grp = 0; grp < m_n_groups ; grp++)
    {
      for(int l_mom = 0; l_mom < m_n_l_mom ; l_mom++)
      {
        phi_new.get_cell_angle_integrated_intensity(cell,grp,l_mom,local_phi_new);
        phi_old.get_cell_angle_integrated_intensity(cell,grp,l_mom,local_phi_old);
        
        for(int el = 0; el < m_n_el ; el++)
        {
          phi_normalizer = std::max(phi_normalizer , fabs(local_phi_new(el) ));
          err = std::max(err , fabs(local_phi_new(el) - local_phi_old(el)));         
        }
      }
    }
  }
  err /= phi_normalizer;
  return err;
}