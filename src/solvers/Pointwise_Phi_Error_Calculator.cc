#include "Pointwise_Phi_Error_Calculator.h"

Pointwise_Phi_Error_Calculator::Pointwise_Phi_Error_Calculator( const Fem_Quadrature& fem_quadrature, 
    const Cell_Data& cell_data,
    const Angular_Quadrature& angular_quadrature)
    :
  V_Phi_Error_Calculator(fem_quadrature, cell_data, angular_quadrature),
  m_phi_small_val_vec(m_n_el , 0.)  
{
  
}

double Pointwise_Phi_Error_Calculator::calculate_phi_error_norm(const Intensity_Moment_Data& phi_new , const Intensity_Moment_Data& phi_old, const int iter)
{
  if(iter == 0)
  {
    phi_new.get_phi_norm(m_phi_small_val_vec);  
  }
  
  Eigen::VectorXd phi_new_vec = Eigen::VectorXd::Zero(m_n_el);
  Eigen::VectorXd phi_old_vec =  Eigen::VectorXd::Zero(m_n_el);
  double loc_err = 0.;
  
  for(int cell = 0; cell < m_n_cell; cell++)
  {
    for(int grp = 0; grp < m_n_groups ; grp++)
    {
      for(int l_mom = 0; l_mom < m_n_l_mom; l_mom++)
      {
        phi_new.get_cell_angle_integrated_intensity(cell,grp,l_mom,phi_new_vec);
        phi_old.get_cell_angle_integrated_intensity(cell,grp,l_mom,phi_old_vec);
        
        for(int el = 0; el < m_n_el ; el++)
        {
          if(fabs(phi_new_vec(el)) > m_phi_small_val_vec[grp])
          {
            /// not dividing by a near zero
            loc_err = std::max(loc_err , fabs( (phi_new_vec(el) - phi_old_vec(el) )/phi_new_vec(el) ) );
          }
          else
          {
            loc_err = std::max(loc_err, fabs( phi_new_vec(el) - phi_old_vec(el) ) );
          }
          
        }
      }
    }
  }
  return loc_err;
}