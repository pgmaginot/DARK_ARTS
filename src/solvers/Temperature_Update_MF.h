#ifndef Temperature_Update_MF_h
#define Temperature_Update_MF_h

#include "V_Temperature_Update.h"

/** @file   Temperature_Update_MF.h
  *   @author pmaginot
  *   @brief Declare the Temperautre_Update class that will update a Temperature_Object given an Intensity_Object
 */

class Temperature_Update_MF : public V_Temperature_Update
{
public:
  Temperature_Update_MF(const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* const material, 
    const Angular_Quadrature& angular_quadrature, const Time_Stepper& time_stepper);
    
    
  ~Temperature_Update_MF(){}

  void update_temperature(const Intensity_Moment_Data& phi, 
    Temperature_Data& t_new, const Temperature_Data& t_star, const Temperature_Data& t_n,
    const K_Temperature& k_t, const int stage, const std::vector<double>& outside_rk,
    const double time, const double dt) override;    
private:  
   
  int m_n_groups;
  
  /// matrix to hold \f$ \sum_{g=1}^G{\mathbf{R}_{\sigma_{a,g}} \mathbf{D}_g } \f$
  Eigen::MatrixXd m_spectrum;
/* ****************************************************
*
*     Protected Functions
*
  **************************************************** */
  void calculate_local_matrices(const int cell , const Eigen::VectorXd& m_t_star ,
    const double dt, const double a_ii , const double time,
    const Intensity_Moment_Data& phi);
};

#endif