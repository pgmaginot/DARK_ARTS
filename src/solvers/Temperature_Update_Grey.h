#ifndef Temperature_Update_Grey_h
#define Temperature_Update_Grey__h

#include "V_Temperature_Update.h"

/** @file   Temperature_Update_Grey.h
  *   @author pmaginot
  *   @brief Declare the Temperautre_Update class that will update a Temperature_Object given an Intensity_Object
 */

class Temperature_Update_Grey: private V_Temperature_Update
{
public:
  Temperature_Update_Grey(const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* material,
    const Angular_Quadrature& angular_quadrature, const Time_Stepper& time_stepper);
    
    
  ~Temperature_Update_Grey(){}

  void update_temperature(const Intensity_Data& intensity, 
    Temperature_Data& t_new, const Temperature_Data& t_star, const Temperature_Data& t_n,
    const K_Temperature& k_t, const int stage, const std::vector<double>& outside_rk_a,
    const double time, const double dt) override;    
  
private:  
    

/* ****************************************************
*
*     Protected Functions
*
  **************************************************** */
  void calculate_local_matrices(const int cell_num, Eigen::VectorXd& t_eval,
    const double dt, const double a_ii, const double time);
    
  void get_planck_vector(const Eigen::VectorXd& t_eval);
  
  void get_planck_derivative_matrix(const Eigen::VectorXd& t_eval);
  
};

#endif