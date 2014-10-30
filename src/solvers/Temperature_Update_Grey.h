#ifndef Temperature_Update_Grey_h
#define Temperature_Update_Grey_h

#include "V_Temperature_Update.h"

/** @file   Temperature_Update_Grey.h
  *   @author pmaginot
  *   @brief Declare the Temperautre_Update class that will update a Temperature_Object given an Intensity_Object
 */

class Temperature_Update_Grey : public V_Temperature_Update
{
public:
  Temperature_Update_Grey(const Fem_Quadrature& fem_quadrature, const Cell_Data& cell_data, Materials& material, 
    const Angular_Quadrature& angular_quadrature, const int n_stages);
    
    
  virtual ~Temperature_Update_Grey(){}

  void update_temperature(const Intensity_Moment_Data& phi, Temperature_Data& t_star, 
    const Temperature_Data& t_n, const K_Temperature& k_t, const double damping, Err_Temperature& err_temperature) override;    
  
  void calculate_k_t(const Temperature_Data& t_star, K_Temperature& k_t, const Intensity_Moment_Data& ard_phi) override;
private:  
    

/* ****************************************************
*
*     Protected Functions
*
  **************************************************** */
  /// just need the cell number to access the material ptr correctly
  /// all local variables (cell_num, t_star_vec, dt, time, a_ii, are already saved in the V_Temperature_Update object
  void calculate_local_matrices(void);
  
};

#endif