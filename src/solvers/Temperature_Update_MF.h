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
  Temperature_Update_MF(const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* material, 
    const Angular_Quadrature& angular_quadrature);
    
    
  ~Temperature_Update_MF(){}

  void update_temperature(const Intensity_Data& intensity, 
    Temperature_Data& t_new, const Temperature_Data& t_star, const Temperature_Data& t_n,
    const K_Temperature& k_t, const int stage, const double time, const double dt) override;    
private:  
   

/* ****************************************************
*
*     Protected Functions
*
  **************************************************** */
  
  
};

#endif