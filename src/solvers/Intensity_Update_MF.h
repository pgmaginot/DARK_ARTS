#ifndef Intensity_Update_MF_h
#define Intensity_Update_MF_h

#include "V_Intensity_Update.h"


/** @file   Intensity_Update_MF.h
  *   @author pmaginot
  *   @brief Calculate radiation intensity for a given temperature iterate for multifrequency radiation
 */

class Intensity_Update_MF : public V_Intensity_Update
{
public:
  Intensity_Update_MF(const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* material, 
    const Angular_Quadrature& angular_quadrature, const int n_stages);
    
  ~Intensity_Update_MF(){}

  void update_intensity(const Temperature_Data& t_star, Intensity_Moment_Data& phi) override;

};

#endif