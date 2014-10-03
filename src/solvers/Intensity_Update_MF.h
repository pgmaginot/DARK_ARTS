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
  Intensity_Update_MF(const Input_Reader& input_reader, const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* material, 
    Angular_Quadrature& angular_quadrature, const int n_stages, 
    const Temperature_Data* const t_old, const Temperature_Data* const t_star, 
    const Intensity_Data* const i_old,
    const K_Temperature* const kt, const K_Intensity* const ki);
    
  ~Intensity_Update_MF(){}

  void update_intensity(const Temperature_Data* const t_star, Intensity_Moment_Data& phi) override;

};

#endif