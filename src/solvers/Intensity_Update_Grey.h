#ifndef Intensity_Update_Grey_h
#define Intensity_Update_Grey_h


#include "V_Intensity_Update.h"

/** @file   Intensity_Update_Grey.h
  *   @author pmaginot
  *   @brief Calculate radiation intensity for a given temperature iterate, for grey radiation
 */

class Intensity_Update_Grey : public V_Intensity_Update
{
public:
  Intensity_Update_Grey(const Input_Reader& input_reader, const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* materials, 
    const Angular_Quadrature& angular_quadrature, const int n_stages);
    
  ~Intensity_Update_Grey(){}

  void update_intensity(const Temperature_Data& t_star, Intensity_Moment_Data& phi) override;

};

#endif