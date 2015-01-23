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
  Intensity_Update_Grey(const Input_Reader& input_reader, 
  const Fem_Quadrature& fem_quadrature, 
  const Cell_Data& cell_data, 
  Materials& materials, 
  const Angular_Quadrature& angular_quadrature, 
  const int n_stages, 
  const Temperature_Data& t_old, 
  const Intensity_Data& i_old,
  const K_Temperature& kt,
  K_Intensity& ki,
  const Temperature_Data& t_star,
  std::vector<double>& phi_ref_norm);
    
  virtual ~Intensity_Update_Grey(){}

  int update_intensity(Intensity_Moment_Data& phi) override;

};

#endif