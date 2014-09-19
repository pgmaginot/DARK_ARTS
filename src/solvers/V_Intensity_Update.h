#ifndef V_Intensity_Update_h
#define V_Intensity_Update_h

#include "Input_Reader.h"
#include "Temperature_Data.h"
#include "Intensity_Moment_Data.h"
#include "Fem_Quadrature.h"
#include "Cell_Data.h"
#include "Materials.h"
#include "Angular_Quadrature.h"


#include <memory>

/** @file   V_Intensity_Update.h
  *   @author pmaginot
  *   @brief Calculate radiation intensity for a given temperature iterate
 */

class V_Intensity_Update
{
public:
  V_Intensity_Update(const Input_Reader& input_reader,
    const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* materials, 
    const Angular_Quadrature& angular_quadrature, const int n_stages);
    
  virtual ~V_Intensity_Update(){}

  virtual void update_intensity(const Temperature_Data& t_star, Intensity_Moment_Data& phi) = 0;
protected:
    double m_dt;
    int m_stage;
};

#endif