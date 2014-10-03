#ifndef V_Intensity_Update_h
#define V_Intensity_Update_h

// #include "Input_Reader.h"
// #include "Temperature_Data.h"
// #include "Intensity_Moment_Data.h"
// #include "Fem_Quadrature.h"
// #include "Cell_Data.h"
// #include "Materials.h"
// #include "Angular_Quadrature.h"

// #include "V_WGRS.h"
#include "WGRS_FP_Sweeps.h"

// #include <memory>

/** @file   V_Intensity_Update.h
  *   @author pmaginot
  *   @brief Calculate radiation intensity for a given temperature iterate
  *    For a given temperature iterate, t_star, calculate an Intensity_Moment_Data object, phi, that can be used by the temperature update equation
  *  this process is different between grey and MF problems, thus update_intensity is a virtual function defined by the concrete instances of 
  *  of V_Intensity_Update, Intensity_Update_Grey and Intensity_Update_MF, respectively.
 */

class V_Intensity_Update
{
public:
  V_Intensity_Update(const Input_Reader& input_reader, const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* materials, 
    Angular_Quadrature& angular_quadrature, const int n_stages,  
    const Temperature_Data* const t_old, const Temperature_Data* const t_star, 
    const Intensity_Data* const i_old,
    const K_Temperature* const kt, const K_Intensity* const ki);
    
  virtual ~V_Intensity_Update(){}

  virtual void update_intensity(const Temperature_Data* const t_star, Intensity_Moment_Data& phi) = 0;
protected:
  double m_dt;
  int m_stage;
  
  std::shared_ptr<V_WGRS> m_within_group_radiation_solver;
};

#endif