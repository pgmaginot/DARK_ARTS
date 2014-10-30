#ifndef Intensity_Update_MF_h
#define Intensity_Update_MF_h

#include "V_Intensity_Update.h"
#include "MF_ARD_Solver_FP_No_LMFGA.h"
// #include "MF_ARD_Solver_FP_LMFGA.h"
// #include "MF_ARD_Solver_Krylov_LMFGA.h"

/** @file   Intensity_Update_MF.h
  *   @author pmaginot
  *   @brief Calculate radiation intensity for a given temperature iterate for multifrequency radiation
 */

class Intensity_Update_MF : public V_Intensity_Update
{
public:
  Intensity_Update_MF(const Input_Reader& input_reader, 
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
    
  virtual ~Intensity_Update_MF(){}

  void update_intensity(Intensity_Moment_Data& phi) override;

private:
  std::shared_ptr<V_MF_ARD_Solver> m_ard_solver;
};

#endif