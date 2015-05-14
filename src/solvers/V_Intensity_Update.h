#ifndef V_Intensity_Update_h
#define V_Intensity_Update_h

#include "WGRS_FP_Sweeps.h"
#include "WGRS_FP_DSA.h"

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
  V_Intensity_Update(const Input_Reader& input_reader, 
    const Fem_Quadrature& fem_quadrature, 
    const Cell_Data& cell_data,
    Materials& materials, 
    const Angular_Quadrature& angular_quadrature, 
    const int n_stages, 
    const Temperature_Data& t_old, 
    const Intensity_Data& i_old,
    const K_Temperature& kt, 
    K_Intensity& ki,
    const Temperature_Data& t_star ,
    std::vector<double>& phi_ref_norm);
    
  virtual ~V_Intensity_Update(){}
  
  virtual void kill_petsc_objects() = 0;

  virtual int update_intensity(Intensity_Moment_Data& phi, bool& update_sucess) = 0;
  
  void set_time_data( const double dt, const double time_stage, const std::vector<double>& rk_a_of_stage_i, const int stage );
  
  /// ard_phi won't change, but can't be labeled as const if it is to use the same transport sweep (and we use V_Solution_Saver as a interface template)
  void calculate_k_i(K_Intensity& k_i, Intensity_Moment_Data& ard_phi);
protected:
  std::shared_ptr<V_WGRS> m_within_group_radiation_solver;
};

#endif