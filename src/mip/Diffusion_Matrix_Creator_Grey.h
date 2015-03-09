#ifndef Diffusion_Matrix_Creator_Grey_h
#define Diffusion_Matrix_Creator_Grey_h

#include "V_Diffusion_Matrix_Creator.h"

/** @file   Diffusion_Matrix_Creator_Grey.h
  *   @author pmaginot
  *   @brief Provide a class that constructs grey TRT MIP matrices with planck lienarization
 */
 

class Diffusion_Matrix_Creator_Grey: public V_Diffusion_Matrix_Creator
{
public:
  /// Only able to initialize if given an Input_Reader object
  /// constructor defined in Fem_Quadrature.cc
  Diffusion_Matrix_Creator_Grey(const Fem_Quadrature& fem_quadrature, Materials& materials,
    const Temperature_Data& t_star);
    
  virtual ~Diffusion_Matrix_Creator_Grey(){}
   
  /** MF needs to update M, r_cv, and spectrium \f$ \sum_{g=0}^G{\mathbf{R}_{\sigma_{a,g}} \mathbf{D}_g}   \f$
    Grey needs to update M, r_cv only
  */
  void set_time_data( const double dt, const double time_stage, const double sdirk_a_of_stage ) override;
    
  void calculate_pseudo_r_sig_a(void) override;
  
  void calculate_pseudo_r_sig_s(void) override;
  
  void evaluate_all_pseudo_d_coefficients(void) override;
  
protected:
  
};

#endif