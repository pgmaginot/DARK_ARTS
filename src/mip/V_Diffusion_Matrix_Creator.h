#ifndef V_Diffusion_Matrix_Creator_h
#define V_Diffusion_Matrix_Creator_h

#include "Inputs_Allowed.h"
#include "Fem_Quadrature.h"
#include "Eigen/Dense"
#include "Matrix_Construction_Exact.h"
#include "Matrix_Construction_Self_Lumping.h"
#include "Matrix_Construction_Trad_Lumping.h"

#include "Cell_Data.h"
#include "Materials.h"
#include "Temperature_Data.h"

/** @file   V_Diffusion_Matrix_Creator.h
  *   @author pmaginot
  *   @brief Provide a base class that constructs MIP matrices, concrete instantiations 
 */
 

class V_Diffusion_Matrix_Creator
{
public:
  /// Only able to initialize if given an Input_Reader object
  /// constructor defined in Fem_Quadrature.cc
  V_Diffusion_Matrix_Creator(const Fem_Quadrature& fem_quadrature, Materials& materials,
  const Temperature_Data& t_star);
    
  virtual ~V_Diffusion_Matrix_Creator(){}
   
  /** MF needs to update M, r_cv, and spectrium \f$ \sum_{g=0}^G{\mathbf{R}_{\sigma_{a,g}} \mathbf{D}_g}   \f$
    Grey needs to update M, r_cv only
  */
  virtual void set_time_data( const double dt, const double time_stage, const double sdirk_a_of_stage ) = 0;
    
  virtual void calculate_pseudo_r_sig_a(void) = 0;
  
  virtual void calculate_pseudo_r_sig_s(void) = 0;
  
  virtual void evaluate_all_pseudo_d_coefficients(void) = 0;
  
protected:
  /** ****************************************************************
    * Variables that are initialzed in the constructor initialization list
    ****************************************************************
   */   
  const int m_np;
  
  const int m_n_integration_pts;
  
  Eigen::MatrixXd m_r_sig_s;
  
  Eigen::MatrixXd m_r_sig_a; 
  
  double m_d_r_cm1;
  double m_d_l_c;
  double m_d_r_c;
  double m_d_l_cp1;
  
  std::vector<double> m_d_at_integration_pts;
  
  const Temperature_Data& m_t_eval;
      
  /// builder/lumper of reaction matrices
  std::shared_ptr<V_Matrix_Construction> m_mtrx_builder;
};

#endif