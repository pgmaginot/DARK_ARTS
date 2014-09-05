#ifndef V_Sweep_Matrix_Creator_h
#define V_Sweep_Matrix_Creator_h

#include "Inputs_Allowed.h"
#include "Fem_Quadrature.h"
#include "Eigen/Dense"
#include "V_Matrix_Construction.h"
#include "Matrix_Construction_Exact.h"
#include "Matrix_Construction_Self_Lumping.h"
#include "Matrix_Construction_Trad_Lumping.h"

#include "Cell_Data.h"
#include "Materials.h"
#include "Temperature_Data.h"
#include "K_Intensity.h"
#include "K_Temperature.h"

#include "Intensity_Data.h"

#include <vector>
#include <memory>

/** @file   V_Matrix_Construction.h
  *   @author pmaginot
  *   @brief Provide a base class that constructs mass, reaction, and gradient matrices, as well as upwind vectors
  *     Concrete cases for  SELF_LUMPING , TRAD_LUMPING , and EXACT integration techniques
 */

class V_Sweep_Matrix_Creator
{
public:
  /// Only able to initialize if given an Input_Reader object
  /// constructor defined in Fem_Quadrature.cc
  V_Sweep_Matrix_Creator(const Fem_Quadrature& fem_quadrature, Materials* const materials,
    Cell_Data* const cell_data, const int n_stages);
    
  virtual ~V_Sweep_Matrix_Creator(){}
  
  void construct_l_matrix(const double mu, Eigen::MatrixXd& l_matrix);
  
  void construct_f_vector(const double mu, Eigen::VectorXd& f_vector);
  
  /// data that changes every time we update radiation intensity
  void set_thermal_iteration_data(const Temperature_Data* t_eval, const Temperature_Data* t_old, 
    const K_Temperature* kt );  

  void set_intensity_iteration_data(const Intensity_Data* i_old, const K_Intensity* ki);
  
  /// data that changes once per stage (per time step)
  
  void set_stage_data(const int stage, const std::vector<double>& rk_a, const double time);

  void set_timestep_data(const double dt);
  
  /** MF needs to update M, r_cv, and spectrium \f$ \sum_{g=0}^G{\mathbf{R}_{\sigma_{a,g}} \mathbf{D}_g}   \f$
    Grey needs to update M, r_cv only
  */
  virtual void update_cell_dependencies(const int cell) = 0;
  
  /**
    calculate \f$ \bar{\bar{\mathbf{R}}}_{\sigma_{t,g}} \f$ , the isotropic components of \f$ \vec{S}_I \f$, and 
    if grey \f$ \mathbf{R}_{\sigma_{s,0}} + \bar{\bar{\nu}}\mathbf{R}_{\sigma_a} \f$ else if MF, calculate
    only \f$ \mathbf{R}_{\sigma_{s,g,0}} \f$ since the scattering source is part of the fixed point iteration    
  */
  virtual void update_group_dependencies(const int grp) = 0;
  
  /// retrieve the already constructed \f$ \bar{\bar{\mathbf{R}}}_{\sigma_{t,g}} \f$
  void get_r_sig_t(Eigen::MatrixXd& r_sig_t);
  
  /** if grey and l_mom = 0, retireve \f$ \mathbf{R}_{\sigma_{s,0}} + \bar{\bar{\nu}}\mathbf{R}_{\sigma_a} \f$, else retrieve
    \f$ \mathbf{R}_{\sigma_{s,g,\text{l_mom}}} \f$
  */
  void get_r_sig_s(Eigen::MatrixXd& r_sig_s, const int l_mom);
  
  /** add in the directional dependent components into \f$ \vec{S}_I} \f$ regardless of grey or MF, the directional components are:
    \f[
      \frac{1}{c\Delta a_{ii} } \mathbf{M} \vec{I}_{n,dir,g} + \frac{1}{c a_{ii} }\mathbf{M} 
        \sum_{j=1}^{\text{stage} - 1}{a_{ij} \vec{k}_{I,j,dir,grp} }
    \f]
  */  
  void get_s_i(Eigen::VectorXd& s_i, const int dir);
    
private:
  const MATRIX_INTEGRATION m_matrix_type; 
  
  const int m_np;
    /// \f$ \mathbf{R}_{\sigma_{a,g}} \f$ evaluated at t_star
  Eigen::MatrixXd m_r_sig_a; 
  
    /// \f$ \mathbf{R}_{\sigma_{s,g}} \f$ evaluated at t_star
  Eigen::MatrixXd m_r_sig_s;
  
  /** evaluate \f$ \mathbf{R}_{\sigma_{s,g}} \f$, copy into this this matrix, evaluate \f$ \mathbf{R}_{\sigma_{a,g}} \f$ and store, copy/add
    to this matrix, then add in \f$ \frac{1}{c\Delta t a_{ii}} \mathbf{M} \f$ component 
  */
  Eigen::MatrixXd m_r_sig_t;
  
  /// \f$ \mathbf{I} \f$
  const Eigen::MatrixXd m_identity_matrix;
  
  /// will actually store the inverse, \f$ \mathbf{R}_{C_v}^{-1} \f$
  Eigen::MatrixXd m_r_cv;
  
  /// \f$ \mathbf{D}_g \f$
  Eigen::MatrixXd m_d_matrix;
  
  /**
    \f[
      \text{m_coefficient} = \left[ \mathbf{I} + \text{m_sn_w} \Delta t a_{ii} \mathbf{R}_{C_v}^{-1} \mathbf{R}_{\sigma_{a,g}} \mathbf{D} \right]
    \f]
  */
  Eigen::MatrixXd m_coefficent;
  
  /// \f$ \mathbf{M} \f$
  Eigen::MatrixXd m_mass;
  
  Eigen::MatrixXd m_no_mu_pos_l_matrix;
  Eigen::MatrixXd m_no_mu_neg_l_matrix;
  Eigen::VectorXd m_no_mu_pos_f_vector;
  Eigen::VectorXd m_no_mu_neg_f_vector;
  
  /** can't be const reference because when call for a material property at DFEM integration points, we modify the scratch vector 
    m_ that is a member of the Materials object
    
    but we don't ever want to be able to change this pointer!
  */
  Materials* const m_materials;
  
  /// pointer won't change and we better not change the Cell_Data that m_cell_data points to
  const Cell_Data* const m_cell_data;  
  
  /// time stepping data that will be used often
  const double m_c;
  double m_dt;
  int m_stage;
  std::vector<double> m_rk_a;
    
  /// for evaluating driving sources
  double m_time;
  
    /// better not change any of the data pointed to by the following pointers (but pointers themselves may change)
  const Temperature_Data* m_t_old;
  const Temperature_Data* m_t_star;  
  const K_Temperature* m_kt;
  
  const Intensity_Data*  m_i_old;
  const K_Intensity* m_ki;
  
  /// cell data that is used repeatedly (set in V_Sweep_Matrix_Creator::update_cell_dependencies() )
  double m_dx;
  int m_cell_num;
    
  /// builder/lumper of matrices
  std::shared_ptr<V_Matrix_Construction> m_mtrx_builder;
  
  /// local vectors
  Eigen::VectorXd m_t_old_vec;
  Eigen::VectorXd m_t_star_vec;
  Eigen::VectorXd m_planck_vec;
  Eigen::VectorXd m_ki_vec;
  Eigen::VectorXd m_kt_vec;
  
  /// isotropic component of xi_i
  Eigen::VectorXd m_xi_isotropic;
  
  /// driving source moments (S_T and S_I)
  Eigen::VectorXd m_driving_source;
  
};

#endif