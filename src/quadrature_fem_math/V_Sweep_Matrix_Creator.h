#ifndef V_Sweep_Matrix_Creator_h
#define V_Sweep_Matrix_Creator_h

#include "Inputs_Allowed.h"
#include "Fem_Quadrature.h"
#include "Eigen/Dense"
#include "Matrix_Construction_Exact.h"
#include "Matrix_Construction_Self_Lumping.h"
#include "Matrix_Construction_Trad_Lumping.h"

#include "Cell_Data.h"
#include "Materials.h"
#include "Temperature_Data.h"
#include "K_Intensity.h"
#include "K_Temperature.h"

#include "Intensity_Data.h"
#include "Intensity_Moment_Data.h"

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
  V_Sweep_Matrix_Creator(const Fem_Quadrature& fem_quadrature, 
    Materials& materials,
    const int n_stages, 
    const double sn_w, 
    const int n_l_mom,
    const Temperature_Data& t_old, 
    const Intensity_Data& i_old,
    const K_Temperature& kt, 
    const K_Intensity& ki,
    const Temperature_Data& t_star);
    
  virtual ~V_Sweep_Matrix_Creator(){}
   
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
  
    /** add in the directional dependent components into \f$ \vec{S}_I} \f$ regardless of grey or MF, the directional components are:
    \f[
      \frac{1}{c\Delta a_{ii} } \mathbf{M} \vec{I}_{n,dir,g} + \frac{1}{c a_{ii} }\mathbf{M} 
        \sum_{j=1}^{\text{stage} - 1}{a_{ij} \vec{k}_{I,j,dir,grp} }
    \f]
    /// this could alternatively be called update_direction_dependencies, with a get_s_i() function added in above...
  */  
  virtual void update_direction_dependencies(const int dir) = 0;
  
  void construct_l_matrix(const double mu, Eigen::MatrixXd& l_matrix);
  
  void construct_f_vector(const double mu, Eigen::VectorXd& f_vector);
    
  /// SDRIK relevant time data (dt doesn't change at different stages, but let's keep it simple)
  void set_time_data( const double dt, const double time_stage, const std::vector<double>& rk_a_of_stage_i, const int stage );
  
  /// retrieve the already constructed \f$ \bar{\bar{\mathbf{R}}}_{\sigma_{t,g}} \f$.  Constructed in \fn update_group_dependencies
  void get_r_sig_t(Eigen::MatrixXd& r_sig_t) const;
  
  /** if grey and l_mom = 0, retireve \f$ \mathbf{R}_{\sigma_{s,0}} + \bar{\bar{\nu}}\mathbf{R}_{\sigma_a} \f$, else retrieve
    \f$ \mathbf{R}_{\sigma_{s,g,\text{l_mom}}} \f$ . Matrix created in \fn update_group_dependencies
  */
  void get_r_sig_s(Eigen::MatrixXd& r_sig_s, const int l_mom) const;
  
  void get_s_i(Eigen::VectorXd& s_i) const;
  
  void get_mass_inverse(Eigen::MatrixXd& m_inv) const;
      
  void set_ard_phi_ptr(Intensity_Moment_Data* ard_phi_ptr);
  
  /// used only by Solution_Saver_K_I
  void calculate_k_i_quantities(void);
  void k_i_get_r_sig_a(Eigen::MatrixXd& r_sig_a) const;
  void k_i_get_r_sig_s_zero(Eigen::MatrixXd& r_sig_s_zero) const;
  void k_i_get_r_sig_t(Eigen::MatrixXd& r_sig_t) const;
  void k_i_get_s_i(Eigen::VectorXd& s_i) const;
  void k_i_get_planck_vec(Eigen::VectorXd& planck) const;
  
private:
  const MATRIX_INTEGRATION m_matrix_type; 
  
protected:
  /** ****************************************************************
    * Variables that are initialzed in the constructor initialization list
    ****************************************************************
   */
   
  /// constants needed to build things
  const double m_sn_w; /// needed here to build source terms in linearization
  const int m_n_l_mom;
  const int m_np; /// needed to know the length/size of vectors/matrices
  
  /// matrices/vectors that hold sweep matrices (end results after calculating linearization, time dependence quantities)
  Eigen::MatrixXd m_r_sig_t;
  std::vector<Eigen::MatrixXd> m_r_sig_s;
  Eigen::VectorXd m_s_i;
  
  /// matrices needed when calculating k_i
  Eigen::MatrixXd m_k_i_r_sig_t;
  Eigen::MatrixXd m_k_i_r_sig_s_zero;
    
  /// quantities commonly used in linearization definitions  
  /// \f$ \mathbf{R}_{\sigma_{a,g}} \f$ evaluated at t_star
  Eigen::MatrixXd m_r_sig_a; 
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
  Eigen::MatrixXd m_coefficient;
  /// \f$ \mathbf{M} \f$ without \f$ \frac{\Delta x}{2} \f$ term
  Eigen::MatrixXd m_mass;
  /// scaled mass matrix
  Eigen::MatrixXd m_dx_div_2_mass;  
  /// inverse of scaled mass matrix
  Eigen::MatrixXd m_mass_inv;
  
   /// isotropic component of xi_i
  Eigen::VectorXd m_xi_isotropic;  
  /// holding source moments (S_T and S_I)
  Eigen::VectorXd m_driving_source;
  
  /// local vectors of respective Data objects
  Eigen::VectorXd m_t_old_vec;
  Eigen::VectorXd m_t_star_vec;
  Eigen::VectorXd m_k_vec;
  Eigen::VectorXd m_planck_vec;
  Eigen::VectorXd m_temp_vec; 
  
  /// gradient matrices/vectors that are not changed during the sweep process
  /// Passed to Matrix_Construction_Objects at construction, then not modified after that
  Eigen::MatrixXd m_no_mu_pos_l_matrix;
  Eigen::MatrixXd m_no_mu_neg_l_matrix;
  Eigen::VectorXd m_no_mu_pos_f_vector;
  Eigen::VectorXd m_no_mu_neg_f_vector;
  
  /** can't be const Materials* because when we call for a material property at DFEM integration points, we modify the scratch vector 
    m_ that is a member of the Materials object...  but we don't ever want to be able to change this pointer!
  */
  Materials& m_materials;   
  
  const double m_c;
  
  /// time stepping data that will be used/changed during the search for each intensity 
  double m_dt;
  int m_stage;
  double m_time;
  
  /// pointers that will not change during a problem, and are used regardless of frequency treatment.  
  const Temperature_Data& m_t_old_ref;
  const Intensity_Data& m_i_old_ref;  
  const K_Temperature& m_kt_ref;  
  const K_Intensity& m_ki_ref;
  
  /// cell data that is used repeatedly (set in V_Sweep_Matrix_Creator::update_cell_dependencies() )
  double m_dx;
  int m_cell_num;
  int m_group_num;
  
  /// could potentially change, though unlikely...
  const Temperature_Data& m_t_star_ref;
  /// used only by MF, but needs to be available to any V_Sweep_Matrix_Creator object since it is called by a 
  /// shared_ptr<V_Sweep_Matrix_Creator>
  Intensity_Moment_Data* m_ard_phi_ptr;
  
  /// vector of length n_stages that will old the rk_a data (changed at each time step stage)
  std::vector<double> m_rk_a;
  
  /** ****************************************************************
    * Variables that are initialzed in the constructor body
    **************************************************************** */
    
  /// builder/lumper of matrices
  std::shared_ptr<V_Matrix_Construction> m_mtrx_builder;

};

#endif