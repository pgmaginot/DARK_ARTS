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
  
  void get_r_sig_t(Eigen::MatrixXd& r_sig_t);
  void get_r_sig_s(Eigen::MatrixXd& r_sig_s);
  void get_r_sig_s_higher_moment(const int l_mom, Eigen::MatrixXd& r_sig_s_ho);
  void contstruct_sweep_matrices(const int cell, const int grp);  
  
  void construct_l_matrix(const double mu, Eigen::MatrixXd& l_matrix);
  
  void construct_f_vector(const double mu, Eigen::VectorXd& f_vector);
  
  /// data that changes every time we update radiation intensity
  void set_thermal_iteration_data(const Temperature_Data* t_eval, const Temperature_Data* t_old, 
    const K_Temperature* kt );
  

  void set_intensity_iteration_data(const Intensity_Data* i_old, const K_Intensity* ki);
  
  /// data that changes once per stage (per time step)
  
  void set_stage_data(const int stage, const std::vector<double>& rk_a, const double time);

  void set_timestep_data(const double dt);
  
  void get_cell_size(const int cell);
  
private:
  virtual void construct_r_sig_t(void) = 0;
  
  virtual void construct_r_sig_s(void) = 0;
  
  virtual void construct_s_i(void) = 0;

  const MATRIX_INTEGRATION m_matrix_type; 
  
  const int m_np;
  
  Eigen::MatrixXd m_r_sig_t;
  
  const Eigen::MatrixXd m_identity_matrix;
  
  Eigen::MatrixXd m_r_cv;
  Eigen::MatrixXd m_r_sig_a;
  Eigen::MatrixXd m_d_matrix;
  Eigen::MatrixXd m_coefficent;
  
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
  double m_time;
  int m_stage;
  
    /// better not change any of the data pointed to by the following pointers
  const Temperature_Data* m_t_old;
  const Temperature_Data* m_t_star;  
  const K_Temperature* m_kt;
  
  const Intensity_Data*  m_i_old;
  const K_Intensity* m_ki;
  
  std::vector<double> m_rk_a;
  
  double m_dx;
  
  std::shared_ptr<V_Matrix_Construction> m_mtrx_builder;
};

#endif