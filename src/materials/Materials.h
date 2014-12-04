#ifndef Materials_h
#define Materials_h

#include <vector>
#include <memory>
#include <stdlib.h>

#include "Input_Reader.h"
#include "Fem_Quadrature.h"

#include "Cell_Data.h"

#include <Eigen/Dense>

#include "Absorption_Opacity_Rational.h"
#include "Scattering_Opacity_Rational.h"
#include "Absorption_Opacity_Constant.h"
#include "Scattering_Opacity_Constant.h"
#include "Absorption_Opacity_Table.h"
#include "Scattering_Opacity_Table.h"
#include "Absorption_Opacity_Polynomial_Space.h"
#include "Scattering_Opacity_Polynomial_Space.h"
#include "Source_I_Constant.h"
#include "Source_T_Constant.h"
#include "Source_I_MMS.h"
#include "Source_T_MMS.h"
#include "Cv_Constant.h"
#include "Cv_Rational.h"

#include "XS_Treatment_SLXS.h"
#include "XS_Treatment_Moment_Preserving.h"
#include "XS_Treatment_Interpolating.h"

#include "Planck.h"

class Materials
{

public:
  Materials(const Input_Reader& input_reader, 
    const Fem_Quadrature& fem_quadrature, 
    const Cell_Data& cell_data,
    const Angular_Quadrature& angular_quadrature);
    
  virtual ~Materials(){}

  /** The update_sigma_a, update_sigma_s, and update_cv functions interact with the 
        VAbsorption_Opacity, VScattering_Opacity, and VCv objects, respectively.
        
      They then constrcut data of length Fem_Quadrature::m_n_integration_points and store within the Materials object.    

      Materials needs to have building functions that properly translates the Fem_Quadrature::m_n_xs_evaluation_points
      to Fem_Quadrature::m_n_integration_points.  This depends on whether opacity treatment is SLXS, MOMENT_PRESERVING, or INTERPOLATING.
      Solution: create a virtual class with concrete instantiations e.g. 
        V_XS_Treatment, XS_Treatment_SLXS, XS_Treatment_Moment_Preserving, XS_Treatment_Interpolating
  */
  
  void calculate_local_temp_and_position(const int cell_num, const Eigen::VectorXd& temperature);
  
  /// calculate the appropriate data, which the Materials object has stored for the given cell
  void get_sigma_a(const int grp, std::vector<double>& sig_a);
  void get_sigma_s(const int grp, const int l_mom, std::vector<double>& sig_s);
  void get_cv(std::vector<double>& cv);
  
  /// get a vector of unknowns equal to the number of source moment quadrature points
  /// pass this vector to matrix constructors, for consistency!
  void get_temperature_source(const double time, std::vector<double>& t_source_evals);
  void get_intensity_source(const double time, const int grp, const int dir, std::vector<double>& i_source_evals);
  
  void get_sigma_a_boundary(const int grp, std::vector<double>& sig_a);
  void get_sigma_s_boundary(const int grp, const int l_mom, std::vector<double>& sig_s);
  void get_cv_boundary(std::vector<double>& cv);
   
  double get_mf_planck(const double t_eval, const int grp);
  double get_grey_planck(const double t_eval);
   
  void get_mf_planck(const Eigen::VectorXd& t_eval_vec, const int grp, Eigen::VectorXd& planck);
  void get_grey_planck(const Eigen::VectorXd& t_eval_vec, Eigen::VectorXd& planck);
  
  void get_mf_planck_derivative(const Eigen::VectorXd& t_eval_vec, const int grp, Eigen::MatrixXd& d_planck);
  void get_grey_planck_derivative(const Eigen::VectorXd& t_eval_vec, Eigen::MatrixXd& d_planck);
  
  double get_c(void);  

  double get_cell_width(void);
  
private:
  Planck m_planck;

  void load_materials(const Input_Reader& input_reader, const Angular_Quadrature& angular_quadrature);

private:
  /** Vectors of pointers to the various opacity, cv, and source objects we will create in
   * the Materials class constructor (based off what we read in the input file)
  */
  std::vector<std::shared_ptr<VAbsorption_Opacity>> m_abs_opacities;
  std::vector<std::shared_ptr<VScattering_Opacity>> m_scat_opacities;
  std::vector<std::shared_ptr<VCv>> m_cv_obj;
  std::vector<std::shared_ptr<VSource_T>> m_source_t;
  std::vector<std::shared_ptr<VSource_I>> m_source_i;
    
  /// "scratch" vectors to hold a matieral property evaluations at the m_n_xs_evaluation quadrature points
  std::vector<double> m_mat_property_evals;
  
  /// total number of materials
  const int m_num_materials;
  
  /// XS evaluation points
  const int m_n_xs_quad;
  
  /// number of DFEM unknwons per cell
  const int m_n_el_cell;
  
    /// cell width that needs to be accesible to matrix constructors
  double m_dx;
  
  /// pointer to Cell_Data
  const Cell_Data& m_cell_data;
  
  /// number of quadrature points used to evaluate driving source moments
  const int m_n_source_pts;
  
  /// what material the current cell is in
  int m_current_material=-1;
  
  /// translator from material evaluations at xs quadrature points to xs at DFEM integration points
  std::shared_ptr<V_XS_Treatment> m_xs_treatment;
  

  /// resized by fem_quadrature
  std::vector<double> m_xs_eval_quad;
  /// only initialized if used (moment_preserving)
  std::vector<double> m_xs_eval_weights;
  
  /// temperature at xs evaluation points
  std::vector<double> m_t_at_xs_eval_points;
  
  /// temperature at the left and right edges of the cell
  double m_t_left_bound = 0.;
  double m_t_right_bound = 0.;
  
  double m_x_left = 0.;
  double m_x_right = 0.;
  
  /// Used in moment_preserving and interpolating schemes only
  std::vector<double> m_dfem_integration_points;
  
  /// dfem functions at selected points
  std::vector<double> m_dfem_at_xs;
  std::vector<double> m_dfem_at_left_edge;
  std::vector<double> m_dfem_at_right_edge;

  /// physical position at material property evaluation points
  std::vector<double> m_xs_position;
  
  /// energy group bounds
  std::vector<double> m_grp_e_min;
  std::vector<double> m_grp_e_max;
  
  
  /// variables needed to evaluate driving sources at driving source quadrature points
  std::vector<double> m_source_quad;
  std::vector<double> m_position_at_source_quad;    
};

#endif