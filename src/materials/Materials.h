#ifndef Materials_h
#define Materials_h

#include <vector>
#include <memory>
#include <stdlib.h>

#include "Input_Reader.h"
#include "Fem_Quadrature.h"
#include "Temperature_Data.h"
#include "Cell_Data.h"

#include "VAbsorption_Opacity.h"
#include "VScattering_Opacity.h"
#include "VCv.h"
#include "VSource_T.h"
#include "VSource_I.h"
#include "Absorption_Opacity_Rational.h"
#include "Scattering_Opacity_Rational.h"
#include "Absorption_Opacity_Constant.h"
#include "Scattering_Opacity_Constant.h"
#include "Absorption_Opacity_Table.h"
#include "Scattering_Opacity_Table.h"
#include "Source_I_Constant.h"
#include "Source_T_Constant.h"
#include "Cv_Constant.h"
#include "V_XS_Treatment.h"
#include "XS_Treatment_SLXS.h"

class Materials
{

public:
  Materials(const Input_Reader& input_reader, const Fem_Quadrature& fem_quadrature, Cell_Data* cell_ptr);
  ~Materials(){}

  /** The update_sigma_a, update_sigma_s, and update_cv functions interact with the 
        VAbsorption_Opacity, VScattering_Opacity, and VCv objects, respectively.
        
      They then constrcut data of length Fem_Quadrature::m_n_integration_points and store within the Materials object.    

      Materials needs to have building functions that properly translates the Fem_Quadrature::m_n_xs_evaluation_points
      to Fem_Quadrature::m_n_integration_points.  This depends on whether opacity treatment is SLXS, MOMENT_PRESERVING, or INTERPOLATING.
      Solution: create a virtual class with concrete instantiations e.g. 
        V_XS_Treatment, XS_Treatment_SLXS, XS_Treatment_Moment_Preserving, XS_Treatment_Interpolating
  */
  /// calculate \f$ \sigma_a \f$ for all DFEM integration points and groups for cell cell_num
  void update_sigma_a(void);
  /// calculate \f$ \sigma_s \f$ for all DFEM integration points, groups, and scattering moments for cell cell_num
  void update_sigma_s(void);
  /// calculate \f$ C_v \f$ for all DFEM integration points and groups for cell cell_num
  void update_cv(void);
  
  void calculate_local_temp_and_position(const int cell_num, const Temperature_Data& temperature);
  
  /// just return the appropriate data, which the Materials object has stored for the given cell
  // void get_sigma_a(const int cell, const int grp, std::vector<double>& sig_a);
  // void get_sigma_s(const int cell, const int grp, const int l_mom, std::vector<double>& sig_s);
  // void get_cv(const int cell, std::vector<double>& cv);
  
private:
  void load_materials(const Input_Reader& input_reader);

private:
  /** Vectors of pointers to the various opacity, cv, and source objects we will create in
   * the Materials class constructor (based off what we read in the input file)
  */
  std::vector<std::shared_ptr<VAbsorption_Opacity>> m_abs_opacities;
  std::vector<std::shared_ptr<VScattering_Opacity>> m_scat_opacities;
  std::vector<std::shared_ptr<VCv>> m_cv_obj;
  std::vector<std::shared_ptr<VSource_T>> m_source_t;
  std::vector<std::shared_ptr<VSource_I>> m_source_i;
  
  /// holders of evaluated material property data at DFEM integration points
  std::vector<double> m_sig_a;
  std::vector<double> m_sig_s;
  std::vector<double> m_cv;
  
  /** holders of evaluated material property data at cell edges
    All data will be stored as follows:
   \f$ \sigma_{a,L,g=0} \sigma_{a,R,g=0} , \sigma_{a,L,g=1} \dots \f$
   \f$ \sigma_{s,L,g=0,l_mom=0} \sigma_{s,R,g=0,l_mom=0} , 
      \sigma_{s,L,g=0,l_mom=1}  {s,R,g=0,l_mom=1} \dots 
     \f$
   \f$ [ C_{v,L} , C_{v,R} ] \f$
  */
  std::vector<double> m_sig_a_boundary;
  std::vector<double> m_sig_s_boundary;
  std::vector<double> m_cv_boundary;
  
  /// "scratch" vectors to hold a matieral property evaluations at the m_n_xs_evaluation quadrature points
  std::vector<double> m_mat_property_evals;
  
  /// total number of materials
  int m_num_materials = -1;
  
  /// translator from material evaluations at xs quadrature points to xs at DFEM itnegration points
  std::shared_ptr<V_XS_Treatment> m_xs_treatment;
  
  /// XS evaluation points
  int m_n_xs_quad = -1;
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
  
  /// pointer to Cell_Data
  Cell_Data* cell_data_ptr;
  
  /// DFEM integration points
  int m_n_dfem_integration_points = -1;
  /// Used in moment_preserving and interpolating schemes only
  std::vector<double> m_dfem_integration_points;
  
  /// basis functions at selected points
  std::vector<double> m_dfem_at_xs;
  std::vector<double> m_dfem_at_left_edge;
  std::vector<double> m_dfem_at_right_edge;
  
  /// number of DFEM unknwons per cell
  int m_n_el_cell = -1;
  
  /// what material the current cell is in
  int m_current_material = -1;
  
  std::vector<double> m_xs_position;
  
  /// data needed for arranging material property in memory
  int m_n_groups = -1;
  int m_n_l_mom = -1;
  
  std::vector<double> m_mat_mapped;
};

#endif