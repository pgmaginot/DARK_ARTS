#ifndef Input_Reader_h
#define Input_Reader_h

#include "tinyxml.h"
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include "Inputs_Allowed.h"
#include "Dark_Arts_Exception.h"

/** @file   Input_Reader.h
  *   @author pmaginot
  *   @brief Declare Input_Reader class
  *   @class Input_Reader
 */
class Input_Reader
{
public:
  Input_Reader(){}
  virtual ~Input_Reader(){}
  
  //! read the supplied input file, start populating data objects
  void read_xml(std::string xmlFile);
  
  /// Output data
  
  
  /// Restart data
  void get_data_dump_path(std::string& dump_path) const { dump_path = m_data_dump_str; return;}
  bool is_mesh_refinement(void) const { return m_is_mesh_refinement; }
  bool is_restart(void) const {return m_is_restart;}
  int get_time_step_restart(void) const {return m_restart_step;}
  void get_initial_input_filename(std::string& base_filename) const { base_filename = m_initial_input_str; return;}
  int get_refinement_factor(void) const {return m_refinement_factor; } 
  
  /// Functions called by Fem_Quadrature class to access input data
  int get_dfem_degree(void) const {  return m_dfem_trial_space_degree; }
  QUADRATURE_TYPE get_dfem_interpolation_point_type(void) const { return m_dfem_interpolation_point_type; } 
  int get_opacity_degree(void) const { return m_opacity_polynomial_degree; }
  OPACITY_TREATMENT get_opacity_treatment(void) const { return m_opacity_treatment; }
  QUADRATURE_TYPE get_opacity_interpolation_point_type(void) const { return m_opacity_interpolation_point_type; }
  MATRIX_INTEGRATION get_integration_method(void) const { return m_integration_type; }
  
  /// Functions called by Cell_Data class to access input data
  int get_n_regions(void) const { return m_number_regions; }
  void get_cells_per_region_vector(std::vector<int>& cells_per_region) const {cells_per_region = m_cells_per_region; return;}
  double get_region_left_bound(int reg_num) const;
  double get_region_right_bound(int reg_num) const;
  GRID_SPACING get_region_spacing(int reg_num) const;
  int get_region_material_number(int reg_num) const;
  double get_min_cell_size(int reg_num) const;
  double get_r_factor(int reg_num) const;
  
  /// called by Time_Data class
  double get_t_start(void) const;
  double get_t_end(void) const;
  double get_dt_min(void) const;
  double get_dt_max(void) const;
  TIME_SOLVER get_time_solver(void) const;
  STARTING_METHOD get_starting_time_method(void) const;
  /// if STARTING_METHOD==EXPONENTIAL, need to know the ratio between time step sizes
  double get_time_start_exponential_ratio(void) const;
  /// if STARTING_METHOD == VECTOR, need vector of time step sizes and vector of the number of steps for each time step size
  void get_time_start_vectors(std::vector<double>& step_size_in_vector_stage, std::vector<int>& steps_in_vector_stage) const;
  int get_number_of_vector_stages(void) const;
  /// if STARTING_METHOD == RAMP, need number of time steps to do before hitting full time step
  int get_number_of_ramp_steps(void) const;
  
  /// Additional functions needed by Angular_Quadrature
  int get_number_of_groups(void) const;
  int get_number_of_angles(void) const;
  ANGULAR_QUADRATURE_TYPE get_angular_quadrature_type(void) const;
  int get_number_of_legendre_moments(void) const;
  void get_lower_energy_bounds(std::vector<double>& low_bounds) const;
  void get_upper_energy_bounds(std::vector<double>& upper_bounds) const;
  
  /// Functions for the Materials class and related objects
  int get_number_of_materials(void) const;
  OPACITY_TYPE get_abs_opacity_type(const int mat_num) const;
  OPACITY_TYPE get_scat_opacity_type(const int mat_num) const;
  CV_TYPE get_cv_type(const int mat_num) const;
  FIXED_SOURCE_TYPE get_temperature_source_type(const int mat_num) const {return m_material_temperature_source_type[mat_num];  }
  FIXED_SOURCE_TYPE get_radiation_source_type(const int mat_num) const { return m_material_radiation_source_type[mat_num];   }
  double get_abs_double_constant_1(const int mat_num) const;
  double get_abs_double_constant_2(const int mat_num) const;
  double get_scat_double_constant_1(const int mat_num) const;
  double get_scat_double_constant_2(const int mat_num) const;
  int get_abs_int_constant(const int mat_num) const;
  int get_scat_int_constant(const int mat_num) const;
  void get_abs_file_str(const int mat_num, std::string& filename) const;
  void get_scat_file_str(const int mat_num, std::string& filename) const;  
  double get_cv_constant(const int mat_num) const;  
  bool use_weird_units(void) const;
  UNITS_TYPE get_units_type(void) const;
  void get_scattering_poly_coeff(const int mat_num , std::vector<double>& coeff) const;
  void get_absorption_poly_coeff(const int mat_num , std::vector<double>& coeff) const;
  int  get_rational_cv_power(const int mat_num) const { return m_cv_rational_powers[mat_num] ; }
  double get_rational_cv_offset(const int mat_num) const { return m_cv_rational_offsets[mat_num] ;}
  /// MMS access functions
  RADIATION_ANGLE_MMS get_mms_radiation_angle_dependence(void) const {return m_mms_rad_angle;}
  RADIATION_SPACE_MMS get_mms_radiation_space_dependence(void) const {return m_mms_rad_space;}
  TEMPERATURE_SPACE_MMS get_mms_temperature_space_dependence(void) const {return m_mms_temp_space;}
  TIME_MMS_TYPE get_mms_time_dependence(void) const {return m_mms_time;}
  void get_mms_time_coeff(std::vector<double>& time_coeff) const {time_coeff = m_time_mms_const ; return;}
  void get_mms_radiation_space_coeff(std::vector<double>& space_coeff) const { space_coeff = m_rad_space_mms_const; return;}
  void get_mms_temperature_space_coeff(std::vector<double>& space_coeff) const { space_coeff = m_temp_space_mms_const; return;}
  void get_mms_angle_coeff(std::vector<double>& angle_coeff) const { angle_coeff = m_mms_angle_coeff; return;}
  
  /// Solver tolerances 
  WG_SOLVE_TYPE get_within_group_solve_type(void) const;
  double get_within_group_solve_tolerance(void) const;
  double get_between_group_solve_tolerance(void) const;
  int get_max_number_sweeps(void) const;
  int get_max_ard_iterations(void) const;
  ARD_SOLVE_TYPE get_ard_solve_type(void) const;
  double get_thermal_tolerance(void) const;
  
  /// IC_BC block
  /**
    Currently assuming one initial condition type for the whole problem.  Could change to region specific
  */
  TEMPERATURE_IC_TYPE get_temperature_ic_type(void) const;
  double get_region_temperature(const int reg_num) const;
  double get_region_radiation_temperature(const int reg_num) const;
  RADIATION_IC_TYPE get_radiation_ic_type(void) const;
  RADIATION_BC_TYPE get_radiation_bc_type_left(void) const;
  RADIATION_BC_TYPE get_radiation_bc_type_right(void) const;
  
  BC_ANGLE_DEPENDENCE get_left_bc_angle_dependence(void) const;
  BC_ANGLE_DEPENDENCE get_right_bc_angle_dependence(void) const;
  
  double get_left_bc_start_time(void) const;
  double get_left_bc_end_time(void) const;
  double get_right_bc_start_time(void) const;
  double get_right_bc_end_time(void) const;

  double get_left_bc_constant(void) const;
  double get_right_bc_constant(void) const;
  
  INCIDENT_BC_VALUE_TYPE get_left_bc_value_type(void) const;
  INCIDENT_BC_VALUE_TYPE get_right_bc_value_type(void) const;

  BC_ENERGY_DEPENDENCE get_left_bc_energy_dependence(void) const;
  BC_ENERGY_DEPENDENCE get_right_bc_energy_dependence(void) const;
  
  
protected:
  /** variables that will be used to store data from input file
    this data will then be used by other class initializers **/
  
  /// Output block
  int m_time_steps_per_dump = -1;
  OUTPUT_TYPE m_output_type = INVALID_OUTPUT_TYPE;
  bool m_is_restart = false;
  int m_restart_step = -1;
  
  /// Restart type
  RESTART_TYPE m_restart_type = INVALID_RESTART_TYPE;
  int m_refinement_factor = -1;
  std::string m_initial_input_str;
  bool m_is_mesh_refinement = false;
  std::string m_data_dump_str;
  
  /// regions input block
  int m_number_regions = -1;
  int m_number_materials = -1;
  std::vector<int> m_cells_per_region;
  std::vector<int> m_region_material_numbers;
  std::vector<GRID_SPACING> m_region_spacing_type;
  std::vector<double> m_region_spacing_constant;
  std::vector<double> m_region_min_size;
  std::vector<double> m_region_left_bounds;
  std::vector<double> m_region_right_bounds;
  
  /// materials input block
  std::vector<OPACITY_TYPE> m_material_scattering_opacity_type;
  std::vector<OPACITY_TYPE> m_material_absorption_opacity_type;
  std::vector<CV_TYPE> m_material_cv_type;
  std::vector<FIXED_SOURCE_TYPE> m_material_radiation_source_type;
  std::vector<FIXED_SOURCE_TYPE> m_material_temperature_source_type;
  std::vector<std::string> m_scat_opacity_str;
  std::vector<std::string> m_abs_opacity_str;
  std::vector<int> m_abs_opacity_integer_constants;
  std::vector<double> m_abs_opacity_double_constants_1;
  std::vector<double> m_abs_opacity_double_constants_2;
  std::vector<int> m_scat_opacity_integer_constants;
  std::vector<double> m_scat_opacity_double_constants_1;
  std::vector<double> m_scat_opacity_double_constants_2;
  std::vector<double> m_cv_constants;
  bool m_weird_units = false;
  UNITS_TYPE m_units_type = INVALID_UNITS_TYPE;
  std::vector< std::vector<double> > m_scat_opacity_poly;
  std::vector< std::vector<double> > m_abs_opacity_poly;
  std::vector<int> m_cv_rational_powers;
  std::vector<double> m_cv_rational_offsets;
  
  /// MMS constants
  RADIATION_SPACE_MMS m_mms_rad_space = INVALID_RADIATION_SPACE_MMS;
  TEMPERATURE_SPACE_MMS m_mms_temp_space = INVALID_TEMPERATURE_SPACE_MMS;  
  TIME_MMS_TYPE m_mms_time = INVALID_TIME_MMS_TYPE;
  RADIATION_ANGLE_MMS m_mms_rad_angle = INVALID_RADIATION_ANGLE_MMS;
  std::vector<double> m_rad_space_mms_const;
  std::vector<double> m_temp_space_mms_const;
  std::vector<double> m_time_mms_const;
  std::vector<double> m_mms_angle_coeff;
  
  
  /// spatial discretization input block
  int m_dfem_trial_space_degree = -1;
  MATRIX_INTEGRATION m_integration_type = INVALID_MATRIX_INTEGRATION;
  QUADRATURE_TYPE m_dfem_interpolation_point_type = INVALID_QUADRATURE_TYPE;
  QUADRATURE_TYPE m_opacity_interpolation_point_type = INVALID_QUADRATURE_TYPE;
  OPACITY_TREATMENT m_opacity_treatment = INVALID_OPACITY_TREATMENT;
  int m_opacity_polynomial_degree = -1;
  
  //! angular discretization data
  ANGULAR_QUADRATURE_TYPE m_angular_quadrature_type = INVALID_ANGULAR_QUADRATURE_TYPE;
  int m_number_angles = -1;
  int m_number_groups = -1;
  int m_n_leg_moments = -1;
  std::vector<double> m_group_lower_bounds;
  std::vector<double> m_group_upper_bounds;
  
  //! time discretization data
  double m_dt_min = -1.;
  double m_dt_max = -2.;
  double m_t_start = -3.;
  double m_t_end = -4.;
  TIME_SOLVER m_time_step_scheme = INVALID_TIME_SOLVER;
  STARTING_METHOD m_time_starting_method = INVALID_STARTING_METHOD;
  std::vector<double> m_vector_start_sizes;
  std::vector<int> m_vector_start_step_numbers;
  double m_exponential_ratio = -1.;
  int m_ramp_steps = -1;
  int m_num_vec_stages = -1;
  
  /// radiation solver type data
  /// within group solve type
  WG_SOLVE_TYPE m_wg_solve_type = INVALID_WG_SOLVE_TYPE;
  /// within group scattering tolerance
  double m_wg_tolerance = 3.;
  /// between group scattering/absorption/re-emission tolerance
  double m_bg_tolerance = 2.;
  /// maximum number of sweeps per within group radiation solve
  int m_max_num_sweeps = -1;
  /** Things only necessary for multi-frequency  */
  ARD_SOLVE_TYPE m_ard_solve_type = INVALID_ARD_SOLVE_TYPE;
  int m_max_ard_iterations = -1;
  double m_thermal_tolerance = 1.;
  
  /// BC_IC block variables
  TEMPERATURE_IC_TYPE m_temperature_ic_type = INVALID_TEMPERATURE_IC_TYPE;
  std::vector<double> m_region_temperature;
  std::vector<double> m_region_radiation_temperature;
  RADIATION_IC_TYPE m_radiation_ic_type  = INVALID_RADIATION_IC_TYPE;

  RADIATION_BC_TYPE m_rad_bc_right = INVALID_RADIATION_BC_TYPE;
  double m_right_bc_value = -1.;
  INCIDENT_BC_VALUE_TYPE m_rad_bc_right_value_type = INVALID_INCIDENT_BC_VALUE_TYPE;
  BC_ANGLE_DEPENDENCE m_right_bc_angle_dependence  = INVALID_BC_ANGLE_DEPENDENCE;
  BC_ENERGY_DEPENDENCE m_right_bc_energy_dependence = INVALID_BC_ENERGY_DEPENDENCE;
  BC_TIME_DEPENDENCE m_right_bc_time_dependence = INVALID_BC_TIME_DEPENDENCE;
  double m_bc_right_start_time = -1.;
  double m_bc_right_end_time = -2.;
  
  RADIATION_BC_TYPE m_rad_bc_left  = INVALID_RADIATION_BC_TYPE;
  double m_left_bc_value = -1.;
  INCIDENT_BC_VALUE_TYPE m_rad_bc_left_value_type = INVALID_INCIDENT_BC_VALUE_TYPE;
  BC_ANGLE_DEPENDENCE m_left_bc_angle_dependence = INVALID_BC_ANGLE_DEPENDENCE;
  BC_ENERGY_DEPENDENCE m_left_bc_energy_dependence = INVALID_BC_ENERGY_DEPENDENCE;
  BC_TIME_DEPENDENCE m_left_bc_time_dependence = INVALID_BC_TIME_DEPENDENCE;
  double m_bc_left_start_time = -1.;
  double m_bc_left_end_time = -2.;

  void load_restart_problem(TiXmlElement* restart_elem);
  void load_from_scratch_problem(TiXmlDocument& doc);
  
  //! readers for the specific xml blocks
  int load_region_data(TiXmlElement* region_element);
  int load_material_data(TiXmlElement* mat_elem);
  int load_time_stepping_data(TiXmlElement* time_elem);
  int load_spatial_discretization_data(TiXmlElement* spatial_element);
  int load_angular_discretization_data(TiXmlElement* angle_element);
  int load_solver_data(TiXmlElement* solver_element);
  int load_bc_ic_data(TiXmlElement* bc_ic_element);
  
  /**loaders for MMS types
    I(x,mu,t) = (mu_fun)*(x_fun_rad)*(t_fun_rad)
    T(x,t) = (x_fun_temp)*(t_fun_temp)
  */
  void load_mms_poly_constants(TiXmlElement* mms_element, std::vector<double>& poly_constants);
  void load_mms_cos_constants(TiXmlElement* mms_element, std::vector<double>& cos_constants);  

};


#endif