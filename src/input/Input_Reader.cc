/** @file   Input_Reader.cc
  *   @author pmaginot
  *   @brief Implement the Input_Reader class that reads XML files with TinyXml
 */
         
#include "Input_Reader.h"
// ##########################################################
// Public functions 
// ##########################################################

bool Input_Reader::read_xml(std::string xmlFile)
{
  bool good_exit = true;
  
  TiXmlDocument doc( xmlFile.c_str() );
  bool load_success = doc.LoadFile();
  
  if(load_success){
    std::cout << "Found input file: " << xmlFile << std::endl << std::endl;
  }
  else{
    std::cerr << "Could not open file: " << xmlFile << std::endl;
    exit(EXIT_FAILURE);
  }
  
  TiXmlElement* inp_block = doc.FirstChildElement( "INPUT_FILE" );
  
  if( !inp_block )
  {
    std::cerr << "INPUT_FILE block not found.  Error. \n" ;
    exit(EXIT_FAILURE);
  }
  
  //! Every Input File Is Required To Have These Blocks
  TiXmlElement* reg_elem = inp_block->FirstChildElement( "REGIONS" );
  TiXmlElement* mat_elem = inp_block->FirstChildElement( "MATERIALS" );
  TiXmlElement* time_elem = inp_block->FirstChildElement( "TIME" );
  TiXmlElement* discr_elem = inp_block->FirstChildElement( "SPATIAL_DISCRETIZATION" );
  TiXmlElement* angle_elem = inp_block->FirstChildElement( "ANGULAR_DISCRETIZATION" );
  TiXmlElement* solver_elem = inp_block->FirstChildElement( "SOLVER" );
  TiXmlElement* bc_ic_elem = inp_block->FirstChildElement( "BC_IC" );
  
  if(!reg_elem)
  {
    std::cerr << "REGIONS block not found.\n";
    exit(EXIT_FAILURE);
  }
  
  if(!mat_elem)
  {
    std::cerr << "MATERIALS block not found.  Error \n";
    exit(EXIT_FAILURE);
  } 
  
  if(!time_elem)
  {
    std::cerr << "TIME block not found.  Error \n";
    exit(EXIT_FAILURE);
  }
  
  if(!discr_elem)
  {
    std::cerr << "SPATIAL_DISCRETIZATION block not found.  Error \n";
    exit(EXIT_FAILURE);
  }
  
  if(!angle_elem)
  {
    std::cerr << "ANGULAR_DISCRETIZATION block not found.  Error \n";
    exit(EXIT_FAILURE);
  }
  
  if(!solver_elem)
  {
    std::cerr << "SOLVER block not found.  Error \n";
    exit(EXIT_FAILURE);
  }
  
  if(!bc_ic_elem)
  {
    std::cerr << "BC_IC block not found.  Error \n";
    exit(EXIT_FAILURE);
  }
  
  int region_return = -1;
  int material_return = -1;
  int time_return = -1;
  int spatial_return = -1;
  int angle_return = -1;
  int solver_return = -1;
  int bc_ic_return = -1;
  
  //! Load Data Appropriately from each input block
  region_return = load_region_data(reg_elem);
  std::cout << "Region block read\n\n" ;
  
  material_return = load_material_data(mat_elem);
  std::cout << "Material block read\n\n" ;
  
  time_return = load_time_stepping_data(time_elem);
  std::cout << "Time block read\n\n" ;
  
  spatial_return = load_spatial_discretization_data(discr_elem);
  std::cout << "Spatial discretization block read\n\n" ;
  
  angle_return = load_angular_discretization_data(angle_elem);
  std::cout << "Angle block read\n\n" ;
  
  solver_return = load_solver_data(solver_elem);
  std::cout << "Solver block read \n\n";
  
  bc_ic_return = load_bc_ic_data(bc_ic_elem);
  std::cout << "BC_IC block read \n\n";
  
  if( (region_return < 0) || (material_return < 0) || (time_return < 0) || 
      (spatial_return < 0) || (angle_return < 0) || (solver_return < 0) ||
      (bc_ic_return < 0) )
  {
    good_exit = false;
  }  
  
  return good_exit;
}

// ##########################################################
// Const (get) Public functions 
// ##########################################################

int Input_Reader::get_dfem_degree(void) const
{
  return m_dfem_trial_space_degree;
}

QUADRATURE_TYPE Input_Reader::get_dfem_interpolation_point_type(void) const
{
  return m_dfem_interpolation_point_type;
}

int Input_Reader::get_opacity_degree(void) const 
{
  return m_opacity_polynomial_degree;
}

OPACITY_TREATMENT Input_Reader::get_opacity_treatment(void) const
{
  return m_opacity_treatment;
}

QUADRATURE_TYPE Input_Reader::get_opacity_interpolation_point_type(void) const
{
  return m_opacity_interpolation_point_type;
}

MATRIX_INTEGRATION Input_Reader::get_integration_method(void) const
{
  return m_integration_type;
}

int Input_Reader::get_n_regions(void) const
{
  return m_number_regions;
}

void Input_Reader::get_cells_per_region_vector(std::vector<int>& cell_per_region) const
{
  cell_per_region.resize(m_number_regions,0);
  for(int i=0;i<m_number_regions;i++)
    cell_per_region[i] = m_cells_per_region[i];
    
  return;
}

GRID_SPACING Input_Reader::get_region_spacing(int reg_num) const
{   
  return m_region_spacing_type[reg_num];
}

int Input_Reader::get_region_material_number(int reg_num) const
{    
  return m_region_material_numbers[reg_num];
}

double Input_Reader::get_region_left_bound(int reg_num) const
{
  return m_region_left_bounds[reg_num];
}

double Input_Reader::get_region_right_bound(int reg_num) const
{
  return m_region_right_bounds[reg_num];
}

double Input_Reader::get_min_cell_size(int reg_num) const
{
  return m_region_min_size[reg_num];
}

double Input_Reader::get_r_factor(int reg_num) const
{
  return m_region_spacing_constant[reg_num];
}

double Input_Reader::get_t_start(void) const
{
  return m_t_start;
}

double Input_Reader::get_t_end(void) const
{
  return m_t_end;
}

double Input_Reader::get_dt_min(void) const
{
  return m_dt_min;
}

double Input_Reader::get_dt_max(void) const
{
  return m_dt_max;
}

TIME_SOLVER Input_Reader::get_time_solver(void) const
{
  return m_time_step_scheme;
}

STARTING_METHOD Input_Reader::get_starting_time_method(void) const
{
  return m_time_starting_method;
}

double Input_Reader::get_time_start_exponential_ratio(void) const
{
  return m_exponential_ratio;
}
  
void Input_Reader::get_time_start_vectors(std::vector<double>& step_size_in_vector_stage, std::vector<int>& steps_in_vector_stage) const
{
  step_size_in_vector_stage = m_vector_start_sizes;
  steps_in_vector_stage = m_vector_start_step_numbers;
  return;
}

int Input_Reader::get_number_of_ramp_steps(void) const
{
  return m_ramp_steps;
}

int Input_Reader::get_number_of_groups(void) const
{
  return m_number_groups;
}

int Input_Reader::get_number_of_angles(void) const
{
  return m_number_angles;
}

ANGULAR_QUADRATURE_TYPE Input_Reader::get_angular_quadrature_type(void) const
{
  return m_angular_quadrature_type;
}

int Input_Reader::get_number_of_legendre_moments(void) const
{
  return m_n_leg_moments;
}

void Input_Reader::get_lower_energy_bounds(std::vector<double>& low_bounds) const
{
  low_bounds = m_group_lower_bounds;
  return;
}

void Input_Reader::get_upper_energy_bounds(std::vector<double>& upper_bounds) const
{
  upper_bounds = m_group_upper_bounds;
  return;
}

int Input_Reader::get_number_of_materials(void) const
{
  return m_number_regions;
}

OPACITY_TYPE Input_Reader::get_abs_opacity_type(const int mat_num) const
{
  return m_material_absorption_opacity_type[mat_num];
}

OPACITY_TYPE Input_Reader::get_scat_opacity_type(const int mat_num) const
{
  return m_material_scattering_opacity_type[mat_num];
}

CV_TYPE Input_Reader::get_cv_type(const int mat_num) const
{
  return m_material_cv_type[mat_num];
}

FIXED_SOURCE_TYPE Input_Reader::get_temperature_source_type(const int mat_num) const
{
  return m_material_temperature_source_type[mat_num];
}

FIXED_SOURCE_TYPE Input_Reader::get_radiation_source_type(const int mat_num) const
{
  return m_material_radiation_source_type[mat_num]; 
}  

double Input_Reader::get_abs_double_constant_1(const int mat_num) const
{
  return m_abs_opacity_double_constants_1[mat_num];
}

double Input_Reader::get_abs_double_constant_2(const int mat_num) const
{
  return m_abs_opacity_double_constants_2[mat_num];
}

double Input_Reader::get_scat_double_constant_1(const int mat_num) const
{
  return m_scat_opacity_double_constants_1[mat_num];
}

double Input_Reader::get_scat_double_constant_2(const int mat_num) const
{
  return m_scat_opacity_double_constants_2[mat_num];;
}

int Input_Reader::get_abs_int_constant(const int mat_num) const
{
  return m_abs_opacity_integer_constants[mat_num];
}

int Input_Reader::get_scat_int_constant(const int mat_num) const
{
  return m_scat_opacity_integer_constants[mat_num];
}

void Input_Reader::get_abs_file_str(const int mat_num, std::string& filename) const
{
  filename = m_abs_opacity_str[mat_num];
  return;
}

void Input_Reader::get_scat_file_str(const int mat_num, std::string& filename) const
{
  filename = m_scat_opacity_str[mat_num];
  return;
}

double Input_Reader::get_cv_constant(const int mat_num) const
{
  return m_cv_constants[mat_num];
}

WG_SOLVE_TYPE Input_Reader::get_within_group_solve_type(void) const
{
  return m_wg_solve_type;
}

double Input_Reader::get_within_group_solve_tolerance(void) const
{
  return m_wg_tolerance;
}

double Input_Reader::get_between_group_solve_tolerance(void) const
{
  return m_bg_tolerance;
}

int Input_Reader::get_max_number_sweeps(void) const
{
  return m_max_num_sweeps;
}

int Input_Reader::get_max_ard_iterations(void) const
{
  return m_max_ard_iterations;
}
  
ARD_SOLVE_TYPE Input_Reader::get_ard_solve_type(void) const
{
  return m_ard_solve_type;
}

bool Input_Reader::use_weird_units(void) const
{
  return m_weird_units;
}

UNITS_TYPE Input_Reader::get_units_type(void) const
{
  return m_units_type;
}

double Input_Reader::get_thermal_tolerance(void) const
{
  return m_thermal_tolerance;
}

TEMPERATURE_IC_TYPE Input_Reader::get_temperature_ic_type(void) const
{
  return m_temperature_ic_type;
}

double Input_Reader::get_region_temperature(const int reg_num) const
{
  if(reg_num > m_number_regions)
  {
    std::cerr << "Asking for a region temperature in a region greater than n_regions" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  return m_region_temperature[reg_num];
}

BC_ANGLE_DEPENDENCE Input_Reader::get_left_bc_angle_dependence(void) const
{
  return m_left_bc_angle_dependence;
}

BC_ANGLE_DEPENDENCE Input_Reader::get_right_bc_angle_dependence(void) const
{
  return m_right_bc_angle_dependence;
}

double Input_Reader::get_bc_start_time(void) const
{
  return m_bc_start_time;
}

double Input_Reader::get_bc_end_time(void) const
{
  return m_bc_end_time;
}

double Input_Reader::get_region_radiation_temperature(const int reg_num) const
{
  std::cerr << "Asking for a region radiation temperature in a region greater than n_regions" << std::endl;
  exit(EXIT_FAILURE);
  
  return m_region_radiation_temperature[reg_num];
}

RADIATION_IC_TYPE Input_Reader::get_radiation_ic_type(void) const
{
  return m_radiation_ic_type;
}

RADIATION_BC_TYPE Input_Reader::get_radiation_bc_type_left(void) const
{
  return m_rad_bc_left;
}

RADIATION_BC_TYPE Input_Reader::get_radiation_bc_type_right(void) const
{
  return m_rad_bc_right;
}

  double Input_Reader::get_left_bc_constant(void) const
  {
    return m_left_bc_value;
  }
  
  double Input_Reader::get_right_bc_constant(void) const
  {
    return m_right_bc_value;
  } 
  

/* ***************************************************
 *
 *  Protected Functions
 *
 * ************************************************** */

int Input_Reader::load_region_data(TiXmlElement* region_element)
{  
  /// Get the number of regions in the problem
  TiXmlElement* num_regions_elem = region_element->FirstChildElement( "Number_of_regions" );
  if(num_regions_elem)
    m_number_regions = atoi(num_regions_elem->GetText());
  else
  {
    std::cerr << "Error.  Missing Number of regions element";
    exit(EXIT_FAILURE);
  }
  
  TiXmlElement* num_materials_elem = region_element->FirstChildElement( "Number_of_materials");
  if(num_regions_elem)
    m_number_materials = atoi(num_materials_elem->GetText());
  else
  {
    std::cerr << "Error.  Missing Number of materials element";
    exit(EXIT_FAILURE);
  }
  
  m_cells_per_region.resize(m_number_regions,0);
  m_region_material_numbers.resize(m_number_regions,0);
  m_region_spacing_type.resize(m_number_regions,INVALID_GRID_SPACING);
  m_region_left_bounds.resize(m_number_regions,-1.0);
  m_region_right_bounds.resize(m_number_regions,-2.0);
  m_region_spacing_constant.resize(m_number_regions,0.0);
  m_region_min_size.resize(m_number_regions,0.0);
  
  TiXmlElement* region_id = region_element->FirstChildElement( "Region");
  if(!region_id)
  {
    std::cerr << "Missing Individual Region Elements" << std::endl;
    exit(EXIT_FAILURE);
  } 
  /// loop over the specified number of regions
  for(int i=0; i<m_number_regions ; i++)
  {
    /// Get all of the elements required for a valid region definition
    int reg_num = atoi(region_id->GetText());
    if(reg_num != i)
    {
      std::cerr << "Error.  Expected Region: " << i << " Got Region: " << reg_num << std::endl;
      exit(EXIT_FAILURE);
    }
    TiXmlElement* n_cells = region_id->FirstChildElement( "N_cells" );
    TiXmlElement* x_left = region_id->FirstChildElement("Left_bound");
    TiXmlElement* x_right = region_id->FirstChildElement("Right_bound");
    TiXmlElement* spacing = region_id->FirstChildElement("Spacing");
    TiXmlElement* mat_number = region_id->FirstChildElement("Material_number");
    
    if(!n_cells || !x_left || !x_right || !spacing || !mat_number)
    {
      std::cerr << "Error.  Region: " << reg_num << " is missing required elements" << std::endl;
      exit(EXIT_FAILURE);
    }
    
    /// store all region specific data    
    /// check region number, demand that all regions be input in order from 0 ... n-1
    m_cells_per_region[i] = atoi(n_cells->GetText());
    if(m_cells_per_region[i] < 0)
    {
      std::cerr << "Error.  Region " << i << " Number of Cells must be a positive integer" << std::endl;
      exit(EXIT_FAILURE);
    }
    
    /// check material number, must be in range: 0 ... n_mat - 1
    m_region_material_numbers[i] = atoi(mat_number->GetText());
    if( (m_region_material_numbers[i] < 0) ||(m_region_material_numbers[i] > (m_number_materials-1)))
    {
      std::cerr << "Error.  Region " << i << " Material number out of range" << std::endl;
      exit(EXIT_FAILURE);
    }
    
    /// get grid spacing type and check against the supported values
    std::string spacing_str = spacing->GetText();
    // change to all upper case 
    transform(spacing_str.begin() , spacing_str.end() , spacing_str.begin() , toupper);
    if(spacing_str == "EQUAL")
      m_region_spacing_type[i] = EQUAL;
    else if(spacing_str == "LOG")
      m_region_spacing_type[i] = LOG;

    if(m_region_spacing_type[i] == INVALID_GRID_SPACING)
    {
      std::cerr << "Error.  Region " << i << " Invalid Grid Spacing" << std::endl;
      exit(EXIT_FAILURE);
    }
    
    /// load in scaling factor for logarithmic grid spacing and check that it is greter than 1
    if( (m_region_spacing_type[i] == LOG) )
    {
      TiXmlElement* space_factor = spacing->FirstChildElement("Log_space_factor");
      TiXmlElement* min_dx = spacing->FirstChildElement("Min_cell_size");
      if(!space_factor)
      {
        std::cerr << "Error.  Region " << i << " Missing size factor for logarithmic grid spacing" << std::endl;
        exit(EXIT_FAILURE);        
      }
      if(!min_dx)
      {
        std::cerr << "Error.  Region " << i << " Missing sminimum cell size for logarithmic spacing" << std::endl;
        exit(EXIT_FAILURE);        
      }
      
      m_region_spacing_constant[i] = atof( space_factor->GetText() );
      if(m_region_spacing_constant[i] < 0.)
      {
        std::cerr << "Error.  Region " << i << " Log spacing factors must be > 0" << std::endl;
        exit(EXIT_FAILURE);
      }
      m_region_min_size[i] = atof( min_dx->GetText() );
      if(m_region_min_size[i] < 0.)
      {
        std::cerr << "Error.  Region " << i << " Minimum cell spacing must be > 0" << std::endl;
        exit(EXIT_FAILURE);
      }
      
    }
    
    /// get left and right values of region, check that x_L < x_R
    m_region_left_bounds[i] = atof( x_left->GetText() ) ;
    m_region_right_bounds[i] = atof(x_right->GetText() ) ;
    if(m_region_left_bounds[i] > m_region_right_bounds[i] )
    {
      std::cerr << "Error.  Region " << i << " x_L > x_R" << std::endl;
      exit(EXIT_FAILURE);
    }
      
    /// go to the next region
    region_id = region_id->NextSiblingElement( "Region");
    if( (region_id == 0) && (i< m_number_regions - 1))
    {
      std::cerr << "Insufficient number of Region elements for the number of regions" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
    
  return 0;
}

int Input_Reader::load_material_data(TiXmlElement* mat_elem)
{
  m_material_absorption_opacity_type.resize(m_number_materials,INVALID_OPACITY_TYPE);  
  m_material_scattering_opacity_type.resize(m_number_materials,INVALID_OPACITY_TYPE);  
  m_material_cv_type.resize(m_number_materials,INVALID_CV_TYPE);
  m_material_radiation_source_type.resize(m_number_materials,INVALID_FIXED_SOURCE_TYPE);
  m_material_temperature_source_type.resize(m_number_materials,INVALID_FIXED_SOURCE_TYPE);
  
  m_scat_opacity_str.resize(m_number_materials);
  m_abs_opacity_str.resize(m_number_materials);
  m_abs_opacity_integer_constants.resize(m_number_materials);
  m_abs_opacity_double_constants_1.resize(m_number_materials);
  m_abs_opacity_double_constants_2.resize(m_number_materials);
  m_scat_opacity_integer_constants.resize(m_number_materials);
  m_scat_opacity_double_constants_1.resize(m_number_materials);
  m_scat_opacity_double_constants_2.resize(m_number_materials);
  m_cv_constants.resize(m_number_materials);
  
  TiXmlElement* units_elem = mat_elem->FirstChildElement("Units");
  if(!units_elem)
  {
    std::cerr << "Missing required Units block in MATERIALS.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    std::string units_str = units_elem->GetText();    
    transform(units_str.begin() , units_str.end() , units_str.begin() , toupper);
    
    if( units_str == "CM_SH_KEV")
    {
      m_units_type = CM_SH_KEV;
      m_weird_units = false;
      std::cout << "Not using any funny units. \n";
      std::cout << "Lengths in cm \n";
      std::cout << "Time in shakes \n";
      std::cout << "Temperatures / Energies in keV\n"; 
      
    }
    else if( units_str == "UNITY" )
    {
      m_units_type = UNITY;
      m_weird_units = true;
      std::cout << "Using unity / dimensionless numbers\n";
      std::cout << "a=c=1" << std::endl;
    }
    else
    {
      std::cerr << "Invalid units type.\n";
      exit(EXIT_FAILURE);
    }
  }
  
  TiXmlElement* mat_descr = mat_elem->FirstChildElement("Material");
  for(int mat_cnt = 0; mat_cnt < m_number_materials ; mat_cnt++)
  {    
    if(!mat_descr)
    {
      std::cerr << "Error.  Missing Material element in MATERIALS block.  Expected: " << m_number_materials 
                << "found: " << mat_cnt << std::endl;
      exit(EXIT_FAILURE);
    }    
    int mat_num = atoi( mat_descr->GetText() );
    if(mat_num != mat_cnt)
    {
      std::cerr << "Error.  Materials not entered in order" << std::endl;
      exit(EXIT_FAILURE);
    }
    
    TiXmlElement* scat_opacity_type = mat_descr->FirstChildElement( "Scattering_Opacity_type");
    TiXmlElement* abs_opacity_type = mat_descr->FirstChildElement( "Absorption_Opacity_type");
    TiXmlElement* cv_type = mat_descr->FirstChildElement( "Cv_type");
    TiXmlElement* rad_source_type = mat_descr->FirstChildElement( "Radiation_Fixed_Source_type");
    TiXmlElement* temp_source_type = mat_descr->FirstChildElement( "Temperature_Fixed_Source_type");
    if(!scat_opacity_type)
    {
      std::cerr << "Missing scattering opacity type in material " << mat_num << std::endl;
      exit(EXIT_FAILURE);
    }
    if(!abs_opacity_type)
    {
      std::cerr << "Missing absorption opacity type in material " << mat_num << std::endl;
      exit(EXIT_FAILURE);
    }
    if(!cv_type)
    {
      std::cerr << "Missing cv type in material " << mat_num << std::endl;
      exit(EXIT_FAILURE);
    }
    if(!rad_source_type)
    {
      std::cerr << "Missing radiation source type in material " << mat_num << std::endl;
      exit(EXIT_FAILURE);
    }
    if(!temp_source_type)
    {
      std::cerr << "Missing temperautre source type in material " << mat_num << std::endl;
      exit(EXIT_FAILURE);
    }
        
    std::string scat_opacity_str = scat_opacity_type->GetText();
    std::string abs_opacity_str = abs_opacity_type->GetText();
    std::string cv_str = cv_type->GetText();
    std::string rad_source_str = rad_source_type->GetText();
    std::string temp_source_str = temp_source_type->GetText();
    
    transform(scat_opacity_str.begin() , scat_opacity_str.end() , scat_opacity_str.begin() , toupper);
    transform(abs_opacity_str.begin() , abs_opacity_str.end() , abs_opacity_str.begin() , toupper);
    transform(cv_str.begin() , cv_str.end() , cv_str.begin() , toupper);
    transform(rad_source_str.begin() , rad_source_str.end() , rad_source_str.begin() , toupper);
    transform(temp_source_str.begin() , temp_source_str.end() , temp_source_str.begin() , toupper);
    
    /// set-up / scan for absorption opacity data
    if(abs_opacity_str == "CONSTANT_XS")
    {
      m_material_absorption_opacity_type[mat_num] = CONSTANT_XS;
      TiXmlElement* const_val = abs_opacity_type->FirstChildElement( "Constant_value" );
      if(!const_val)
      {
        std::cerr << "Missing constant_value tag for material: " << mat_num << " absorption opacity\n" ;
        exit(EXIT_FAILURE);
      }
      m_abs_opacity_double_constants_1[mat_num] = atof( const_val->GetText() ); 
    }
    else if(abs_opacity_str == "RATIONAL")
    {
      m_material_absorption_opacity_type[mat_num] = RATIONAL;
      TiXmlElement* mult_val = abs_opacity_type->FirstChildElement( "Multiplier" );
      TiXmlElement* denom_power_val = abs_opacity_type->FirstChildElement( "Denominator_power" );
      TiXmlElement* denom_offset_val = abs_opacity_type->FirstChildElement( "Denominator_offset" );
      if(!mult_val)
      {
        std::cerr << "Missing Multiplier tag for material " << mat_num << " RATIONAL absorption opacity\n";
        exit(EXIT_FAILURE);
      }
      if(!denom_power_val)
      {
        std::cerr << "Missing Denominator_power tag for material " << mat_num << " RATIONAL absorption opacity\n";
        exit(EXIT_FAILURE);
      }
      if(!denom_offset_val)
      {
        std::cerr << "Missing Denominator_offset tag for material " << mat_num << " RATIONAL absorption opacity\n";
        exit(EXIT_FAILURE);
      }
      m_abs_opacity_double_constants_1[mat_num] = atof(mult_val->GetText() );
      m_abs_opacity_double_constants_2[mat_num] = atof(denom_offset_val->GetText() );
      m_abs_opacity_integer_constants[mat_num] = atoi(denom_power_val->GetText() );
    }
    else if(abs_opacity_str == "TABLE_LOOKUP")
    {
      m_material_absorption_opacity_type[mat_num] = TABLE_LOOKUP;
      TiXmlElement* abs_op_file = abs_opacity_type->FirstChildElement( "File_name" );
      if(!abs_op_file)
      {
        std::cerr << "Missing File_name tag for material " << mat_num << " TABLE_LOOKUP absorption opacity " << std::endl;
        exit(EXIT_FAILURE);
      }
      m_abs_opacity_str[mat_num] = abs_op_file->GetText();
    }
    
    if(m_material_absorption_opacity_type[mat_num] == INVALID_OPACITY_TYPE)
    {
      std::cerr << "Error.  Invalid absorption opacity type for material " << mat_num << std::endl;
      exit(EXIT_FAILURE);    
    }
    
    /// set-up / scan for scattering opacity data
    if(scat_opacity_str == "CONSTANT_XS")
    {
      m_material_scattering_opacity_type[mat_num] = CONSTANT_XS;
      TiXmlElement* const_val = scat_opacity_type->FirstChildElement( "Constant_value" );
      if(!const_val)
      {
        std::cerr << "Missing constant_value tag for material: " << mat_num << " scattering opacity\n" ;
        exit(EXIT_FAILURE);
      }
      m_scat_opacity_double_constants_1[mat_num] = atof( const_val->GetText() ); 
    }
    else if(scat_opacity_str == "RATIONAL")
    {
      m_material_scattering_opacity_type[mat_num] = RATIONAL;
      TiXmlElement* mult_val = scat_opacity_type->FirstChildElement( "Multiplier" );
      TiXmlElement* denom_power_val = scat_opacity_type->FirstChildElement( "Denominator_power" );
      TiXmlElement* denom_offset_val = scat_opacity_type->FirstChildElement( "Denominator_offset" );
      if(!mult_val)
      {
        std::cerr << "Missing Multiplier tag for material " << mat_num << " RATIONAL scattering opacity\n";
        exit(EXIT_FAILURE);
      }
      if(!denom_power_val)
      {
        std::cerr << "Missing Denominator_power tag for material " << mat_num << " RATIONAL scattering opacity\n";
        exit(EXIT_FAILURE);
      }
      if(!denom_offset_val)
      {
        std::cerr << "Missing Denominator_offset tag for material " << mat_num << " RATIONAL scattering opacity\n";
        exit(EXIT_FAILURE);
      }
      m_scat_opacity_double_constants_1[mat_num] = atof(mult_val->GetText() );
      m_scat_opacity_double_constants_2[mat_num] = atof(denom_offset_val->GetText() );
      m_scat_opacity_integer_constants[mat_num] = atoi(denom_power_val->GetText() );
    }
    else if(scat_opacity_str == "TABLE_LOOKUP")
    {
      m_material_scattering_opacity_type[mat_num] = TABLE_LOOKUP;
      TiXmlElement* scat_op_file = scat_opacity_type->FirstChildElement( "File_name" );
      if(!scat_op_file)
      {
        std::cerr << "Missing File_name tag for material " << mat_num << " TABLE_LOOKUP scattering opacity " << std::endl;
        exit(EXIT_FAILURE);
      }
      m_scat_opacity_str[mat_num] = scat_op_file->GetText();
    }
    
    if(m_material_scattering_opacity_type[mat_num] == INVALID_OPACITY_TYPE)
    {
      std::cerr << "Error.  Invalid absorption opacity type for material " << mat_num << std::endl;
      exit(EXIT_FAILURE);    
    }
    
    if(rad_source_str == "NO_SOURCE")
      m_material_radiation_source_type[mat_num] = NO_SOURCE;
   
    if(m_material_radiation_source_type[mat_num] == INVALID_FIXED_SOURCE_TYPE)
    {
      std::cerr << "Error.  Invalid radiation source type for material " << mat_num << std::endl;
      exit(EXIT_FAILURE);    
    }
    
    if(temp_source_str == "NO_SOURCE")
      m_material_temperature_source_type[mat_num] = NO_SOURCE;
   
    if(m_material_temperature_source_type[mat_num] == INVALID_FIXED_SOURCE_TYPE)
    {
      std::cerr << "Error.  Invalid temperature source type for material " << mat_num << std::endl;
      exit(EXIT_FAILURE);    
    }
    
    if(cv_str == "CONSTANT_CV")
    {
      m_material_cv_type[mat_num] = CONSTANT_CV;
      TiXmlElement* cv_const = cv_type->FirstChildElement( "Cv_constant" );
      if(!cv_const)
      {
        std::cerr << "Error.  Missing Cv_constant tag in material " << mat_num << std::endl;
        exit(EXIT_FAILURE);
      }
      m_cv_constants[mat_num] = atof(cv_const->GetText() );
      if(m_cv_constants[mat_num] < 0. )
      {
        std::cerr << "Invalid Cv in material "<< mat_num << " values must be positive floats\n";
        exit(EXIT_FAILURE);
      }
    }
    
    if(m_material_cv_type[mat_num] == INVALID_CV_TYPE)
    {
      std::cerr << "Error.  Invalid cv type for material " << mat_num << std::endl;
      exit(EXIT_FAILURE);    
    }
    
    mat_descr = mat_descr->NextSiblingElement("Material");
  }
  
  return 0;
}

int Input_Reader::load_time_stepping_data(TiXmlElement* time_elem)
{
  //! All of the required elements in the time block
  TiXmlElement* dt_min_elem = time_elem->FirstChildElement( "Dt_min");
  TiXmlElement* dt_max_elem = time_elem->FirstChildElement( "Dt_max");
  TiXmlElement* t_start_elem = time_elem->FirstChildElement( "T_start");
  TiXmlElement* t_end_elem = time_elem->FirstChildElement( "T_end");
  TiXmlElement* solver_elem = time_elem->FirstChildElement( "Time_solver");
  TiXmlElement* start_meth_elem = time_elem->FirstChildElement( "Starting_method");
  
  if(!dt_min_elem)
  {
    std::cerr << "Error.  Missing Dt_min element in TIME input block" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(!dt_max_elem)
  {
    std::cerr << "Error.  Missing Dt_max element in TIME input block" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(!t_start_elem)
  {
    std::cerr << "Error.  Missing T_start element in TIME input block" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(!t_end_elem)
  {
    std::cerr << "Error.  Missing T_end element in TIME input block" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(!solver_elem)
  {
    std::cerr << "Error.  Missing Time_solver element in TIME input block" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(!start_meth_elem)
  {
    std::cerr << "Error.  Missing Starting_method element in TIME input block" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  m_dt_min = atof(dt_min_elem->GetText() );
  m_dt_max = atof(dt_max_elem->GetText() );
  m_t_start = atof(t_start_elem->GetText() );
  m_t_end = atof(t_end_elem->GetText() );
  
  if( (m_dt_min < 0.) || (m_dt_max < 0.) )
  {
    std::cerr << "Error.  Time step sizes must be positive floats" << std::endl;
    exit(EXIT_FAILURE);
  }
  if( m_dt_min > m_dt_max )
  {
    std::cerr << "Error.  Minimum time step must be greater than maximum time step" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(m_t_start > m_t_end)
  {
    std::cerr << "Error.  t_end must be greater than t_start " << std::endl;
    exit(EXIT_FAILURE);
  }
  
  std::string stepper_str = solver_elem->GetText();
  std::string starter_str = start_meth_elem->GetText();
  transform(stepper_str.begin() , stepper_str.end() , stepper_str.begin() , toupper);
  transform(starter_str.begin() , starter_str.end() , starter_str.begin() , toupper);

  if(stepper_str == "IMPLICIT_EULER")
    m_time_step_scheme = IMPLICIT_EULER;
  
  if(starter_str == "EXPONENTIAL")
    m_time_starting_method = EXPONENTIAL;
  else if(starter_str == "VECTOR")
    m_time_starting_method = VECTOR;
  else if(starter_str == "RAMP")
    m_time_starting_method = RAMP;
    
  if(m_time_step_scheme == INVALID_TIME_SOLVER)
  {
    std::cerr << "Error.  Invalid time stepping scheme" << std::endl;
    exit(EXIT_FAILURE);  
  }
  
  switch (m_time_starting_method)
  {
    case INVALID_STARTING_METHOD:
    {
      std::cerr << "Error.  Invalid time iteration starting scheme" << std::endl;
      exit(EXIT_FAILURE);  
    
      break;
    }
    case EXPONENTIAL:
    {
      /// Get increase factor
      TiXmlElement* exp_incr_fact = time_elem->FirstChildElement( "Increase_factor");
      if(!exp_incr_fact)
      {
        std::cerr << "Error.  Missing Increase_factor element.  Required for Starting_method EXPONENTIAL" << std::endl;
        exit(EXIT_FAILURE);
      }
     
      m_exponential_ratio = atof( exp_incr_fact->GetText() );
      if(m_exponential_ratio < 1.)
      {
        std::cerr << "Error.  Starting time step increase must be greater than 1.0" << std::endl;
        exit(EXIT_FAILURE);
      }
    
      break;
    }
    case RAMP:
    {
      TiXmlElement* ramp_steps = time_elem->FirstChildElement( "Ramp_steps" );
      if(!ramp_steps)
      {
        std::cerr << "Error.  RAMP time starter requires Ramp_steps element\n" ;
        exit(EXIT_FAILURE);
      }
      
      m_ramp_steps = atoi(ramp_steps->GetText() );
      
      if(m_ramp_steps < 1)
      {
        std::cerr << "Error.  RAMP starter requires at least one time step before full time step can be taken\n";
        exit(EXIT_FAILURE);
      }      
      break;
    }
    case VECTOR:
    {
      TiXmlElement* n_stages = time_elem->FirstChildElement( "N_vector_stages" );
      if(!n_stages)
      {
        std::cerr << "Error.  VECTOR time starter requires N_stages element" << std::endl;
        exit(EXIT_FAILURE);
      }
      int num_vec_stages = atoi(n_stages->GetText() );
      if(num_vec_stages < 1)
      {
        std::cerr << "Error.  Must have at least 1 stage of small time steps with VECTOR time starter" << std::endl;
        exit(EXIT_FAILURE);
      }
      
      m_vector_start_sizes.resize(num_vec_stages , 0.);
      m_vector_start_step_numbers.resize(num_vec_stages , 0);
      
      TiXmlElement* vector_stage = time_elem->FirstChildElement("Vector_stage");

      for(int i=0; i<num_vec_stages ; i++)
      {
        if(!vector_stage)
        {
          std::cerr << "Error.  Expected " << num_vec_stages << " Vector_stage elements" << std::endl;
          std::cerr << "Only found: " << i << " Vector_stage elements" << std::endl;
          exit(EXIT_FAILURE);
        }
        
        int curr_stage = atoi(vector_stage->GetText() );
        if(curr_stage < 0 || curr_stage == num_vec_stages)
        {
          std::cerr << "Error.  Vector_stage number must be between 0 and n_stages - 1" << std::endl;
          exit(EXIT_FAILURE);
        }
      
        TiXmlElement* stage_steps = vector_stage->FirstChildElement("Stage_steps");
        TiXmlElement* stage_factor = vector_stage->FirstChildElement("Stage_size");

        if(!stage_steps || !stage_factor)
        {
          std::cerr << "Expect a Stage_steps and Stage_size element in each Vector_stage element" << std::endl;
          exit(EXIT_FAILURE);
        }
        
        int stages = atoi( stage_steps->GetText() );
        double size_factor = atof( stage_factor->GetText() );
        
        if(stages < 1)
        {
          std::cerr << "Error. Stage_steps must be an integer greater than 1" << std::endl;
          exit(EXIT_FAILURE);       
        }
        if(size_factor <= 0.  || size_factor >= 1.)
        {
          std::cerr << "Error.  Stage_size must be a double between 0 and 1 " << std::endl;
          exit(EXIT_FAILURE);      
        }
        m_vector_start_step_numbers[curr_stage] = stages;
        m_vector_start_sizes[curr_stage] = size_factor;      
        
        vector_stage = vector_stage->NextSiblingElement( "Vector_stage");
      }
      break;
    }
  }
    
  return 0;
}

int Input_Reader::load_spatial_discretization_data(TiXmlElement* spatial_element)
{
  /// all of the required elements in the spatial_discretization input block
  TiXmlElement* trial_space_elem = spatial_element->FirstChildElement( "DFEM_degree");
  TiXmlElement* int_type_elem = spatial_element->FirstChildElement( "Integration_type");
  TiXmlElement* dfem_interp_point_type_elem = spatial_element->FirstChildElement( "DFEM_interpolation_point_type");
  TiXmlElement* opacity_treatment_elem = spatial_element->FirstChildElement( "Opacity_treatment");
  TiXmlElement* opacity_interp_point_type_elem = spatial_element->FirstChildElement( "Opacity_interpolation_point_type");

  if(!trial_space_elem)
  {
    std::cerr << "Error.  Missing DFEM_degree element from SPATIAL_DISCRETIZATION block" << std::endl;
    exit(EXIT_FAILURE);  
  }
  if(!int_type_elem)
  {
    std::cerr << "Error.  Missing Integration_type element from SPATIAL_DISCRETIZATION block" << std::endl;
    exit(EXIT_FAILURE);  
  }  
  if(!dfem_interp_point_type_elem)
  {
    std::cerr << "Error.  Missing DFEM_interpolation_point_type element from SPATIAL_DISCRETIZATION block" << std::endl;
    exit(EXIT_FAILURE);  
  }  
  if(!opacity_treatment_elem)
  {
    std::cerr << "Error.  Missing Opacity_treatment element from SPATIAL_DISCRETIZATION block" << std::endl;
    exit(EXIT_FAILURE);  
  }  
  if(!opacity_interp_point_type_elem)
  {
    std::cerr << "Error.  Missing Opacity_interp_point_type_elem element from SPATIAL_DISCRETIZATION block" << std::endl;
    exit(EXIT_FAILURE);  
  }
  
  m_dfem_trial_space_degree = atoi( trial_space_elem->GetText() );
  if(m_dfem_trial_space_degree < 1)
  {
    std::cerr << "Error.  Trial space degree must be a positive integrer greater than 0" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  std::string int_type_str = int_type_elem->GetText();
  transform(int_type_str.begin() , int_type_str.end() , int_type_str.begin() , toupper);
  if(int_type_str == "SELF_LUMPING")
    m_integration_type = SELF_LUMPING;
  else if(int_type_str == "TRAD_LUMPING")
    m_integration_type = TRAD_LUMPING;
  else if(int_type_str == "EXACT")
    m_integration_type = EXACT;
    
  if( m_integration_type == INVALID_MATRIX_INTEGRATION)
  {
    std::cerr << "Invalid matrix integration strategy" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  std::string dfem_intep_type_str = dfem_interp_point_type_elem->GetText();
  transform(dfem_intep_type_str.begin() , dfem_intep_type_str.end() , dfem_intep_type_str.begin() , toupper);
  if(dfem_intep_type_str == "GAUSS")
    m_dfem_interpolation_point_type = GAUSS;
  else if(dfem_intep_type_str == "LOBATTO")
    m_dfem_interpolation_point_type = LOBATTO;
  else if(dfem_intep_type_str == "EQUAL_SPACED")
    m_dfem_interpolation_point_type = EQUAL_SPACED;
    
  if( m_dfem_interpolation_point_type == INVALID_QUADRATURE_TYPE)
  {
    std::cerr << "Invalid DFEM interpolation point type" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  std::string opacity_intep_type_str = opacity_interp_point_type_elem->GetText();
  transform(opacity_intep_type_str.begin() , opacity_intep_type_str.end() , opacity_intep_type_str.begin() , toupper);
  if(opacity_intep_type_str == "GAUSS")
    m_opacity_interpolation_point_type = GAUSS;
  else if(opacity_intep_type_str == "LOBATTO")
    m_opacity_interpolation_point_type = LOBATTO;
  else if(opacity_intep_type_str == "EQUAL_SPACED")
    m_opacity_interpolation_point_type = EQUAL_SPACED;
    
  if( m_dfem_interpolation_point_type == INVALID_QUADRATURE_TYPE)
  {
    std::cerr << "Invalid DFEM interpolation point type" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  std::string opacity_treat_str = opacity_treatment_elem->GetText();
  transform(opacity_treat_str.begin() , opacity_treat_str.end() , opacity_treat_str.begin() , toupper);
  if(opacity_treat_str == "MOMENT_PRESERVING")
    m_opacity_treatment = MOMENT_PRESERVING;
  else if(opacity_treat_str == "SLXS" )
    m_opacity_treatment = SLXS;
  else if(opacity_treat_str == "INTERPOLATING")
    m_opacity_treatment = INTERPOLATING;
    
  if(m_opacity_treatment == INVALID_OPACITY_TREATMENT)
  {
    std::cerr << "Errror. Invalid opacity treatment" << std::endl;
    exit(EXIT_FAILURE);    
  }
  
  /// if not using self lumping treatment, need to get degree of opacity spatial treatment
  if(m_opacity_treatment != SLXS)
  {
    TiXmlElement* opacity_degree_elem = spatial_element->FirstChildElement( "Opacity_polynomial_degree");
    if(!opacity_degree_elem)
    {
      std::cerr << "Error.  Missing Opacity_degree_elem element from SPATIAL_DISCRETIZATION block" << std::endl;
      exit(EXIT_FAILURE);  
    }
    
    m_opacity_polynomial_degree = atoi( opacity_degree_elem->GetText() );
    if(m_opacity_polynomial_degree < 0)
    {
      std::cerr << "Error.  Invalid opacity polynomial" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  

  
  return 0;  
}

  int Input_Reader::load_angular_discretization_data(TiXmlElement* angle_element)
{
  TiXmlElement* n_angle_elem = angle_element->FirstChildElement( "Number_of_angles");
  TiXmlElement* n_group_elem = angle_element->FirstChildElement( "Number_of_groups");
  TiXmlElement* quad_type_elem = angle_element->FirstChildElement( "Quadrature_type");
  TiXmlElement* n_legendre_mom_elem = angle_element->FirstChildElement( "Number_of_legendre_moments");
  
  
  if(!n_angle_elem)
  {
    std::cerr << "Error.  Missing Number_of_angles element" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  if(!n_group_elem)
  {
    std::cerr << "Error.  Missing Number_of_groups element" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  if(!quad_type_elem)
  {
    std::cerr << "Error.  Missing Quadrature_type element" << std::endl;
    exit(EXIT_FAILURE);
  }
  
   if(!n_legendre_mom_elem)
  {
    std::cerr << "Error.  Missing Number_of_legendre_moments element" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  m_number_angles = atoi( n_angle_elem->GetText() );
  if( m_number_angles < 0 || ((m_number_angles % 2) != 0) )
  {
    std::cerr << "Error.  Number of Angles must be a positive, even integer" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  m_number_groups = atoi(n_group_elem->GetText() ) ;
  if( m_number_groups < 1)
  {
    std::cerr << "Error.  Number of groups must be an integer > 0" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  m_group_lower_bounds.resize(m_number_groups,0.);
  m_group_upper_bounds.resize(m_number_groups,0.);
  
  /// load group bounds
  if( m_number_groups > 1)
  {
    TiXmlElement* grp_bounds = angle_element->FirstChildElement("Group_Boundaries");
    for(int edge_cnt = 0; edge_cnt <= m_number_groups ; edge_cnt++)
    {    
      if(!grp_bounds)
      {
        std::cerr << "Error.  Missing Group_Boundaries element in ANGULAR_DISCRETIZATION block.  Expected: " << m_number_groups + 1 
                  << " found: " << edge_cnt << std::endl;
        exit(EXIT_FAILURE);
      }    
      int edge_num = atoi( grp_bounds->GetText() );
      if(edge_num != edge_cnt)
      {
        std::cerr << "Error.  Group Bounds not entered in order" << std::endl;
        exit(EXIT_FAILURE);
      }
      TiXmlElement* grp_edge_val = grp_bounds->FirstChildElement("Edge_Value");
      double edge = atof( grp_edge_val->GetText() );
      
      /// the lower the group number, the higher the average frequency group energy (by convention)
      if(edge_cnt < m_number_groups)
        m_group_upper_bounds[edge_cnt] = edge;
        
      if(edge_cnt > 0)
        m_group_lower_bounds[edge_cnt - 1] = edge;      
    
      grp_bounds = grp_bounds->NextSiblingElement("Group_Boundaries");
    }
    
    /// check that these are logical values
    for(int g = 1; g<m_number_groups ; g++)
    {
      if( m_group_lower_bounds[g] > m_group_lower_bounds[g-1] )
      {
        std::cerr << "Frequency groups lower energy bounds are not in descending order\n" ;
        std::cerr << "Problem lies between group: " << g << " and group: " << g-1 << std::endl ;
        exit(EXIT_FAILURE);
      }
      if( m_group_upper_bounds[g] >  m_group_upper_bounds[g-1] )
      {
        std::cerr << "Frequency group upper energy bounds are not in descending order\n";
        std::cerr << "Problem lies between group: " << g << " and group: " << g-1 << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  }
  
  std::string ang_quad_str = quad_type_elem->GetText();
  // change to all upper case 
  transform(ang_quad_str.begin() , ang_quad_str.end() , ang_quad_str.begin() , toupper);
  if(ang_quad_str == "GAUSS_ANGLE")
    m_angular_quadrature_type = GAUSS_ANGLE;
  else if(ang_quad_str == "LOBATTO_ANGLE")
    m_angular_quadrature_type = LOBATTO_ANGLE;

  if(m_angular_quadrature_type == INVALID_ANGULAR_QUADRATURE_TYPE)
  {
    std::cerr << "Error. Invalid Quadrature_type in ANGULAR_DISCRETIZATION" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  m_n_leg_moments = atoi( n_legendre_mom_elem->GetText() );
  if(m_n_leg_moments < 1)
  {
    std::cerr << "Error.  Number of legendre moments must be >= 1\n";
    exit(EXIT_FAILURE);
  }
    
  return 0;
}

int Input_Reader::load_solver_data(TiXmlElement* solver_element)
{
  TiXmlElement* solver_type_elem = solver_element->FirstChildElement( "WG_Solver_type");
  TiXmlElement* wg_tolerance_elem = solver_element->FirstChildElement( "WG_Tolerance");  
  
  if(!solver_type_elem)
  {
    std::cerr << "Missing WG_Solver_type element\n" ;
    exit(EXIT_FAILURE);  
  }
  
  if(!wg_tolerance_elem)
  {
    std::cerr << "Missing WG_Tolerance element\n" ;
    exit(EXIT_FAILURE);  
  }
   
  m_wg_tolerance = atof( wg_tolerance_elem->GetText() );
  if( (m_wg_tolerance < 1.E-15) || (m_wg_tolerance> 1.E-4))
  {
    std::cerr << "Invalid within group tolerance.  Must be greater 1E-15 and less than 1E-4\n";
    exit(EXIT_FAILURE);
  }
    
  std::string solve_type_str = solver_type_elem->GetText();
  // change to all upper case 
  transform(solve_type_str.begin() , solve_type_str.end() , solve_type_str.begin() , toupper);
  if(solve_type_str == "FP_SWEEPS")
    m_wg_solve_type = FP_SWEEPS;
  else if (solve_type_str == "FP_DSA")
    m_wg_solve_type = FP_SWEEPS;
  else if (solve_type_str == "KRYLOV_SWEEPS")
    m_wg_solve_type = KRYLOV_SWEEPS;
  else if (solve_type_str == "KRYLOV_DSA")
    m_wg_solve_type = KRYLOV_DSA;
  
  if(m_wg_solve_type == INVALID_WG_SOLVE_TYPE)
  {
    std::cerr << "Invalid within group radiation solver type\n";
    exit(EXIT_FAILURE);
  }  
  
  if( (m_wg_solve_type == FP_SWEEPS) || (m_wg_solve_type == FP_DSA))
  {
    TiXmlElement* num_sweep_elem = solver_type_elem->FirstChildElement( "Max_Within_Group_Sweeps");
    if(!num_sweep_elem)
    {
      std::cerr << "Missing Max_Within_Group_Sweeps element.\n" ;
      exit(EXIT_FAILURE);
    }
    m_max_num_sweeps = atoi( num_sweep_elem->GetText() );
    if(m_max_num_sweeps < 1)
    {
      std::cerr << "Must allow at least one sweep per within group solve\n";
      exit(EXIT_FAILURE);
    }
  }
  
  /**
      Things that are only required / neeeded for multi frequency problems
  */
  if(m_number_groups > 1)
  {
    TiXmlElement* bg_tolerance_elem = solver_element->FirstChildElement( "BG_Tolerance");
    
    /// only require the between group solver tolerance if number of groups is greater than 1
    if(!bg_tolerance_elem)
    {
      std::cerr << "MF problems require specification of between group tolerance, BG_Tolerance.  Element not found\n";
      exit(EXIT_FAILURE);
    }
    
    m_bg_tolerance = atof(bg_tolerance_elem->GetText() );
    if( m_bg_tolerance < m_wg_tolerance)
    {
      std::cerr << "Invalid between group tolerance.  Must be greater than the within group scattering tolerance.\n";
      exit(EXIT_FAILURE);
    }
    if( m_bg_tolerance > 1.E-4)
    {
      std::cerr << "Invalid between group tolerance.  Must be less than 1.0E-4.\n";
      exit(EXIT_FAILURE);
    }
    
    TiXmlElement* mf_ard_solve_elem = solver_element->FirstChildElement( "MF_Solver_Type" );
    if(!mf_ard_solve_elem)
    {
      std::cerr<< "MF_Solver_Type element not found in multi-frequency problem\n";
      exit(EXIT_FAILURE);
    } 
    
    std::string ard_solve_type_str = mf_ard_solve_elem->GetText();
    // change to all upper case 
    transform(ard_solve_type_str.begin() , ard_solve_type_str.end() , ard_solve_type_str.begin() , toupper);
    if(ard_solve_type_str == "FP_NO_ACCEL")
      m_ard_solve_type = FP_NO_ACCEL;
    else if (ard_solve_type_str == "FP_LMFGA")
      m_ard_solve_type = FP_LMFGA;
    else if (ard_solve_type_str == "KRYLOV_LMFGA")
      m_ard_solve_type = KRYLOV_LMFGA;
    
    if(m_ard_solve_type == INVALID_ARD_SOLVE_TYPE)
    {
      std::cerr << "Invalid ARD_SOLVE_TYPE found\n";
      exit(EXIT_FAILURE);
    }
    
    if( (m_ard_solve_type == FP_NO_ACCEL) || (m_ard_solve_type == FP_LMFGA) )
    {
      TiXmlElement* fp_ard_iter_elem = mf_ard_solve_elem->FirstChildElement("Max_Iterations");
      
      if(!fp_ard_iter_elem)
      {
        std::cerr << "Fixed point ARD solvers require Max_Iterations element\n";
        exit(EXIT_FAILURE);
      }
      
      m_max_ard_iterations = atoi( fp_ard_iter_elem->GetText() );
      if(m_max_ard_iterations < 1)
      {
        std::cerr << "Require at least 1 iteration for FP ARD solver schemes\n";
        exit(EXIT_FAILURE);
      }
    }
  }
  
  /// get thermal convergence tolerance
  TiXmlElement* thermal_tolerance_elem = solver_element->FirstChildElement( "Thermal_Tolerance");
  if(!thermal_tolerance_elem)
  {
    std::cerr << "Thermal_Tolerance element required \n";
    exit(EXIT_FAILURE);
  }
  m_thermal_tolerance = atof( thermal_tolerance_elem->GetText() );
  if(m_thermal_tolerance < 1.E-14)
  {
    std::cerr << "Too small of thermal tolerance. \n" ; 
    exit(EXIT_FAILURE);    
  }
  if( m_number_groups > 1)
  {
    if( m_bg_tolerance < m_wg_tolerance)
    {
      std::cerr << "Between group tolerance must be greater within group tolerance \n";
      exit(EXIT_FAILURE);
    }
    if( m_thermal_tolerance < m_bg_tolerance)
    {
      std::cerr << "Thermal tolerance must be greater than between group absorption/re-emission tolerance \n";
      exit(EXIT_FAILURE);
    }
  }
  else
  {
    if(m_thermal_tolerance < m_wg_tolerance)
    {
      std::cerr << "Grey problem.  Within group tolerance must be smaller than thermal iteration tolerance\n";
      exit(EXIT_FAILURE);    
    }
  }
  
  return 0;
}

int Input_Reader::load_bc_ic_data(TiXmlElement* bc_ic_element)
{
    // TEMPERATURE_IC_TYPE m_temp_ic_type = INVALID_TEMPERATURE_IC_TYPE;
  // std::vector<double> m_region_temperature;
  // RADIATION_IC_TYPE m_radiation_ic_type  = INVALID_RADIATION_IC_TYPE;
  // RADIATION_BC_TYPE m_rad_bc_left  = INVALID_RADIATION_BC_TYPE;
  // RADIATION_BC_TYPE m_rad_bc_right = INVALID_RADIATION_BC_TYPE;
  TiXmlElement* temp_ic_type_elem = bc_ic_element->FirstChildElement( "Temperature_IC_Type");
  TiXmlElement* rad_ic_type_elem = bc_ic_element->FirstChildElement( "Radiation_IC_Type");
  TiXmlElement* rad_left_bc_type_elem = bc_ic_element->FirstChildElement( "Left_Radiation_BC_Type");
  TiXmlElement* rad_right_bc_type_elem = bc_ic_element->FirstChildElement( "Right_Radiation_BC_Type");
  
  if(!temp_ic_type_elem)
  {
    std::cerr << "Missing Temperature_IC_Type element\n" ;
    exit(EXIT_FAILURE);
  }
  
  if(!rad_ic_type_elem)
  {
    std::cerr << "Missing Radiation_IC_Type element\n" ;
    exit(EXIT_FAILURE);
  }
  
  if(!rad_left_bc_type_elem)
  {
    std::cerr << "Missing Left_Radiation_BC_Type element\n" ;
    exit(EXIT_FAILURE);
  }
  
  if(!rad_right_bc_type_elem)
  {
    std::cerr << "Missing _Right_Radiation_BC_Type element\n" ;
    exit(EXIT_FAILURE);
  }
  

  
  /// Get IC strings and set IC types
  std::string temp_ic_str = temp_ic_type_elem->GetText();
  transform(temp_ic_str.begin() , temp_ic_str.end() , temp_ic_str.begin() , toupper);
  std::string rad_ic_str = rad_ic_type_elem->GetText();
  transform(rad_ic_str.begin() , rad_ic_str.end() , rad_ic_str.begin() , toupper);
  
  if(temp_ic_str == "CONSTANT_TEMPERATURE_IC")
  {
    m_temperature_ic_type = CONSTANT_TEMPERATURE_IC;
  }
  
  if(rad_ic_str == "PLANCKIAN_IC")
  {
    m_radiation_ic_type = PLANCKIAN_IC;
  }  
  
  /// Get BC strings and set BC types
  std::string rad_bc_left_str = rad_left_bc_type_elem->GetText();
  transform(rad_bc_left_str.begin() , rad_bc_left_str.end() , rad_bc_left_str.begin() , toupper);
  std::string rad_bc_right_str = rad_right_bc_type_elem->GetText();
  transform(rad_bc_right_str.begin() , rad_bc_right_str.end() , rad_bc_right_str.begin() , toupper);
  
  if(rad_bc_left_str == "VACUUM")
  {
    m_rad_bc_left = VACUUM;
  }
  else if(rad_bc_left_str == "PLANCKIAN_BC")
  {
    m_rad_bc_left = PLANCKIAN_BC;
  }
  else if(rad_bc_left_str == "REFLECTIVE")
  {
    m_rad_bc_left = REFLECTIVE;
  }
  
  if(rad_bc_right_str == "VACUUM")
  {
    m_rad_bc_right = VACUUM;
  }
  else if(rad_bc_right_str == "PLANCKIAN_BC")
  {
    m_rad_bc_right = PLANCKIAN_BC;
  }
  else if(rad_bc_right_str == "REFLECTIVE")
  {
    m_rad_bc_right = REFLECTIVE;
  }
  
  
  /// Handle left radiation boundary condition
  if( m_rad_bc_left == INVALID_RADIATION_BC_TYPE)
  {
    std::cerr << "Left radiation BC type not recognized\n";
    exit(EXIT_FAILURE);
  }
  else if( m_rad_bc_left == VACUUM)
  {
    m_left_bc_value = 0.; 
  }
  else if( m_rad_bc_left == PLANCKIAN_BC )
  {
    TiXmlElement* rad_left_bc_value_elem = rad_left_bc_type_elem->FirstChildElement( "BC_Temperature");
    if(!rad_left_bc_value_elem)
    {
      std::cerr << "Missing BC_Temperature element in Left Radiation Planckian BC Specification\n";
      exit(EXIT_FAILURE);
    }
    m_left_bc_value = atof( rad_left_bc_value_elem->GetText() ); 
    if(m_left_bc_value < 0.)
    {
      std::cerr << "Invalid value for Left BC temperature.  Must be >= 0\n";
      exit(EXIT_FAILURE);
    }
  }
  else if( m_rad_bc_left== REFLECTIVE)
  {
    /// don't need any additional data
  }
  
  /// Right radiation boundary condition
  if( m_rad_bc_right == INVALID_RADIATION_BC_TYPE)
  {
    std::cerr << "Right radiation BC type not recognized\n";
    exit(EXIT_FAILURE);
  }
  else if( m_rad_bc_right == VACUUM)
  {
    m_right_bc_value = 0.; 
  }
  else if( m_rad_bc_right == PLANCKIAN_BC )
  {
    TiXmlElement* rad_right_bc_value_elem = rad_right_bc_type_elem->FirstChildElement( "BC_Temperature");
    if(!rad_right_bc_value_elem)
    {
      std::cerr << "Missing BC_Temperature element in Right Radiation Planckian BC Specification\n";
      exit(EXIT_FAILURE);
    }
    m_right_bc_value = atof( rad_right_bc_value_elem->GetText() ); 
    if(m_right_bc_value < 0.)
    {
      std::cerr << "Invalid value for Right BC temperature.  Must be >= 0\n";
      exit(EXIT_FAILURE);
    }
  }
  else if( m_rad_bc_right == REFLECTIVE)
  {
    std::cerr << "Cannot use reflective boundary condition on right edge.  Only left edge\n";
    exit(EXIT_FAILURE);
  }
    
  /// get radiation initial conditions
  if( m_radiation_ic_type == INVALID_RADIATION_IC_TYPE)
  {
    std::cerr << "Radiation IC type not recognized\n";
    exit(EXIT_FAILURE);
  }
  else if( m_radiation_ic_type == PLANCKIAN_IC )
  {
    m_region_radiation_temperature.resize(m_number_regions,0.);
    TiXmlElement* rad_ic_reg_elem = rad_ic_type_elem->FirstChildElement( "Region");
    for(int reg = 0; reg < m_number_regions; reg++)
    {
      if(!rad_ic_reg_elem)
      {
        std::cerr << "Missing a region radiation temperature for every region.  Found: " << reg << " Need: " << m_number_regions << std::endl;
        exit(EXIT_FAILURE);
      }
      if( atoi(rad_ic_reg_elem->GetText()) != reg)
      {
        std::cerr << "Temperature IC for regions not in order\n";
        exit(EXIT_FAILURE);
      }
      
      
      TiXmlElement* rad_temp_value = rad_ic_reg_elem->FirstChildElement("Radiation_Temperature");
      if(!rad_temp_value)
      {
        std::cerr << "Missing required Radiation_Temperature element for region: " << reg << std::endl;
        exit(EXIT_FAILURE);
      }
      m_region_radiation_temperature[reg] = atof( rad_temp_value->GetText() );
      if(m_region_radiation_temperature[reg] < 0.)
      {
        std::cerr << "Invalid radiation temperature in region: " << reg << " , must be >= 0. \n";
        exit(EXIT_FAILURE);
      }
      
      rad_ic_reg_elem->NextSiblingElement("Region");
    }
  }  
  /// Get temperature initial conditions
  if( m_temperature_ic_type == INVALID_TEMPERATURE_IC_TYPE)
  {
    std::cerr << "Temperature IC type not recognized\n";
    exit(EXIT_FAILURE);
  }
  else if( m_temperature_ic_type == CONSTANT_TEMPERATURE_IC )
  {
    m_region_temperature.resize(m_number_regions,0.);
    TiXmlElement* temp_ic_reg_elem = temp_ic_type_elem->FirstChildElement( "Region");
    for(int reg = 0; reg < m_number_regions; reg++)
    {
      if(!temp_ic_reg_elem)
      {
        std::cerr << "Missing a region material temperature for every region.  Found: " << reg << " Need: " << m_number_regions << std::endl;
        exit(EXIT_FAILURE);
      }
      if( atoi(temp_ic_reg_elem->GetText()) != reg)
      {
        std::cerr << "Temperature IC for regions not in order\n";
        exit(EXIT_FAILURE);
      }
      
      
      TiXmlElement* temp_value = temp_ic_reg_elem->FirstChildElement("Material_Temperature");
      if(!temp_value)
      {
        std::cerr << "Missing required Material_Temperature element for region: " << reg << std::endl;
        exit(EXIT_FAILURE);
      }
      m_region_temperature[reg] = atof( temp_value->GetText() );
      if(m_region_temperature[reg] < 0.)
      {
        std::cerr << "Invalid material temperature in region: " << reg << " , must be >= 0. \n";
        exit(EXIT_FAILURE);
      }
      
      temp_ic_reg_elem->NextSiblingElement("Region");
    }
  }
  
  // m_bc_start_time; m_bc_end_time;
  /// check and see if dirichlet boundary conditions are only for part of time
  if( (m_rad_bc_left == PLANCKIAN_BC ) || (m_rad_bc_right == PLANCKIAN_BC) )
  {
    TiXmlElement* bc_time_dependent_elem = bc_ic_element->FirstChildElement( "BC_Time_Dependent");
    if(!bc_time_dependent_elem )
    {
      std::cerr << "Missing BC_Time_Dependence block, required\n";
      exit(EXIT_FAILURE);
    }
    
    std::string bc_time_str = bc_time_dependent_elem->GetText();
    transform(bc_time_str.begin() , bc_time_str.end() , bc_time_str.begin() , toupper);
    if(bc_time_str == "BC_CONSTANT")
    {
      m_bc_time_dependence = BC_CONSTANT;
    }
    else if(bc_time_str == "BC_BURST")
    {
      m_bc_time_dependence = BC_BURST;
    }
    
    if(m_bc_time_dependence == INVALID_BC_TIME_DEPENDENCE)
    {
      std::cerr << "Invalid Planckian BC time dependence\n";
      exit(EXIT_FAILURE);
    }    
    else if(m_bc_time_dependence == BC_CONSTANT)
    {
      /// assume dirichlet conditions last forever
      m_bc_start_time = m_t_start - 0.01*(m_t_end - m_t_start);
      m_bc_end_time = m_t_end + 0.01*(m_t_end - m_t_start);  
    }
    else if(m_bc_time_dependence == BC_BURST)
    {
      TiXmlElement* bc_start_elem = bc_time_dependent_elem->FirstChildElement( "BC_Turn_On" );
      TiXmlElement* bc_end_elem = bc_time_dependent_elem->FirstChildElement( "BC_Turn_Off" );
      
      if(!bc_start_elem || !bc_end_elem)
      {
        std::cerr << "BC_Time_Dependent block requires BC_Turn_On and BC_Turn_Off elements\n";
        exit(EXIT_FAILURE);
      }
      
      m_bc_start_time = atof( bc_start_elem->GetText() );
      m_bc_end_time = atof( bc_end_elem->GetText() );
      
      if(m_bc_end_time < m_bc_start_time)
      {
        std::cerr << "BC_Turn_Off must be later (in time) than BC_Turn_Off\n";
        exit(EXIT_FAILURE);
      }      
    }
  }
  /// get angular dependence for dirichlet angular BC
  if( m_rad_bc_left == PLANCKIAN_BC )
  {
    TiXmlElement* bc_left_incidence_type = rad_left_bc_type_elem->FirstChildElement("BC_Incidence");
    if(!bc_left_incidence_type)
    {
      std::cerr << "Missing Angular incidence type for Planckian BC on left edge\n";
      exit(EXIT_FAILURE);
    }
    std::string left_bc_incidence_str = bc_left_incidence_type->GetText();
    transform(left_bc_incidence_str.begin() , left_bc_incidence_str.end() , left_bc_incidence_str.begin() , toupper);
    if(left_bc_incidence_str == "BC_ISOTROPIC")
    {
      m_left_bc_angle_dependence = BC_ISOTROPIC;
    }
    else if(left_bc_incidence_str == "BC_GLANCE")
    {
      m_left_bc_angle_dependence = BC_GLANCE;
    }
    else if(left_bc_incidence_str == "BC_NORMAL")
    {
      m_left_bc_angle_dependence = BC_NORMAL;
    }
    if(m_left_bc_angle_dependence == INVALID_BC_ANGLE_DEPENDENCE)
    {
      std::cerr << "Must specify angular dependence of Planckian BC on left edge\n";
      exit(EXIT_FAILURE);
    }
  }
  
  if(m_rad_bc_right == PLANCKIAN_BC)
  {
    TiXmlElement* bc_right_incidence_type = rad_right_bc_type_elem->FirstChildElement("BC_Incidence");
    if(!bc_right_incidence_type)
    {
      std::cerr << "Missing Angular incidence type for Planckian BC on right edge\n";
      exit(EXIT_FAILURE);
    }
    std::string right_bc_incidence_str = bc_right_incidence_type->GetText();
    transform(right_bc_incidence_str.begin() , right_bc_incidence_str.end() , right_bc_incidence_str.begin() , toupper);
    if(right_bc_incidence_str == "BC_ISOTROPIC")
    {
      m_right_bc_angle_dependence = BC_ISOTROPIC;
    }
    else if(right_bc_incidence_str == "BC_GLANCE")
    {
      m_right_bc_angle_dependence = BC_GLANCE;
    }
    else if(right_bc_incidence_str == "BC_NORMAL")
    {
      m_right_bc_angle_dependence = BC_NORMAL;
    }
    if(m_right_bc_angle_dependence == INVALID_BC_ANGLE_DEPENDENCE)
    {
      std::cerr << "Must specify angular dependence of Planckian BC on right edge\n";
      exit(EXIT_FAILURE);
    }
  }
  
  
  /// do not allow refelective condition and krylov solving
  if( (m_rad_bc_left == REFLECTIVE) 
    && ( (m_wg_solve_type == KRYLOV_SWEEPS) || (m_wg_solve_type == KRYLOV_DSA) ) )
  {
    std::cerr << "A refelective radiation boundary condition is not allowed with Krylov WGRS\n";
    exit(EXIT_FAILURE);  
  }
  
  return 0;
}