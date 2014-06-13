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
  
  int region_return = -1;
  int material_return = -1;
  int time_return = -1;
  int spatial_return = -1;
  int angle_return = -1;
  
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
  
  if( (region_return < 0) || (material_return < 0) || (time_return < 0) || 
      (spatial_return < 0) || (angle_return < 0) )
  {
    good_exit = false;
  }  
  
  return good_exit;
}

int Input_Reader::get_dfem_degree(void)
{
  return m_dfem_trial_space_degree;
}

QUADRATURE_TYPE Input_Reader::get_dfem_interpolation_point_type(void)
{
  return m_dfem_interpolation_point_type;
}

int Input_Reader::get_opacity_degree(void)
{
  return m_opacity_polynomial_degree;
}

OPACITY_TREATMENT Input_Reader::get_opacity_treatment(void)
{
  return m_opacity_treatment;
}

QUADRATURE_TYPE Input_Reader::get_opacity_interpolation_point_type(void)
{
  return m_opacity_interpolation_point_type;
}

MATRIX_INTEGRATION Input_Reader::get_integration_method(void)
{
  return m_integration_type;
}

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
    /// only tested if a region has logarithmic spacing
    TiXmlElement* space_factor = region_id->FirstChildElement("Log_space_factor");
    
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
    else if(spacing_str == "LOG_THIN_LEFT")
      m_region_spacing_type[i] = LOG_THIN_LEFT;
    else if(spacing_str == "LOG_THIN_RIGHT")
      m_region_spacing_type[i] = LOG_THIN_RIGHT;

    if(m_region_spacing_type[i] == INVALID_GRID_SPACING)
    {
      std::cerr << "Error.  Region " << i << " Invalid Grid Spacing" << std::endl;
      exit(EXIT_FAILURE);
    }
    
    /// load in scaling factor for logarithmic grid spacing and check that it is greter than 1
    if( (m_region_spacing_type[i] == LOG_THIN_LEFT ) || (m_region_spacing_type[i] == LOG_THIN_RIGHT))
    {
      if(!space_factor)
      {
        std::cerr << "Error.  Region " << i << " Missing size factor for logarithmic grid spacing" << std::endl;
        exit(EXIT_FAILURE);        
      }
      m_region_spacing_constant[i] = atof( space_factor->GetText() );
      if(m_region_spacing_constant[i] <= 1.)
      {
        std::cerr << "Error.  Region " << i << " Log spacing factors must be > 1" << std::endl;
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
  m_material_opacity_type.resize(m_number_materials,INVALID_OPACITY_TYPE);  
  m_material_cv_type.resize(m_number_materials,INVALID_CV_TYPE);
  m_material_source_type.resize(m_number_materials,INVALID_FIXED_SOURCE_TYPE);
  
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
    
    TiXmlElement* opacity_type = mat_descr->FirstChildElement( "Opacity_type");
    TiXmlElement* cv_type = mat_descr->FirstChildElement( "Cv_type");
    TiXmlElement* source_type = mat_descr->FirstChildElement( "Fixed_source_type");
    
    std::string opacity_str = opacity_type->GetText();
    std::string cv_str = cv_type->GetText();
    std::string source_str = source_type->GetText();
    
    transform(opacity_str.begin() , opacity_str.end() , opacity_str.begin() , toupper);
    transform(cv_str.begin() , cv_str.end() , cv_str.begin() , toupper);
    transform(source_str.begin() , source_str.end() , source_str.begin() , toupper);

    if(opacity_str == "CONSTANT_XS")
      m_material_opacity_type[mat_num] = CONSTANT_XS;
    else if(opacity_str == "RATIONAL")
      m_material_opacity_type[mat_num] = RATIONAL;
    else if(opacity_str == "TABLE_LOOKUP")
      m_material_opacity_type[mat_num] = TABLE_LOOKUP;
   
    if(m_material_opacity_type[mat_num] == INVALID_OPACITY_TYPE)
    {
      std::cerr << "Error.  Invalid opacity type for material " << mat_num << std::endl;
      exit(EXIT_FAILURE);    
    }
    
    if(source_str == "NO_SOURCE")
      m_material_source_type[mat_num] = NO_SOURCE;
   
    if(m_material_source_type[mat_num] == INVALID_FIXED_SOURCE_TYPE)
    {
      std::cerr << "Error.  Invalid source type for material " << mat_num << std::endl;
      exit(EXIT_FAILURE);    
    }
    
    if(cv_str == "CONSTANT_CV")
      m_material_cv_type[mat_num] = CONSTANT_CV;
   
    if(m_material_cv_type[mat_num] == INVALID_CV_TYPE)
    {
      std::cerr << "Error.  Invalid cv type for material " << mat_num << std::endl;
      exit(EXIT_FAILURE);    
    }
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
    m_time_start_meth = EXPONENTIAL;
  else if(starter_str == "VECTOR")
    m_time_start_meth = VECTOR;
    
  if(m_time_step_scheme == INVALID_TIME_SOLVER)
  {
    std::cerr << "Error.  Invalid time stepping scheme" << std::endl;
    exit(EXIT_FAILURE);  
  }
  if(m_time_start_meth == INVALID_STARTING_METHOD)
  {
    std::cerr << "Error.  Invalid time iteration starting scheme" << std::endl;
    exit(EXIT_FAILURE);  
  }
  
  if(m_time_start_meth == EXPONENTIAL)
  {
    /// Get increase factor
    TiXmlElement* exp_incr_fact = time_elem->FirstChildElement( "Increase_factor");
    if(!exp_incr_fact)
    {
      std::cerr << "Error.  Missing Increase_factor element.  Required for Starting_method EXPONENTIAL" << std::endl;
      exit(EXIT_FAILURE);
    }
    m_starting_constants.resize(1,0.);
    m_starting_constants[0] = atof( exp_incr_fact->GetText() );
    if(m_starting_constants[0] < 1.)
    {
      std::cerr << "Error.  Starting time step increase must be greater than 1.0" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  else if(m_time_start_meth == VECTOR)
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
    
    m_starting_constants.resize(2*num_vec_stages , 0.);
    
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
      m_starting_constants[2*curr_stage] = double (stages);
      m_starting_constants[2*curr_stage + 1] = size_factor;      
      
      vector_stage = vector_stage->NextSiblingElement( "Vector_stage");
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
  TiXmlElement* opacity_degree_elem = spatial_element->FirstChildElement( "Opacity_polynomial_degree");

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
  if(!opacity_degree_elem)
  {
    std::cerr << "Error.  Missing Opacity_degree_elem element from SPATIAL_DISCRETIZATION block" << std::endl;
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
  else if(opacity_treat_str == "INTERPOLATING")
    m_opacity_treatment = INTERPOLATING;
    
  if(m_opacity_treatment == INVALID_OPACITY_TREATMENT)
  {
    std::cerr << "Errror. Invalid opacity treatment" << std::endl;
    exit(EXIT_FAILURE);    
  }
  
  m_opacity_polynomial_degree = atoi( opacity_degree_elem->GetText() );
  if(m_opacity_polynomial_degree < 0)
  {
    std::cerr << "Error.  Invalid opacity polynomial" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  return 0;  
}

  int Input_Reader::load_angular_discretization_data(TiXmlElement* angle_element)
{
  TiXmlElement* n_angle_elem = angle_element->FirstChildElement( "Number_of_angles");
  TiXmlElement* n_group_elem = angle_element->FirstChildElement( "Number_of_groups");
  TiXmlElement* quad_type_elem = angle_element->FirstChildElement( "Quadrature_type");
  
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
    
  return 0;
}