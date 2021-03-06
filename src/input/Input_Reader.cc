/** @file   Input_Reader.cc
  *   @author pmaginot
  *   @brief Implement the Input_Reader class that reads XML files with TinyXml
 */
         
#include "Input_Reader.h"
// ##########################################################
// Public functions 
// ##########################################################

void Input_Reader::read_xml(std::string xmlFile)
{  
  TiXmlDocument doc( xmlFile.c_str() );
  
  bool loaded = doc.LoadFile();
  
  if( loaded  )
  {
    std::cout << "Found input file: " << xmlFile << std::endl;
  }
  else{
    std::cout << "Looking for file: " << xmlFile << std::endl;
    throw Dark_Arts_Exception( INPUT , "Error reading input file, file non-existent or XML error" );
  }
  
  /// save the name of the input file (without directory path)
  unsigned int found = xmlFile.find_last_of("/");
  m_short_input_file = xmlFile.substr(found+1);
  
  TiXmlElement* restart_elem = doc.FirstChildElement( "RESTART_FILE" );  

  if(restart_elem)
  {
    /// mesh refinement or restart run    
    load_restart_problem(restart_elem);
  }
  else
  {
    /// run beginning from scratch
    load_from_scratch_problem(doc);
  }
  return;
}

// ##########################################################
// Const (get) Public functions 
// ##########################################################
  
void Input_Reader::get_time_start_vectors(std::vector<double>& step_size_in_vector_stage, std::vector<int>& steps_in_vector_stage) const
{
  step_size_in_vector_stage = m_vector_start_sizes;
  steps_in_vector_stage = m_vector_start_step_numbers;
  return;
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

void Input_Reader::get_scattering_poly_coeff(const int mat_num , std::vector<double>& coeff) const
{
  coeff = m_scat_opacity_poly[mat_num];
  return;
}

void Input_Reader::get_absorption_poly_coeff(const int mat_num , std::vector<double>& coeff) const
{
  coeff = m_abs_opacity_poly[mat_num];
  return;
}

double Input_Reader::get_cv_constant(const int mat_num) const
{
  return m_cv_constants[mat_num];
} 

double Input_Reader::get_region_temperature(const int reg_num) const
{
  if(reg_num > m_number_regions)
    throw Dark_Arts_Exception( INPUT ,  "Asking for a region temperature in a region greater than n_regions" );
  
  return m_region_temperature[reg_num];
}

double Input_Reader::get_region_radiation_temperature(const int reg_num) const
{
  if(reg_num > m_number_regions)
    throw Dark_Arts_Exception( INPUT ,  "Asking for a region radiation temperature in a region greater than n_regions" );
    
  return m_region_radiation_temperature[reg_num];
}

RADIATION_IC_TYPE Input_Reader::get_radiation_ic_type(void) const
{
  return m_radiation_ic_type;
}
  

/* ***************************************************
 *
 *  Protected Functions
 *
 * ************************************************** */

 void Input_Reader::load_restart_problem(TiXmlElement* restart_elem)
{
  TiXmlElement* type_elem = restart_elem->FirstChildElement("Restart_type");
  if(!type_elem)
    throw Dark_Arts_Exception(INPUT, "Missing Restart_type element");
    
  std::string type_str = type_elem->GetText();
  transform(type_str.begin() , type_str.end() , type_str.begin() , toupper);
  
  if(type_str == "MESH_REFINEMENT")
  {
    m_is_mesh_refinement = true;
    m_restart_type = MESH_REFINEMENT;
  }
  else if(type_str == "RESTART")
  {
    m_is_restart = true;
    m_restart_type = RESTART;
  }
  
  TiXmlElement* data_dump_path = restart_elem->FirstChildElement( "Data_dump_path" );
  if(!data_dump_path)
    throw Dark_Arts_Exception(INPUT, "Every RESTART_FILE element must have Data_dump_path element ");
    
  m_data_dump_str = data_dump_path->GetText();
  
  TiXmlElement* initial_input = type_elem->FirstChildElement( "Initial_inputfile");
  if(!initial_input)
    throw Dark_Arts_Exception(INPUT, "Must have Initial_inputfile block for any restart type!");
  m_initial_input_str = initial_input->GetText();
  
  switch(m_restart_type)
  {
    case MESH_REFINEMENT:
    {
      /// load everything basically the same as before, but change the number of cells per region from the initial file
      TiXmlElement* refine_factor = type_elem->FirstChildElement( "Refinement_factor");
      TiXmlElement* time_refinement_factor = type_elem->FirstChildElement( "Time_refinement_factor" );
      TiXmlElement* trial_space_degree_elem = type_elem->FirstChildElement( "Trial_space_degree" );
      
      TiXmlDocument doc( m_initial_input_str.c_str() );  
      bool loaded = doc.LoadFile();  
      if( loaded  )
      {
        std::cout << "Found initial input file: " << m_initial_input_str << std::endl ;
      }
      else{
        std::cout << "Failed to find initial input file: " << m_initial_input_str << std::endl;
        throw Dark_Arts_Exception( INPUT , "Error reading initial input file, file non-existent or XML error" );
      }      
      
      /// get the initial input file
      load_from_scratch_problem(doc);
      
      /// change things about the initial input file
      /// number of spatial cells
      if(refine_factor)
      {        
        m_refinement_factor = atoi( refine_factor->GetText() );
        if(m_refinement_factor < 1)
          throw Dark_Arts_Exception(INPUT , "Spatial refinement factor must be an integer greater than 1");     
        
        for(int i =0 ; i < m_number_regions ; i++) 
          m_cells_per_region[i] *= m_refinement_factor;      
      }      
      
      /// time step size
      if(time_refinement_factor)
      {
        double time_factor = atof(time_refinement_factor->GetText() );
        if(time_factor < 1.)
          throw Dark_Arts_Exception(INPUT, "Time refinement factor must be greater than 1.");
          
        m_dt_min /= time_factor;
        m_dt_max /= time_factor;
      }
      
      /// trial space degree
      if(trial_space_degree_elem)
      {
        int trial_space_degree = atoi(trial_space_degree_elem->GetText() );
        if( trial_space_degree < 1)
          throw Dark_Arts_Exception(INPUT , "Invalid trial space degree for MESH_REFINEMENT");
        else
          m_dfem_trial_space_degree = trial_space_degree;        
      }
      
      if( (!trial_space_degree_elem) && (!time_refinement_factor) && (!refine_factor))
        throw Dark_Arts_Exception(INPUT, "MESH_REFINEMENT input without any mesh refinement parameters");
      
      break;
    }
    case RESTART:
    {
      /** what do we need to do for this?
        last time step completed, time in simulation (get new t_start)
        time that we want to end (new t_end)
        time stepping parameters, excluding time stepping scheme (MUST remain the same)
        
        path where data files that hold I, T , and phi,
        time step where we want to resume
        initial input file all other data will remain the same / can be populated with already existing routines        
      */
      TiXmlElement* step_element = type_elem->FirstChildElement( "Last_complete_time_step_number" );      
      
      if(!step_element)
        throw Dark_Arts_Exception(INPUT, "RESTART requires previous ending time step (integer)");      
        
      m_restart_step = atoi(step_element->GetText());
      if(m_restart_step < 1)
        throw Dark_Arts_Exception(INPUT, "Invalid previous time step number.  Must be greater than 0");
            
      TiXmlDocument doc( m_initial_input_str.c_str() );  
      bool loaded = doc.LoadFile();  
      if( loaded  )
      {
        std::cout << "Found initial input file: " << m_initial_input_str << std::endl ;
      }
      else{
        std::cout << "Failed to find initial input file: " << m_initial_input_str << std::endl;
        throw Dark_Arts_Exception( INPUT , "Error reading initial input file, file non-existent or XML error" );
      }      
      
      load_from_scratch_problem(doc);
      
      /// Time stepping data      
      TiXmlElement* new_time_element = type_elem->FirstChildElement( "New_time_stepping_data");      
      TIME_SOLVER old_time_step_scheme = m_time_step_scheme;
      load_time_stepping_data(new_time_element);
      if(old_time_step_scheme != m_time_step_scheme)
        throw Dark_Arts_Exception(INPUT, "Restart simulation must use the same time stepping scheme as old simulation");     
      
      break;
    }
    case INVALID_RESTART_TYPE:
    {
      throw Dark_Arts_Exception(INPUT, "RESTART_FILE inputs require a valid entry in Restart_type element");      
    }
  }
  
  return;
}

void Input_Reader::load_from_scratch_problem(TiXmlDocument& doc)
{
  TiXmlElement* inp_block = doc.FirstChildElement( "INPUT_FILE" );
  
  if( !inp_block )
    throw Dark_Arts_Exception( INPUT ,  "INPUT_FILE block not found. " );
  
  //! Every Input File Is Required To Have These Blocks
  TiXmlElement* reg_elem = inp_block->FirstChildElement( "REGIONS" );
  TiXmlElement* mat_elem = inp_block->FirstChildElement( "MATERIALS" );
  TiXmlElement* time_elem = inp_block->FirstChildElement( "TIME" );
  TiXmlElement* discr_elem = inp_block->FirstChildElement( "SPATIAL_DISCRETIZATION" );
  TiXmlElement* angle_elem = inp_block->FirstChildElement( "ANGULAR_DISCRETIZATION" );
  TiXmlElement* solver_elem = inp_block->FirstChildElement( "SOLVER" );
  TiXmlElement* bc_ic_elem = inp_block->FirstChildElement( "BC_IC" );  
  TiXmlElement* output_elem = inp_block->FirstChildElement( "OUTPUT");
  
  if(!reg_elem)
    throw Dark_Arts_Exception( INPUT , "REGIONS block not found.");
  
  if(!mat_elem)
    throw Dark_Arts_Exception( INPUT , "MATERIALS block not found.");
    
  if(!time_elem)
    throw Dark_Arts_Exception( INPUT ,  "TIME block not found. ");
  
  if(!discr_elem)
    throw Dark_Arts_Exception( INPUT , "SPATIAL_DISCRETIZATION block not found.");
 
  if(!angle_elem)
    throw Dark_Arts_Exception( INPUT , "ANGULAR_DISCRETIZATION block not found.");
  
  if(!solver_elem)
    throw Dark_Arts_Exception( INPUT ,  "SOLVER block not found. ");
 
  if(!bc_ic_elem)
    throw Dark_Arts_Exception( INPUT , "BC_IC block not found. ");  
    
  if(!output_elem)
    throw Dark_Arts_Exception(INPUT , "OUTPUT block not found");
     
  //! Load Data Appropriately from each input block
  load_region_data(reg_elem);  
  // std::cout << "Past reader block 1" << std::endl;
  load_angular_discretization_data(angle_elem);
  // std::cout << "Past reader block 2" << std::endl;
  /// load_material needs angular data to already be loaded
  load_material_data(mat_elem);
  // std::cout << "Past reader block 3" << std::endl;
  load_time_stepping_data(time_elem);
  // std::cout << "Past reader block 4" << std::endl;
  load_spatial_discretization_data(discr_elem);
  // std::cout << "Past reader block 5" << std::endl;
  load_solver_data(solver_elem);
  // std::cout << "Past reader block 6" << std::endl;
  load_bc_ic_data(bc_ic_elem);
  // std::cout << "Past reader block 7" << std::endl;
  load_output_data(output_elem);
  // std::cout << "Past reader block 8" << std::endl;
  
  
  return;
}  

int Input_Reader::load_output_data(TiXmlElement* output_element)
{
  /**
    Need to load the following data, output type, checkpointing frequency, filename base if 
  */
  TiXmlElement* type_elem = output_element->FirstChildElement( "Output_type" );
  TiXmlElement* check_point_elem = output_element->FirstChildElement( "Checkpoint_frequency" );
  TiXmlElement* output_dir_elem = output_element->FirstChildElement( "Output_directory" );
  if(!type_elem)
    throw Dark_Arts_Exception(INPUT, "Missing Output_type element");
    
  if(!check_point_elem)
    throw Dark_Arts_Exception(INPUT, "In OUTPUT element, missing Checkpoint_frequency element"); 
    
  if(!output_dir_elem)
    throw Dark_Arts_Exception(INPUT , "In OUTPUT block, missing Output_directory element");
    
  TiXmlElement* suppress_elem = output_element->FirstChildElement( "Suppress_data_dumps" );
  if( suppress_elem )
    m_suppress_output_dumps = true;
    
  m_output_directory = output_dir_elem->GetText() ;
    
  std::string type_str  = type_elem->GetText();
  transform(type_str.begin() , type_str.end() , type_str.begin() , toupper);
  if(type_str == "DUMP")
  {
    m_output_type = DUMP;
    
    TiXmlElement* multi_time_elem = type_elem->FirstChildElement("Multi_time");
    if(multi_time_elem)
    {
      m_n_data_outputs = atoi(multi_time_elem->GetText() );
      if(m_n_data_outputs < 2)
      {
        throw Dark_Arts_Exception(INPUT , "In Multi_time element number of dumps must be >= 2");
      }
      m_times_to_dump.resize(m_n_data_outputs,0.);
      TiXmlElement* dump_time_elem = multi_time_elem->FirstChildElement("Dump_time");
      for(int i = 0; i < m_n_data_outputs ; i++)
      {
        if(!dump_time_elem)
        {
          std::stringstream err;
          err << "Missing Dump_time: " << i << " in Multi_time dump";
          throw Dark_Arts_Exception(INPUT , err );
        }
        m_times_to_dump[i] = atof( dump_time_elem->GetText() );
        dump_time_elem = dump_time_elem->NextSiblingElement("Dump_time");
      }      
      
      for(int i = 0 ; i < (m_n_data_outputs-1) ; i++)
      {
        if( m_times_to_dump[i] >= m_times_to_dump[i+1] )
          throw Dark_Arts_Exception(INPUT , "Dump_time elements must be in ascending order");
      }
    }
        
  }
  else if(type_str =="SPACE_TIME_ERROR")
  {
    m_output_type = SPACE_TIME_ERROR;
  }
  else if(type_str=="END_SPACE_ERROR")
  {
    m_output_type = END_SPACE_ERROR;
  } 
  else if(type_str == "BOTH_ERRORS")
  {
    m_output_type = BOTH_ERRORS;
  }
  
  if(m_output_type == INVALID_OUTPUT_TYPE)
    throw Dark_Arts_Exception(INPUT , "Invalid output type");
    
  if( m_output_type == END_SPACE_ERROR) 
    m_end_space_error = true;
  else if(m_output_type == SPACE_TIME_ERROR)
    m_space_time_error = true;
  else if(m_output_type == BOTH_ERRORS)
  {
    m_end_space_error = true;
    m_space_time_error = true;
  }
  

  
  if(m_end_space_error || m_space_time_error)
  {
    TiXmlElement* filename_elem = type_elem->FirstChildElement( "Results_filename_base" );
    if(!filename_elem)
      throw Dark_Arts_Exception(INPUT , "Must have Results_filename_base for END_SPACE_ERROR or SPACE_TIME_ERROR or BOTH_ERRORS types");
      
    m_results_file_base = filename_elem->GetText();
    
    /// check to make sure this is a MMS problem, since this is the only type of problem we have comparator solutions for
    if(m_material_radiation_source_type[0] != MMS_SOURCE)
      throw Dark_Arts_Exception(INPUT, "Can only be an MMS problem to calculate spatial errors!");
    
  }  
    
  m_restart_frequency = atoi(check_point_elem->GetText() ) ;
  if(m_restart_frequency < 1)
    throw Dark_Arts_Exception(INPUT , "Restart frequency must be a positive integer >= 1");
  
  return 0;
}
 
int Input_Reader::load_region_data(TiXmlElement* region_element)
{  
  /// Get the number of regions in the problem
  TiXmlElement* num_regions_elem = region_element->FirstChildElement( "Number_of_regions" );
  if(num_regions_elem)
    m_number_regions = atoi(num_regions_elem->GetText());
  else
    throw Dark_Arts_Exception( INPUT , "In REGIONS element: Missing Number of regions element");
  
  TiXmlElement* num_materials_elem = region_element->FirstChildElement( "Number_of_materials");
  if(num_regions_elem)
    m_number_materials = atoi(num_materials_elem->GetText());
  else
    throw Dark_Arts_Exception( INPUT , "In REGIONS element:  Missing Number of materials element");
  
  m_cells_per_region.resize(m_number_regions,0);
  m_region_material_numbers.resize(m_number_regions,0);
  m_region_spacing_type.resize(m_number_regions,INVALID_GRID_SPACING);
  m_region_left_bounds.resize(m_number_regions,-1.0);
  m_region_right_bounds.resize(m_number_regions,-2.0);
  m_region_spacing_constant.resize(m_number_regions,0.0);
  m_region_min_size.resize(m_number_regions,0.0);
  
  TiXmlElement* region_id = region_element->FirstChildElement( "Region");
  if(!region_id)
    throw Dark_Arts_Exception( INPUT , "In REGIONS element: Missing Individual Region Elements");
    
  /// loop over the specified number of regions
  for(int i=0; i<m_number_regions ; i++)
  {
    /// Get all of the elements required for a valid region definition
    int reg_num = atoi(region_id->GetText());
    if(reg_num != i)
    {
      std::stringstream err;
      err << "In REGIONS element: Expected Region: " <<  i  << " Got Region: " <<  reg_num ;
      throw Dark_Arts_Exception( INPUT , err );
    }
      
    TiXmlElement* n_cells = region_id->FirstChildElement( "N_cells" );
    TiXmlElement* x_left = region_id->FirstChildElement("Left_bound");
    TiXmlElement* x_right = region_id->FirstChildElement("Right_bound");
    TiXmlElement* spacing = region_id->FirstChildElement("Spacing");
    TiXmlElement* mat_number = region_id->FirstChildElement("Material_number");
    
    if(!n_cells || !x_left || !x_right || !spacing || !mat_number)
    {
      std::stringstream err; 
      err << "In REGIONS block, Region: " << reg_num << " is missing required elements";
      throw Dark_Arts_Exception( INPUT , err );      
    }
    
    /// store all region specific data    
    /// check region number, demand that all regions be input in order from 0 ... n-1
    m_cells_per_region[i] = atoi(n_cells->GetText());
    if(m_cells_per_region[i] < 0)
    {
      std::stringstream err;
      err      << "In REGIONS element: Region " << i << " Number of Cells must be a positive integer";
      throw Dark_Arts_Exception( INPUT , err );
    }
    
    /// check material number, must be in range: 0 ... n_mat - 1
    m_region_material_numbers[i] = atoi(mat_number->GetText());
    if( (m_region_material_numbers[i] < 0) ||(m_region_material_numbers[i] > (m_number_materials-1)))
    {
      std::stringstream err;
      err      << "Error.  Region " << i << " Material number out of range" ;
      throw Dark_Arts_Exception( INPUT , err );
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
      std::stringstream err;
      err << "In REGIONS element: Region " << i << " Invalid Grid Spacing";
      throw Dark_Arts_Exception( INPUT , err );
    }
    
    /// load in scaling factor for logarithmic grid spacing and check that it is greter than 1
    if( (m_region_spacing_type[i] == LOG) )
    {
      TiXmlElement* space_factor = spacing->FirstChildElement("Log_space_factor");
      TiXmlElement* min_dx = spacing->FirstChildElement("Min_cell_size");
      if(!space_factor)
      {
        std::stringstream err;
        err << "In REGIONS element: Region " << i << " Missing size factor for logarithmic grid spacing" ;
        throw Dark_Arts_Exception( INPUT , err );       
      }
      if(!min_dx)
      {
        std::stringstream err;
        err << "In REGIONS element:  Region " << i << " Missing minimum cell size for logarithmic spacing" ;
        throw Dark_Arts_Exception( INPUT , err );        
      }
      
      m_region_spacing_constant[i] = atof( space_factor->GetText() );
      if(m_region_spacing_constant[i] < 0.)
      {
        std::stringstream err;
        err << "In REGIONS element:   Region " << i << " Log spacing factors must be > 0" ;
        throw Dark_Arts_Exception( INPUT , err );
      }
      m_region_min_size[i] = atof( min_dx->GetText() );
      if(m_region_min_size[i] < 0.)
      {
        std::stringstream err;
        err << "In REGIONS element: Region " << i << " Minimum cell spacing must be > 0" ;
        throw Dark_Arts_Exception( INPUT , err );
      }      
    }
    
    /// get left and right values of region, check that x_L < x_R
    m_region_left_bounds[i] = atof( x_left->GetText() ) ;
    m_region_right_bounds[i] = atof(x_right->GetText() ) ;
    if(m_region_left_bounds[i] > m_region_right_bounds[i] )
    {
      std::stringstream err;
      err << "In REGIONS element:  Region " << i << " x_L > x_R" ;
      throw Dark_Arts_Exception( INPUT , err );
    }
      
    /// go to the next region
    region_id = region_id->NextSiblingElement( "Region");
    if( (region_id == 0) && (i< m_number_regions - 1))
      throw Dark_Arts_Exception( INPUT , "In REGIONS block: Insufficient number of Region elements for the number of regions" );
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
  m_cv_rational_powers.resize(m_number_materials);
  m_cv_rational_offsets.resize(m_number_materials);
  
  m_scat_opacity_poly.resize(m_number_materials);
  m_abs_opacity_poly.resize(m_number_materials);
  
  m_cv_polynomial_coeff.resize(m_number_materials);
  m_cv_poly_max_power.resize(m_number_materials,0);
  
  m_rad_source_start_time.resize(m_number_materials,0.);
  m_rad_source_end_time.resize(m_number_materials,0.);
  m_rad_source_output.resize(m_number_materials,0.);
  
  m_temp_source_start_time.resize(m_number_materials,0.);
  m_temp_source_end_time.resize(m_number_materials,0.);
  m_temp_source_output.resize(m_number_materials,0.);
    
  TiXmlElement* units_elem = mat_elem->FirstChildElement("Units");
  if(!units_elem)
  {
    throw Dark_Arts_Exception( INPUT , "In MATERIALS block: Missing required Units block in MATERIALS.");
  }
  else
  {
    std::string units_str = units_elem->GetText();    
    transform(units_str.begin() , units_str.end() , units_str.begin() , toupper);
    
    if( units_str == "CM_SH_KEV")
    {
      m_units_type = CM_SH_KEV;
      std::cout << "Not using any funny units. \n";
      std::cout << "Lengths in cm \n";
      std::cout << "Time in shakes \n";
      std::cout << "Temperatures / Energies in keV\n"; 
      
    }
    else if( units_str == "UNITY" )
    {
      m_units_type = UNITY;
      std::cout << "Using unity / dimensionless numbers\n";
      std::cout << "a=c=1" << std::endl;
    }
    else
    {
      throw Dark_Arts_Exception( INPUT , "In MATERIALS block:Invalid units type.");
    }
  }
  
  TiXmlElement* mat_descr = mat_elem->FirstChildElement("Material");
  for(int mat_cnt = 0; mat_cnt < m_number_materials ; mat_cnt++)
  { 
    std::stringstream err;
    if(!mat_descr)
    {      
      err << "Missing Material element in MATERIALS block.  Expected: " << m_number_materials  << "found: " << mat_cnt;
      throw Dark_Arts_Exception( INPUT , err );
    }    
    int mat_num = atoi( mat_descr->GetText() );
    if(mat_num != mat_cnt)
      throw Dark_Arts_Exception( INPUT , "In MATERIALS block: Materials not entered in order" );
    
    TiXmlElement* scat_opacity_type = mat_descr->FirstChildElement( "Scattering_opacity_type");
    TiXmlElement* abs_opacity_type = mat_descr->FirstChildElement( "Absorption_opacity_type");
    TiXmlElement* cv_type = mat_descr->FirstChildElement( "Cv_type");
    TiXmlElement* rad_source_type = mat_descr->FirstChildElement( "Radiation_fixed_source_type");
    TiXmlElement* temp_source_type = mat_descr->FirstChildElement( "Temperature_fixed_source_type");
    if(!scat_opacity_type)
    {
      err      << "In MATERIALS block: Missing scattering opacity type in material " << mat_num;
      throw Dark_Arts_Exception( INPUT , err );
    }
    if(!abs_opacity_type)
    {
      err      << "In MATERIALS block: Missing absorption opacity type in material " << mat_num ;
      throw Dark_Arts_Exception( INPUT , err );
    }
    if(!cv_type)
    {
      err      << "In MATERIALS block: Missing cv type in material " << mat_num;
      throw Dark_Arts_Exception( INPUT , err );
    }
    if(!rad_source_type)
    {
      err      << "In MATERIALS block: Missing radiation source type in material " << mat_num ;
      throw Dark_Arts_Exception( INPUT , err );
    }
    if(!temp_source_type)
    {
      err      << "In MATERIALS block: Missing temperautre source type in material " << mat_num ;
      throw Dark_Arts_Exception( INPUT , err );
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
        err      << "In MATERIALS block:Missing constant_value tag for material: " << mat_num << " absorption opacity" ;
        throw Dark_Arts_Exception( INPUT , err );
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
        err  << "In MATERIALS block: Missing Multiplier tag for material " << mat_num << " RATIONAL absorption opacity";
        throw Dark_Arts_Exception( INPUT , err );
      }
      if(!denom_power_val)
      {
        err  << "In MATERIALS block: Missing Denominator_power tag for material " << mat_num << " RATIONAL absorption opacity";
        throw Dark_Arts_Exception( INPUT , err );
      }
      if(!denom_offset_val)
      {
        err  << "In MATERIALS block: Missing Denominator_offset tag for material " << mat_num << " RATIONAL absorption opacity";
        throw Dark_Arts_Exception( INPUT , err );
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
        err  << "In MATERIALS block: Missing File_name tag for material " << mat_num << " TABLE_LOOKUP absorption opacity " ;
        throw Dark_Arts_Exception( INPUT , err );
      }
      m_abs_opacity_str[mat_num] = abs_op_file->GetText();
    }
    else if(abs_opacity_str == "POLYNOMIAL_SPACE")
    {
      m_material_absorption_opacity_type[mat_num] = POLYNOMIAL_SPACE;
      TiXmlElement* abs_poly = abs_opacity_type->FirstChildElement( "Highest_polynomial_degree" );
      if(!abs_poly)
      {
        err  << "In MATERIALS block: Missing Highest_polynomial_degree tag for material " << mat_num << " POLYNOMIAL_SPACE absorption opacity " ;
        throw Dark_Arts_Exception( INPUT , err );
      }
      m_abs_opacity_integer_constants[mat_num] = atoi( abs_poly->GetText() );
      if(m_abs_opacity_integer_constants[mat_num] < 1)
      {
        err  << "In MATERIALS block: Highest_polynomial_degree tag for material " << mat_num << " POLYNOMIAL_SPACE absorption opacity less than 1" ;
        throw Dark_Arts_Exception( INPUT , err );
      }
      
      m_abs_opacity_poly[mat_num].resize(m_abs_opacity_integer_constants[mat_num] + 1,0.);
      
      TiXmlElement* poly_coeff = abs_opacity_type->FirstChildElement( "Degree_coefficient" );
      for(int p = 0 ; p < m_abs_opacity_integer_constants[mat_num] + 1 ; p++)
      {
        if(!poly_coeff)
        {
          err << "In MATERIALS block, missing polynomial coefficent for degree " << p << " term in material " << mat_num;
          throw Dark_Arts_Exception( INPUT , err );
        }
        
        if( atoi( poly_coeff->GetText() ) != p)
        {
          err << "In MATERIALS block, missing polynomial coefficents for material " << mat_num << "absorption opacity out of order";
          throw Dark_Arts_Exception( INPUT , err );
        }
        else
        {
          TiXmlElement* poly_val = poly_coeff->FirstChildElement( "Coefficient_value" );
          if(!poly_val)
          {
            err << "Missing value element in " << mat_num << " Polynomial aborption opactity term " << p << " Degree_coefficient block";
            throw Dark_Arts_Exception(INPUT, err );
          }
          else
          {
            m_abs_opacity_poly[mat_num][p] = atof( poly_val->GetText() );            
          }
        }
        
        poly_coeff = poly_coeff->NextSiblingElement( "Degree_coefficient" );
      }
    }
    
    if(m_material_absorption_opacity_type[mat_num] == INVALID_OPACITY_TYPE)
    {
      err  << "In MATERIALS block: Invalid absorption opacity type for material " << mat_num ;
      throw Dark_Arts_Exception( INPUT , err ); 
    }
    
    /// set-up / scan for scattering opacity data
    if(scat_opacity_str == "CONSTANT_XS")
    {
      m_material_scattering_opacity_type[mat_num] = CONSTANT_XS;
      TiXmlElement* const_val = scat_opacity_type->FirstChildElement( "Constant_value" );
      if(!const_val)
      {
        err  << "In MATERIALS block: Missing constant_value tag for material: " << mat_num << " scattering opacity" ;
        throw Dark_Arts_Exception( INPUT , err );
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
        err  << "In MATERIALS block: Missing Multiplier tag for material " << mat_num << " RATIONAL scattering opacity";
        throw Dark_Arts_Exception( INPUT , err );
      }
      if(!denom_power_val)
      {
        err  << "In MATERIALS block: Missing Denominator_power tag for material " << mat_num << " RATIONAL scattering opacity";
        throw Dark_Arts_Exception( INPUT , err );
      }
      if(!denom_offset_val)
      {
        err  << "In MATERIALS block: Missing Denominator_offset tag for material " << mat_num << " RATIONAL scattering opacity";
        throw Dark_Arts_Exception( INPUT , err );
      }
      m_scat_opacity_double_constants_1[mat_num] = atof(mult_val->GetText() );
      m_scat_opacity_double_constants_2[mat_num] = atof(denom_offset_val->GetText() );
      m_scat_opacity_integer_constants[mat_num] = atoi(denom_power_val->GetText() );
    }
    else if(scat_opacity_str == "TABLE_LOOKUP")
    {
      err  << "In MATERIALS block: Scattering opacities cannot be table look-up.  In material " << mat_num << " TABLE_LOOKUP scattering opacity " ;
      throw Dark_Arts_Exception( INPUT , err );      
    }
    else if(scat_opacity_str == "POLYNOMIAL_SPACE")
    {
      m_material_scattering_opacity_type[mat_num] = POLYNOMIAL_SPACE;
      TiXmlElement* scat_poly = scat_opacity_type->FirstChildElement( "Highest_polynomial_degree" );
      if(!scat_poly)
      {
        err  << "In MATERIALS block: Missing Highest_polynomial_degree tag for material " << mat_num << " POLYNOMIAL_SPACE scattering opacity " ;
        throw Dark_Arts_Exception( INPUT , err );
      }
      m_scat_opacity_integer_constants[mat_num] = atoi( scat_poly->GetText() );
      if(m_scat_opacity_integer_constants[mat_num] < 1)
      {
        err  << "In MATERIALS block: Highest_polynomial_degree tag for material " << mat_num << " POLYNOMIAL_SPACE scattering less than 1" ;
        throw Dark_Arts_Exception( INPUT , err );
      }
      
      m_scat_opacity_poly[mat_num].resize(m_scat_opacity_integer_constants[mat_num] + 1,0.);
      
      TiXmlElement* poly_coeff = scat_opacity_type->FirstChildElement( "Degree_coefficient" );
      for(int p = 0 ; p < m_scat_opacity_integer_constants[mat_num] + 1 ; p++)
      {
        if(!poly_coeff)
        {
          err << "In MATERIALS block, missing polynomial coefficent for degree " << p << " term in material " << mat_num << " scattering opacity";
          throw Dark_Arts_Exception( INPUT , err );
        }
        
        if( atoi( poly_coeff->GetText() ) != p)
        {
          err << "In MATERIALS block, missing polynomial coefficients for material " << mat_num << "scattering opacity out of order";
          throw Dark_Arts_Exception( INPUT , err );
        }
        else
        {
          TiXmlElement* poly_val = poly_coeff->FirstChildElement( "Coefficient_value" );
          if(!poly_val)
          {
            err << "Missing value element in " << mat_num << " Polynomial scattering opactity term " << p << " Degree_coefficient block";
            throw Dark_Arts_Exception(INPUT, err );
          }
          else
          {
            m_scat_opacity_poly[mat_num][p] = atof( poly_val->GetText() );            
          }
        }
        
        poly_coeff = poly_coeff->NextSiblingElement( "Degree_coefficient" );
      }
    }
    
    if(m_material_scattering_opacity_type[mat_num] == INVALID_OPACITY_TYPE)
    {
      err  << "In MATERIALS block:  Invalid absorption opacity type for material " << mat_num ;
      throw Dark_Arts_Exception( INPUT , err );   
    }
    
    if(rad_source_str == "NO_SOURCE")
    {
      m_material_radiation_source_type[mat_num] = NO_SOURCE;
    }
    else if(rad_source_str == "MMS_SOURCE")
    {
      m_material_radiation_source_type[mat_num] = MMS_SOURCE;
    }
    else if(rad_source_str == "CONSTANT_SOURCE")
    {
      m_material_radiation_source_type[mat_num] = CONSTANT_SOURCE;
    }
    
    /// check that we have a valid source type
    if(m_material_radiation_source_type[mat_num] == INVALID_FIXED_SOURCE_TYPE)
    {
      err  << "In MATERIALS block: Invalid radiation source type for material " << mat_num ;
      throw Dark_Arts_Exception( INPUT , err );   
    }
        
    if(m_material_radiation_source_type[mat_num] == MMS_SOURCE)
    {
      /// restrictions on MMS simulations
      if(m_number_materials > 1)
       throw Dark_Arts_Exception(INPUT, "MMS_SOURCE only allowed for single material problems");
       
      if( m_material_absorption_opacity_type[0] == TABLE_LOOKUP )
        throw Dark_Arts_Exception(INPUT, "MMS_SOURCE not allowed with TABLE_LOOKUP opacities");
      
      if( m_number_groups != 1) 
        throw Dark_Arts_Exception(INPUT, "MMS_SOURCE only valid for grey simulations");
        
      if(m_n_leg_moments > 1)
        throw Dark_Arts_Exception(INPUT, "MMS_SOURCE only valid for storing P0 angular moments");
    
      /// require a radiation space element, radiation angle element , temperature space element , and temporal dependence element
      TiXmlElement* rad_space_elem = rad_source_type->FirstChildElement( "Radiation_space");
      TiXmlElement* temp_space_elem = rad_source_type->FirstChildElement( "Temperature_space");
      TiXmlElement* mms_time_elem = rad_source_type->FirstChildElement( "Temporal");
      TiXmlElement* rad_angle_elem = rad_source_type->FirstChildElement( "Radiation_angle");
      
      if(!rad_space_elem || !temp_space_elem || ! mms_time_elem || !rad_angle_elem)
        throw Dark_Arts_Exception(INPUT, "MMS_SOURCE requires Radiation_space, Temperature_space, Temporal, and Radiation_angle elements");
        
      std::string rad_space_str = rad_space_elem->GetText();
      std::string temp_space_str = temp_space_elem->GetText();
      std::string mms_time_str = mms_time_elem->GetText();
      std::string rad_angle_str = rad_angle_elem->GetText();
      
      transform(rad_space_str.begin() , rad_space_str.end() , rad_space_str.begin() , toupper);
      transform(temp_space_str.begin() , temp_space_str.end() , temp_space_str.begin() , toupper);
      transform(mms_time_str.begin() , mms_time_str.end() , mms_time_str.begin() , toupper);
      transform(rad_angle_str.begin() , rad_angle_str.end() , rad_angle_str.begin() , toupper);
      
      
      if(rad_space_str == "RAD_POLY_SPACE")
      {
        m_mms_rad_space = RAD_POLY_SPACE;
      }
      else if(rad_space_str == "RAD_COS_SPACE" )
      {
        m_mms_rad_space = RAD_COS_SPACE;
      }
      
      if(temp_space_str == "TEMP_POLY_SPACE")
      {
        m_mms_temp_space = TEMP_POLY_SPACE;
      }
      else if(temp_space_str == "TEMP_COS_SPACE" )
      {
        m_mms_temp_space = TEMP_COS_SPACE;
      }
      
      if(mms_time_str == "POLY_TIME")
      {
        m_mms_time = POLY_TIME;
      }
      else if(mms_time_str == "COS_TIME")
      {
        m_mms_time = COS_TIME;
      }
      
      if(rad_angle_str == "MMS_ISOTROPIC")
      {
        m_mms_rad_angle = MMS_ISOTROPIC;
      }
      else if( rad_angle_str == "MMS_ANGLE_POLY")
      {
        m_mms_rad_angle = MMS_ANGLE_POLY;
      }
      
      switch(m_mms_rad_space)
      {
        case RAD_POLY_SPACE:
        {
          load_mms_poly_constants( rad_space_elem , m_rad_space_mms_const);
          break;
        }
        case RAD_COS_SPACE:
        {
          load_mms_cos_constants( rad_space_elem, m_rad_space_mms_const );
          break;
        }
        case INVALID_RADIATION_SPACE_MMS:
        {
          throw Dark_Arts_Exception(INPUT, "Invalid Radiation spatial dependence for MMS problems");
          break;
        }
      }
      
      switch(m_mms_temp_space)
      {
        case TEMP_POLY_SPACE:
        {
          load_mms_poly_constants( temp_space_elem , m_temp_space_mms_const);
          break;
        }
        case TEMP_COS_SPACE:
        {
          load_mms_cos_constants( temp_space_elem, m_temp_space_mms_const );
          break;
        }
        case INVALID_TEMPERATURE_SPACE_MMS:
        {
          throw Dark_Arts_Exception(INPUT, "Invalid temperature spatial dependence for MMS problems");
          break;
        }
      }
      
      switch(m_mms_time)
      {
        case COS_TIME:
        {
          load_mms_cos_constants( mms_time_elem , m_time_mms_const);
          break;
        }
        case POLY_TIME:
        {
          load_mms_poly_constants( mms_time_elem , m_time_mms_const);
          break;
        }
        case INVALID_TIME_MMS_TYPE:
        {
          throw Dark_Arts_Exception(INPUT, "Invalid mms time dependence for MMS problems");
          break;
        }
      }
      
      switch(m_mms_rad_angle)
      {
        case MMS_ISOTROPIC:
        {
          break;
        }
        case MMS_ANGLE_POLY:
        {
          load_mms_poly_constants( rad_angle_elem , m_mms_angle_coeff);
          break;
        }
        case INVALID_RADIATION_ANGLE_MMS:
        {
          throw Dark_Arts_Exception(INPUT, "Invalid MMS angle dependence");
          break;
        }
      }
      
      
    }
    else if(m_material_radiation_source_type[mat_num] == CONSTANT_SOURCE)
    {
      TiXmlElement* t_start = rad_source_type->FirstChildElement("Time_start");
      TiXmlElement* t_end = rad_source_type->FirstChildElement("Time_end");
      TiXmlElement* rad_output = rad_source_type->FirstChildElement("Isotropic_output");
      TiXmlElement* temperature_output = rad_source_type->FirstChildElement( "Temperature_source" ); 
      
      if(!t_start)
      {
        err << "In Material: " << mat_num << " radiation source type CONSTANT_SOURCE requires Time_start element";
        throw Dark_Arts_Exception(INPUT , err);
      }
      
      if(!t_end)
      {
        err << "In Material: " << mat_num << " radiation source type CONSTANT_SOURCE requires Time_end element";
        throw Dark_Arts_Exception(INPUT , err);
      }
      
      if( (!rad_output) && (!temperature_output) )
      {
        err << "In Material: " << mat_num << " radiation source type CONSTANT_SOURCE requires Isotropic_output element or Temperature_source";
        throw Dark_Arts_Exception(INPUT , err);
      }
      
      if( rad_output && temperature_output )
      {
        err << "In Material: " << mat_num << " radiation source type CONSTANT_SOURCE requires Isotropic_output element or Temperature_source, not both";
        throw Dark_Arts_Exception(INPUT , err);
      }
      
      std::cout << "Checked for all tags" << std::endl;
      
      m_rad_source_start_time[mat_num] = atof( t_start->GetText() );
      m_rad_source_end_time[mat_num] = atof( t_end->GetText() );
            
      if(m_rad_source_start_time[mat_num] > m_rad_source_end_time[mat_num])
      {
        err << "Material: " << mat_num << " constant rad source start time must be less than end time";
        throw Dark_Arts_Exception(INPUT , err);
      }    
      
      if( rad_output)
      {
        m_rad_source_output[mat_num] = atof( rad_output->GetText() );
        if(m_rad_source_output[mat_num] < 0.)
        {
          err << "Material: " << mat_num << " constant rad source Isotropic_output must be positive";
          throw Dark_Arts_Exception(INPUT , err);
        }      
      }
      
      if( temperature_output)
      {
        m_rad_source_output[mat_num] = atof( temperature_output->GetText() ) ;
        if(m_rad_source_output[mat_num] < 0.)
        {
          err << "Material: " << mat_num << " constant rad source Temperature_source must be positive";
          throw Dark_Arts_Exception(INPUT , err);
        }      
      }
      
    }

    if(temp_source_str == "NO_SOURCE")
    {
      m_material_temperature_source_type[mat_num] = NO_SOURCE;
    }
    else if(temp_source_str == "MMS_SOURCE")
    {
      m_material_temperature_source_type[mat_num] = MMS_SOURCE;
      if( m_material_radiation_source_type[mat_num] != MMS_SOURCE )
        throw Dark_Arts_Exception(INPUT, "For MMS_Source must be defined for temperature and radiation fixed source");
    }
    else if(temp_source_str == "CONSTANT_SOURCE")
    {
      m_material_temperature_source_type[mat_num] = CONSTANT_SOURCE;
      TiXmlElement* t_start = temp_source_type->FirstChildElement( "Time_start");
      TiXmlElement* t_end = temp_source_type->FirstChildElement( "Time_end");
      TiXmlElement* temp_output = temp_source_type->FirstChildElement( "Source_output");
      
      if(!t_start)
      {
        err << "In Material: " << mat_num << " temperature source type CONSTANT_SOURCE requires Time_start element";
        throw Dark_Arts_Exception(INPUT , err);
      }
      
      if(!t_end)
      {
        err << "In Material: " << mat_num << " temperature source type CONSTANT_SOURCE requires Time_end element";
        throw Dark_Arts_Exception(INPUT , err);
      }
      
      if(!temp_output)
      {
        err << "In Material: " << mat_num << " temperature source type CONSTANT_SOURCE requires Source_output element";
        throw Dark_Arts_Exception(INPUT , err);
      }
      
      m_temp_source_start_time[mat_num] = atof(t_start->GetText() );
      m_temp_source_end_time[mat_num] = atof(t_end->GetText() );
      m_temp_source_output[mat_num] = atof(temp_output->GetText() );
      
      if(m_temp_source_start_time[mat_num] >= m_temp_source_end_time[mat_num])
      {
        err << "Material: " << mat_num << " constant temperature source start time must be less than end time";
        throw Dark_Arts_Exception(INPUT , err);
      }
      if(m_temp_source_output[mat_num] < 0.)
      {
        err << "Material: " << mat_num << " constant temperature source must be positive";
        throw Dark_Arts_Exception(INPUT , err);
      }      
    }    
    
    if(m_material_radiation_source_type[mat_num] == MMS_SOURCE)
    {
      if( m_material_temperature_source_type[mat_num] != MMS_SOURCE)      
        throw Dark_Arts_Exception(INPUT, "Must specify MMS_SOURCE for temperature and radiation or neither");
    }
    
    if(m_material_temperature_source_type[mat_num] == INVALID_FIXED_SOURCE_TYPE)
    {
      err  << "In MATERIALS block: Invalid temperature source type for material " << mat_num;
      throw Dark_Arts_Exception( INPUT , err );   
    }
       
    if(cv_str == "CONSTANT_CV")
    {
      m_material_cv_type[mat_num] = CONSTANT_CV;
      TiXmlElement* cv_const = cv_type->FirstChildElement( "Cv_constant" );
      if(!cv_const)
      {
        err  << "In MATERIALS block:   Missing Cv_constant tag in material " << mat_num;
        throw Dark_Arts_Exception( INPUT , err );
      }
      m_cv_constants[mat_num] = atof(cv_const->GetText() );
      if(m_cv_constants[mat_num] < 0. )
      {
        err  << "In MATERIALS block:  Invalid Cv in material "<< mat_num << " values must be positive floats";
        throw Dark_Arts_Exception( INPUT , err );
      }
    }
    else if(cv_str == "RATIONAL_CV")
    {
      m_material_cv_type[mat_num] = RATIONAL_CV;
      TiXmlElement* cv_const_rat = cv_type->FirstChildElement( "Rational_cv_constant" );
      TiXmlElement* cv_offset = cv_type->FirstChildElement( "Rational_cv_offset" );
      TiXmlElement* cv_power = cv_type->FirstChildElement( "Rational_cv_power" );
      
      if(!cv_const_rat || !cv_offset || !cv_power)
      {
        err << "In material: " << mat_num << " Rational_Cv must have Rational_cv_constant, Rational_cv_offset, and Rational_cv_power elements";
        throw Dark_Arts_Exception(INPUT, err );
      }
      
      m_cv_constants[mat_num] = atof(cv_const_rat->GetText() );
      m_cv_rational_powers[mat_num] = atoi(cv_power->GetText() );
      m_cv_rational_offsets[mat_num] = atof( cv_offset->GetText() );
      
      if( m_cv_constants[mat_num]  < 0.)
      {
        err << "In material: " << mat_num << " Rational_Cv must have positive Rational_cv_constant";
        throw Dark_Arts_Exception(INPUT, err );
      }
      
      if( m_cv_rational_powers[mat_num]  < 1)
      {
        err << "In material: " << mat_num << " Rational_Cv must have positive integer Rational_cv_power";
        throw Dark_Arts_Exception(INPUT, err );
      }
      
      if( m_cv_rational_offsets[mat_num]  < 0.)
      {
        err << "In material: " << mat_num << " Rational_Cv must have positive Rational_cv_offset";
        throw Dark_Arts_Exception(INPUT, err );
      }      
    }
    else if( cv_str == "POLYNOMIAL_CV")
    {
      m_material_cv_type[mat_num] = POLYNOMIAL_CV;
      TiXmlElement* cv_poly_degree_elem = cv_type->FirstChildElement( "Polynomial_cv_degree" );
      TiXmlElement* cv_coeff = cv_type->FirstChildElement( "Polynomial_cv_coefficient" );
      m_cv_poly_max_power[mat_num] = atoi(cv_poly_degree_elem->GetText() );
      if( m_cv_poly_max_power[mat_num] < 0)
      {
        err << "In material: " << mat_num << " polynomial Cv requires degree >=0 ";
        throw Dark_Arts_Exception(INPUT, err );      
      }
      m_cv_polynomial_coeff[mat_num].resize(m_cv_poly_max_power[mat_num] + 1, 0.);
      for( int i = 0 ; i <= m_cv_poly_max_power[mat_num]  ; i++)
      {
        if(!cv_coeff)
        {
          err << "Missing coefficent of polynomial degree: " << i << " in material: " << mat_num<< " Cv_polynomial"; 
          throw Dark_Arts_Exception(INPUT, err ) ; 
        }
        if(i != atoi(cv_coeff->GetText() ) )
        {
          err << "Coefficents of Cv polynomial do not match input labeling for degree: " << i << " in material: " << mat_num<< " Cv_polynomial"; 
          throw Dark_Arts_Exception(INPUT, err ) ; 
        }
             
        TiXmlElement* val_elem = cv_coeff->FirstChildElement( "Value" );
        if(!val_elem)
        {
          err << "Coefficient: " << i << " in material " << mat_num << " Cv_Polynomial is missinng Value element";
          throw Dark_Arts_Exception(INPUT, err );
        }
        m_cv_polynomial_coeff[mat_num][i] = atof(val_elem->GetText() );
        
        cv_coeff = cv_coeff->NextSiblingElement("Polynomial_cv_coefficient");
      }      
    }
    
    if(m_material_cv_type[mat_num] == INVALID_CV_TYPE)
    {
      err  << "In MATERIALS block:  Invalid cv type for material " << mat_num ;
      throw Dark_Arts_Exception( INPUT , err );  
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
  TiXmlElement* adaptive_time_step_elem = time_elem->FirstChildElement( "Adaptive_time_step" );
  
  if(!dt_min_elem)
    throw Dark_Arts_Exception( INPUT , "In TIME block: Missing Dt_min element in TIME input block" );
    
  if(!dt_max_elem)
    throw Dark_Arts_Exception( INPUT , "In TIME block: Missing Dt_max element in TIME input block" );
    
  if(!t_start_elem)
    throw Dark_Arts_Exception( INPUT , "In TIME block:   Missing T_start element in TIME input block");
    
  if(!t_end_elem)
    throw Dark_Arts_Exception( INPUT , "In TIME block:  Missing T_end element in TIME input block" );
    
  if(!solver_elem)
    throw Dark_Arts_Exception( INPUT , "In TIME block:  Missing Time_solver element in TIME input block" );
    
  if(!start_meth_elem)
    throw Dark_Arts_Exception( INPUT , "In TIME block:  Missing Starting_method element in TIME input block" );
    
   if(!adaptive_time_step_elem)
   {
    m_adaptive_time_method = NO_ADAPTIVE_CONTROL;
   }
   else
   {
    m_adaptive_time_method = INVALID_ADAPTIVE_TIME_STEP_CONTROL;
    std::string adaptive_str = adaptive_time_step_elem->GetText() ;
    transform(adaptive_str.begin() , adaptive_str.end() , adaptive_str.begin() , toupper);
    if( adaptive_str == "NO_ADAPTIVE_CONTROL" )
    {
      m_adaptive_time_method = NO_ADAPTIVE_CONTROL;
    }
    else if( adaptive_str == "DUMB_POINTWISE_T" )
    {
      m_adaptive_time_method = DUMB_POINTWISE_T;
    }
    else if( adaptive_str == "CHANGE_IN_T" )
    {
      m_adaptive_time_method = CHANGE_IN_T;
    }
    else if( adaptive_str == "CHANGE_IN_T_VOLUMETRIC")
    {
      m_adaptive_time_method = CHANGE_IN_T_VOLUMETRIC;
      TiXmlElement* cells_per_grouping_elem = adaptive_time_step_elem->FirstChildElement( "Cells_per_grouping" );
      if(!cells_per_grouping_elem)
        throw Dark_Arts_Exception(INPUT , "CHANGE_IN_T_VOLUMETRIC requires Cells_per_grouping element");
        
      m_cells_per_adaptive_group = atoi( cells_per_grouping_elem->GetText() );
      if(m_cells_per_adaptive_group < 1 )
        throw Dark_Arts_Exception(INPUT , "Cells_per_grouping must be an integer greater than 0 for CHANGE_IN_T_VOLUMETRIC scheme"); 
    }
    if( (m_adaptive_time_method == CHANGE_IN_T_VOLUMETRIC) || (m_adaptive_time_method == CHANGE_IN_T) || (m_adaptive_time_method == DUMB_POINTWISE_T) )
    {    
      TiXmlElement* offset_elem = adaptive_time_step_elem->FirstChildElement( "Offset_temperature");
      if(!offset_elem)
      {
        m_adaptive_temperature_floor = 0.;
      }
      else
      {
        m_adaptive_temperature_floor = atof( offset_elem->GetText() );
        if(m_adaptive_temperature_floor < 0.)
          throw Dark_Arts_Exception(INPUT , "T_Change_Volumetric floor element requires positive temperature");
      }
      
      TiXmlElement* target_change_elem = adaptive_time_step_elem->FirstChildElement("Target_increase_fraction") ;
      if(!target_change_elem)
        throw Dark_Arts_Exception(INPUT , "CHANGE_IN_T, CHANGE_IN_T_VOLUMETRIC, and DUMB_POITNWISE_T require Target_increase_fraction tag");
        
      m_change_in_t_goal = atof( target_change_elem->GetText() );
      std::cout << "Target increase fraction: " << m_change_in_t_goal << std::endl;
      if( (m_change_in_t_goal < 0.) || ( m_change_in_t_goal > 1.) )
        throw Dark_Arts_Exception(INPUT , "Target_increase_fraction for adaptive time step scheme must be a decimal between 0 and 1");    
    }
    
    if( m_adaptive_time_method == INVALID_ADAPTIVE_TIME_STEP_CONTROL)
      throw Dark_Arts_Exception(INPUT, "Valid Adaptive_time_step control not given");
   }
  
  m_dt_min = atof(dt_min_elem->GetText() );
  m_dt_max = atof(dt_max_elem->GetText() );
  m_t_start = atof(t_start_elem->GetText() );
  m_t_end = atof(t_end_elem->GetText() );
  
  if( (m_dt_min < 0.) || (m_dt_max < 0.) )
    throw Dark_Arts_Exception( INPUT , "In TIME block:  Time step sizes must be positive floats" );
    
  if( m_dt_min > m_dt_max )
    throw Dark_Arts_Exception( INPUT , "In TIME block:  Minimum time step must be greater than maximum time step" );
    
  if(m_t_start > m_t_end)
    throw Dark_Arts_Exception( INPUT , "In TIME block:  t_end must be greater than t_start ");
  
  std::string stepper_str = solver_elem->GetText();
  std::string starter_str = start_meth_elem->GetText();
  transform(stepper_str.begin() , stepper_str.end() , stepper_str.begin() , toupper);
  transform(starter_str.begin() , starter_str.end() , starter_str.begin() , toupper);

  if(stepper_str == "IMPLICIT_EULER")
  {
    m_time_step_scheme = IMPLICIT_EULER;
  }
  else if(stepper_str == "ALEXANDER_2_2")
  {
    m_time_step_scheme = ALEXANDER_2_2;
  }  
  else if(stepper_str == "ALEXANDER_2_2_PLUS")
  {
    m_time_step_scheme = ALEXANDER_2_2_PLUS;
  }  
  else if(stepper_str == "ALEXANDER_3_3")
  {
    m_time_step_scheme = ALEXANDER_3_3;
  }
  
  if(starter_str == "EXPONENTIAL")
  {
    m_time_starting_method = EXPONENTIAL;
  }
  else if(starter_str == "VECTOR")
  {
    m_time_starting_method = VECTOR;
  }
  else if(starter_str == "RAMP")
  {
    m_time_starting_method = RAMP;
   }
   else if(starter_str == "ADAPTIVE")
  {   
      m_time_starting_method = ADAPTIVE;
  }
  
  if(m_time_step_scheme == INVALID_TIME_SOLVER)
    throw Dark_Arts_Exception( INPUT , "In TIME block:  Invalid time stepping scheme" );
  
  switch (m_time_starting_method)
  {
    case INVALID_STARTING_METHOD:
    {
      throw Dark_Arts_Exception( INPUT , "In TIME block:  Invalid time iteration starting scheme" );      
      break;
    }
    case EXPONENTIAL:
    {
      /// Get increase factor
      TiXmlElement* exp_incr_fact = start_meth_elem->FirstChildElement( "Increase_factor");
      if(!exp_incr_fact)
        throw Dark_Arts_Exception( INPUT , "In TIME block: Missing Increase_factor element.  Required for Starting_method EXPONENTIAL" );
     
      m_exponential_ratio = atof( exp_incr_fact->GetText() );
      if(m_exponential_ratio < 1.)
        throw Dark_Arts_Exception( INPUT , "In TIME block: Starting time step increase must be greater than 1.0");
    
      break;
    }
    case RAMP:
    {
      TiXmlElement* ramp_steps = start_meth_elem->FirstChildElement( "Ramp_steps" );
      if(!ramp_steps)
        throw Dark_Arts_Exception( INPUT , "In TIME block:  RAMP time starter requires Ramp_steps element" );
      
      m_ramp_steps = atoi(ramp_steps->GetText() );      
      if(m_ramp_steps < 1)
        throw Dark_Arts_Exception( INPUT , "In TIME block: RAMP starter requires at least one time step before full time step can be taken");
        
      break;
    }
    case VECTOR:
    {
      TiXmlElement* n_stages = start_meth_elem->FirstChildElement( "N_vector_stages" );
      if(!n_stages)
        throw Dark_Arts_Exception( INPUT , "In TIME block:   VECTOR time starter requires N_stages element" );
        
      m_num_vec_stages = atoi(n_stages->GetText() );
      if(m_num_vec_stages < 1)
        throw Dark_Arts_Exception( INPUT , "In TIME block:   Must have at least 1 stage of small time steps with VECTOR time starter" );
      
      m_vector_start_sizes.resize(m_num_vec_stages , 0.);
      m_vector_start_step_numbers.resize(m_num_vec_stages , 0);
      
      TiXmlElement* vector_stage = start_meth_elem->FirstChildElement("Vector_stage");
      for(int i=0; i<m_num_vec_stages ; i++)
      {
        if(!vector_stage)
        {
          std::stringstream err;
          err << "In TIME block: Expected " << m_num_vec_stages << " Vector_stage elements" << std::endl;
          err << "Only found: " << i << " Vector_stage elements" ;
          throw Dark_Arts_Exception( INPUT , err ); 
        }
        
        int curr_stage = atoi(vector_stage->GetText() );
        if(curr_stage < 0 || curr_stage == m_num_vec_stages)
          throw Dark_Arts_Exception( INPUT , "In TIME block: Vector_stage number must be between 0 and n_stages - 1" );
      
        TiXmlElement* stage_steps = vector_stage->FirstChildElement("Stage_steps");
        TiXmlElement* stage_factor = vector_stage->FirstChildElement("Stage_divisor");
        if(!stage_steps || !stage_factor)
          throw Dark_Arts_Exception( INPUT , "In TIME block: Expect a Stage_steps and Stage_divisor element in each Vector_stage element" );
        
        int stages = atoi( stage_steps->GetText() );
        double divisor = atof( stage_factor->GetText() );        
        if(stages < 1)
          throw Dark_Arts_Exception( INPUT , "In TIME block:  Stage_steps must be an integer greater than 1");
          
        if(divisor < 1. )
          throw Dark_Arts_Exception( INPUT , "In TIME block: Stage_divisor must be a double greater than 1. " );
          
        m_vector_start_step_numbers[curr_stage] = stages;
        m_vector_start_sizes[curr_stage] = divisor;      
        
        vector_stage = vector_stage->NextSiblingElement( "Vector_stage");
      }
      break;
    }
    case ADAPTIVE:
    {
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
  

  if(!trial_space_elem)
    throw Dark_Arts_Exception(INPUT, "SPATIAL_DISCRETIZATION Block: Missing DFEM_degree element");
    
  if(!int_type_elem)
    throw Dark_Arts_Exception(INPUT, "SPATIAL_DISCRETIZATION Block: Missing Integration_type element");
  
  if(!dfem_interp_point_type_elem)
    throw Dark_Arts_Exception(INPUT, "SPATIAL_DISCRETIZATION Block:  Missing DFEM_interpolation_point_type element");
    
  if(!opacity_treatment_elem)
    throw Dark_Arts_Exception(INPUT,  "Missing Opacity_treatment element from SPATIAL_DISCRETIZATION block" );
        
  m_dfem_trial_space_degree = atoi( trial_space_elem->GetText() );
  if(m_dfem_trial_space_degree < 1)
    throw Dark_Arts_Exception(INPUT, "SPATIAL_DISCRETIZATION Block: Trial space degree must be a positive integrer greater than 0" );
  
  std::string int_type_str = int_type_elem->GetText();
  transform(int_type_str.begin() , int_type_str.end() , int_type_str.begin() , toupper);
  if(int_type_str == "SELF_LUMPING")
    m_integration_type = SELF_LUMPING;
  else if(int_type_str == "TRAD_LUMPING")
    m_integration_type = TRAD_LUMPING;
  else if(int_type_str == "EXACT")
    m_integration_type = EXACT;
    
  if( m_integration_type == INVALID_MATRIX_INTEGRATION)
    throw Dark_Arts_Exception(INPUT, "SPATIAL_DISCRETIZATION Block:Invalid matrix integration strategy" );
  
  std::string dfem_intep_type_str = dfem_interp_point_type_elem->GetText();
  transform(dfem_intep_type_str.begin() , dfem_intep_type_str.end() , dfem_intep_type_str.begin() , toupper);
  if(dfem_intep_type_str == "GAUSS")
    m_dfem_interpolation_point_type = GAUSS;
  else if(dfem_intep_type_str == "LOBATTO")
    m_dfem_interpolation_point_type = LOBATTO;
  else if(dfem_intep_type_str == "EQUAL_SPACED")
    m_dfem_interpolation_point_type = EQUAL_SPACED;
    
  if( m_dfem_interpolation_point_type == INVALID_QUADRATURE_TYPE)
    throw Dark_Arts_Exception(INPUT, "SPATIAL_DISCRETIZATION Block:Invalid DFEM interpolation point type" );
    
  if( m_dfem_interpolation_point_type == INVALID_QUADRATURE_TYPE)
    throw Dark_Arts_Exception(INPUT, "SPATIAL_DISCRETIZATION Block: Invalid opacity interpolation point type" );
  
  std::string opacity_treat_str = opacity_treatment_elem->GetText();
  transform(opacity_treat_str.begin() , opacity_treat_str.end() , opacity_treat_str.begin() , toupper);
  if(opacity_treat_str == "MOMENT_PRESERVING")
    m_opacity_treatment = MOMENT_PRESERVING;
  else if(opacity_treat_str == "SLXS" )
    m_opacity_treatment = SLXS;
  else if(opacity_treat_str == "INTERPOLATING")
    m_opacity_treatment = INTERPOLATING;
    
  if(m_opacity_treatment == INVALID_OPACITY_TREATMENT)
    throw Dark_Arts_Exception(INPUT, "SPATIAL_DISCRETIZATION Block: Invalid opacity treatment" );
    
  if( (m_opacity_treatment == SLXS) && (m_integration_type != SELF_LUMPING) )    
    throw Dark_Arts_Exception(INPUT, "SPATIAL_DISCRETIZATION Block: SLXS opacity treatment can only be used with SELF_LUMPING integration" );
    
  
  /// if not using self lumping treatment, need to get degree of opacity spatial treatment
  if(m_opacity_treatment != SLXS)
  {   
    if(m_opacity_treatment == INTERPOLATING)
    {
      TiXmlElement* opacity_interp_point_type_elem = opacity_treatment_elem->FirstChildElement( "Opacity_interpolation_point_type");
      if(!opacity_interp_point_type_elem)
        throw Dark_Arts_Exception(INPUT, "SPATIAL_DISCRETIZATION Block: Missing Opacity_interpolation_point_type element ");
        
      std::string opacity_intep_type_str = opacity_interp_point_type_elem->GetText();
      transform(opacity_intep_type_str.begin() , opacity_intep_type_str.end() , opacity_intep_type_str.begin() , toupper);
      if(opacity_intep_type_str == "GAUSS")
        m_opacity_interpolation_point_type = GAUSS;
      else if(opacity_intep_type_str == "LOBATTO")
        m_opacity_interpolation_point_type = LOBATTO;
      else if(opacity_intep_type_str == "EQUAL_SPACED")
        m_opacity_interpolation_point_type = EQUAL_SPACED;
    }
    else if(m_opacity_treatment == MOMENT_PRESERVING)
      m_opacity_interpolation_point_type = GAUSS;
    
    TiXmlElement* opacity_degree_elem = opacity_treatment_elem->FirstChildElement( "Opacity_polynomial_degree");
    if(!opacity_degree_elem)
      throw Dark_Arts_Exception(INPUT, "SPATIAL_DISCRETIZATION Block:  Missing Opacity_degree_elem element");
    
    m_opacity_polynomial_degree = atoi( opacity_degree_elem->GetText() );
    if(m_opacity_polynomial_degree < 0)
      throw Dark_Arts_Exception(INPUT, "SPATIAL_DISCRETIZATION Block: Invalid opacity polynomial" );
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
    throw Dark_Arts_Exception( INPUT, "In ANGULAR_DISCRETIZATION block:  Missing Number_of_angles element" );
  
  if(!n_group_elem)
    throw Dark_Arts_Exception( INPUT, "In ANGULAR_DISCRETIZATION block:Missing Number_of_groups element");
  
  if(!quad_type_elem)
    throw Dark_Arts_Exception( INPUT, "In ANGULAR_DISCRETIZATION block: Missing Quadrature_type element");
  
  if(!n_legendre_mom_elem)
    throw Dark_Arts_Exception( INPUT, "In ANGULAR_DISCRETIZATION block:Missing Number_of_legendre_moments element");
  
  m_number_angles = atoi( n_angle_elem->GetText() );
  if( m_number_angles < 0 || ((m_number_angles % 2) != 0) )
    throw Dark_Arts_Exception( INPUT, "In ANGULAR_DISCRETIZATION block: Number of Angles must be a positive, even integer" );
  
  m_number_groups = atoi(n_group_elem->GetText() ) ;
  if( m_number_groups < 1)
    throw Dark_Arts_Exception( INPUT, "In ANGULAR_DISCRETIZATION block: Number of groups must be an integer > 0");
  
  m_group_lower_bounds.resize(m_number_groups,0.);
  m_group_upper_bounds.resize(m_number_groups,0.);
  
  /// load group bounds
  if( m_number_groups > 1)
  {
    TiXmlElement* grp_bounds = angle_element->FirstChildElement("Group_boundaries");
    for(int edge_cnt = 0; edge_cnt <= m_number_groups ; edge_cnt++)
    {    
      if(!grp_bounds)
      {
        std::stringstream err;
        err << "In ANGULAR_DISCRETIZATION block: Missing Group_Boundaries element " << m_number_groups + 1 
            << " found: " << edge_cnt ;
        throw Dark_Arts_Exception( INPUT, err );
      }
      int edge_num = atoi( grp_bounds->GetText() );
      if(edge_num != edge_cnt)
      throw Dark_Arts_Exception( INPUT, "In ANGULAR_DISCRETIZATION block:  Group Bounds not entered in order");
      
      TiXmlElement* grp_edge_val = grp_bounds->FirstChildElement("Edge_value");
      double edge = atof( grp_edge_val->GetText() );
      
      /// the lower the group number, the higher the average frequency group energy (by convention)
      if(edge_cnt < m_number_groups)
        m_group_upper_bounds[edge_cnt] = edge;
        
      if(edge_cnt > 0)
        m_group_lower_bounds[edge_cnt - 1] = edge;      
    
      grp_bounds = grp_bounds->NextSiblingElement("Group_boundaries");
    }
    
    /// check that these are logical values
    for(int g = 1; g<m_number_groups ; g++)
    {
      if( m_group_lower_bounds[g] > m_group_lower_bounds[g-1] )
      {
        std::stringstream err;
        err << "Frequency groups lower energy bounds are not in descending order\n" ;
        err << "Problem lies between group: " << g << " and group: " << g-1;
        throw Dark_Arts_Exception(INPUT, err );
      }
      if( m_group_upper_bounds[g] >  m_group_upper_bounds[g-1] )
      {
        std::stringstream err;
        err << "Frequency group upper energy bounds are not in descending order\n";
        err << "Problem lies between group: " << g << " and group: " << g-1;
        throw Dark_Arts_Exception(INPUT, err );
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
    throw Dark_Arts_Exception(INPUT, "In ANGULAR_DISCRETIZATION block: Invalid Quadrature_type ");
  
  m_n_leg_moments = atoi( n_legendre_mom_elem->GetText() );
  if(m_n_leg_moments < 1)
    throw Dark_Arts_Exception(INPUT, "In ANGULAR_DISCRETIZATION block:   Number of legendre moments must be >= 1");
    
  return 0;
}

int Input_Reader::load_solver_data(TiXmlElement* solver_element)
{
  TiXmlElement* solver_type_elem = solver_element->FirstChildElement( "WG_solver_type");
  TiXmlElement* wg_tolerance_elem = solver_element->FirstChildElement( "WG_tolerance");  
  TiXmlElement* damping_factor_elem = solver_element->FirstChildElement("Damping_factor");
  TiXmlElement* increase_elem = solver_element->FirstChildElement("Iteration_increase_factor");
  TiXmlElement* iter_before_damp_elem = solver_element->FirstChildElement("Iterations_before_damping");
  TiXmlElement* n_damps_elem = solver_element->FirstChildElement("Max_damping_resets");
  TiXmlElement* n_thermals_elem = solver_element->FirstChildElement("Max_thermal_iterations_per_stage");
  
  TiXmlElement* phi_norm_type_elem = solver_element->FirstChildElement("Convergence_norm_type");
  
  if(!solver_type_elem)
    throw Dark_Arts_Exception(INPUT, "In SOLVER element: Missing WG_Solver_type element") ;
  
  if(!wg_tolerance_elem)
    throw Dark_Arts_Exception(INPUT, "In SOLVER element: Missing WG_Tolerance element" );
  
  if(!damping_factor_elem)
    throw Dark_Arts_Exception(INPUT, "In solver element: missing Damping_factor element");
    
  if(!increase_elem)
    throw Dark_Arts_Exception(INPUT, "In solver element: missing Iteration_increase_factor element");
    
  if(!iter_before_damp_elem)
    throw Dark_Arts_Exception(INPUT, "In solver element: missing Iterations_before_damping element");
       
  if(!n_damps_elem)
    throw Dark_Arts_Exception(INPUT, "In SOLVER element, must give Max_damping_resets element");
    
  if(!n_thermals_elem)
    throw Dark_Arts_Exception(INPUT, "In SOLVER element, must provide a Max_thermal_iterations_per_stage element");
    
  /// default to pointwise error norm
  if(!phi_norm_type_elem)
  {
    m_err_norm_type = POINTWISE;
  }
  else
  {
    std::string norm_str = phi_norm_type_elem->GetText();
    transform(norm_str.begin() , norm_str.end() , norm_str.begin() , toupper);
    if(norm_str == "L1")
    {
      m_err_norm_type = L1;
    }
    else if( norm_str =="L1_RHO" )
    {
      m_err_norm_type = L1_RHO;
    }
    else if( norm_str == "POINTWISE" )
    {
      m_err_norm_type = POINTWISE;
    }
    else
    {
      throw Dark_Arts_Exception(INPUT, "Convergence_norm_type element found, but invalid convergence norm type given");
    }
  }
      
  m_max_damps = atoi(n_damps_elem->GetText() );
  m_iters_before_damp = atoi(iter_before_damp_elem->GetText() );
  m_damping_factor = atof(damping_factor_elem->GetText() );
  m_iter_increase_factor = atof(increase_elem->GetText() );
  m_max_thermals = atoi(n_thermals_elem->GetText() );
  
  if(m_max_thermals < 1)
    throw Dark_Arts_Exception(INPUT, "Max_thermal_iteration_per_stage must be >=1");
  
  if(m_max_damps < 0)
    throw Dark_Arts_Exception(INPUT, "Number of damps must be greater than 0");
      
  if( m_iters_before_damp < 2)
    throw Dark_Arts_Exception(INPUT, "Must have at least 2 thermal iterations before damping");
    
  if( (m_damping_factor < 0.) || ( fabs(m_damping_factor) > 1.) )
    throw Dark_Arts_Exception(INPUT, "Damping factor must be a positive float less than 1.");
    
  if( m_iter_increase_factor < 1.)
    throw Dark_Arts_Exception(INPUT, "Must increase iterations by at least a factor of 2 for damping");
  
    
  m_wg_tolerance = atof( wg_tolerance_elem->GetText() );
  if( (m_wg_tolerance < 1.E-15) || (m_wg_tolerance> 1.E-4))
    throw Dark_Arts_Exception(INPUT, "In SOLVER element: Invalid within group tolerance.  Must be greater 1E-15 and less than 1E-4");
    
  std::string solve_type_str = solver_type_elem->GetText();
  // change to all upper case 
  transform(solve_type_str.begin() , solve_type_str.end() , solve_type_str.begin() , toupper);
  if(solve_type_str == "FP_SWEEPS")
    m_wg_solve_type = FP_SWEEPS;
  else if (solve_type_str == "FP_DSA")
    m_wg_solve_type = FP_DSA;
  else if (solve_type_str == "KRYLOV_SWEEPS")
    m_wg_solve_type = KRYLOV_SWEEPS;
  else if (solve_type_str == "KRYLOV_DSA")
    m_wg_solve_type = KRYLOV_DSA;
  
  if(m_wg_solve_type == INVALID_WG_SOLVE_TYPE)
    throw Dark_Arts_Exception(INPUT, "In SOLVER element:Invalid within group radiation solver type");
  

  if( (m_wg_solve_type == FP_SWEEPS) || (m_wg_solve_type == FP_DSA))
  {
    TiXmlElement* num_sweep_elem = solver_type_elem->FirstChildElement( "Max_within_group_sweeps");
    if(!num_sweep_elem)
      throw Dark_Arts_Exception(INPUT, "In SOLVER element: Missing Max_Within_Group_Sweeps element." );
      
    m_max_num_sweeps = atoi( num_sweep_elem->GetText() );
    if(m_max_num_sweeps < 1)
      throw Dark_Arts_Exception(INPUT, "In SOLVER element: Must allow at least one sweep per within group solve");
  }
  
  if( (m_wg_solve_type == FP_DSA) || (m_wg_solve_type==KRYLOV_DSA) )
  {
    int total_cells = 0;
    if(m_number_regions < 0)
    {
      throw Dark_Arts_Exception(INPUT , "Cant perform this check, need to know the number of cells first");
    }
    else
    {
      for(int i = 0 ; i < m_number_regions ; i++)
        total_cells += m_cells_per_region[i];
    }
    if( total_cells < 2)
    {
      std::cout << "Stop being CUTE.  Need more cells for DSA operator to not die.  No DSA for you" << std::endl;
      
      if(m_wg_solve_type == FP_DSA)
      {
        m_wg_solve_type = FP_SWEEPS;
      }
        
      if(m_wg_solve_type == KRYLOV_DSA)
      {
        m_wg_solve_type=KRYLOV_SWEEPS;
      }
    }
    
    TiXmlElement* mip_options_elem = solver_type_elem->FirstChildElement( "MIP_solve_options" );
    if(!mip_options_elem)
      throw Dark_Arts_Exception(INPUT, "FP_DSA or KRYLOV_DSA require a MIP_solve_options block in WG_solver_type element");
      
    TiXmlElement* z_mip_elem = mip_options_elem->FirstChildElement( "Z_mip" );
    if(!z_mip_elem)
      throw Dark_Arts_Exception(INPUT , "Missing Z_mip element from WG MIP_solve_options_block");
      
    m_wg_z_mip = atof( z_mip_elem->GetText() );
    
    if(m_wg_z_mip < 2.)
      throw Dark_Arts_Exception(INPUT , "Z_mip must be a double greater than 1");
  }
  
  /**
      Things that are only required / neeeded for multi frequency problems
  */
  if(m_number_groups > 1)
  {
    TiXmlElement* bg_tolerance_elem = solver_element->FirstChildElement( "BG_tolerance");
    
    /// only require the between group solver tolerance if number of groups is greater than 1
    if(!bg_tolerance_elem)
      throw Dark_Arts_Exception(INPUT, "In SOLVER element: MF problems require specification of between group tolerance, BG_Tolerance.  Element not found");
        
    m_bg_tolerance = atof(bg_tolerance_elem->GetText() );
    if( m_bg_tolerance < m_wg_tolerance)
      throw Dark_Arts_Exception(INPUT, "In SOLVER element: Invalid between group tolerance.  Must be greater than the within group scattering tolerance.");
    
    if( m_bg_tolerance > 1.E-4)
      throw Dark_Arts_Exception(INPUT, "In SOLVER element: Invalid between group tolerance.  Must be less than 1.0E-4");
        
    TiXmlElement* mf_ard_solve_elem = solver_element->FirstChildElement( "MF_solver_type" );
    if(!mf_ard_solve_elem)
      throw Dark_Arts_Exception(INPUT, "In SOLVER element:MF_Solver_Type element not found in multi-frequency problem");
    
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
      throw Dark_Arts_Exception(INPUT, "In SOLVER element:Invalid ARD_SOLVE_TYPE found");
    
    if( (m_ard_solve_type == FP_NO_ACCEL) || (m_ard_solve_type == FP_LMFGA) )
    {
      TiXmlElement* fp_ard_iter_elem = mf_ard_solve_elem->FirstChildElement("Max_iterations");
      
      if(!fp_ard_iter_elem)
        throw Dark_Arts_Exception(INPUT, "In SOLVER element:Fixed point ARD solvers require Max_Iterations element");
      
      m_max_ard_iterations = atoi( fp_ard_iter_elem->GetText() );
      if(m_max_ard_iterations < 1)
        throw Dark_Arts_Exception(INPUT, "In SOLVER element:Require at least 1 iteration for FP ARD solver schemes");
    }
    
    if( (m_ard_solve_type == FP_LMFGA) || (m_ard_solve_type == KRYLOV_LMFGA) )
    {
      /// need to see if we want to collapse groups or not
      TiXmlElement* lmfga_type_elem = mf_ard_solve_elem->FirstChildElement("LMFGA_structure");
      if(!lmfga_type_elem)
        throw Dark_Arts_Exception(INPUT, "LMFGA_structure element missing from MF_solver_type element");
        
      std::string lmfga_structure_str = lmfga_type_elem->GetText();
      transform(lmfga_structure_str.begin() , lmfga_structure_str.end() , lmfga_structure_str.begin() , toupper);
      if(lmfga_structure_str == "NO_COLLAPSE")
      {
        m_lmfga_structure = NO_COLLAPSE;
        
        TiXmlElement* lmfga_order_elem = lmfga_type_elem->FirstChildElement("LMFGA_ordering");
        if(!lmfga_order_elem)
          throw Dark_Arts_Exception(INPUT, "LMFGA No_collapse requires LMFGA_ordering element in LMFGA_structure element");
          
        std::string lmfga_ordering_str = lmfga_order_elem->GetText();
        transform(lmfga_ordering_str.begin() , lmfga_ordering_str.end() , lmfga_ordering_str.begin() , toupper);
        
        if( lmfga_ordering_str == "GROUP_OUTER")
        {
          m_lmfga_ordering = GROUP_OUTER;
        }
        else if(lmfga_ordering_str == "CELL_OUTER")
        {
          m_lmfga_ordering = CELL_OUTER;
        }
        if(m_lmfga_ordering ==INVALID_LMFGA_ORDERING)
          throw Dark_Arts_Exception(INPUT, "Invalid LMFGA no collapse ordering");
      }
      else if(lmfga_structure_str == "GROUP_COLLAPSE")
      {
        m_lmfga_structure = GROUP_COLLAPSE;
      }
      
      if(m_lmfga_structure == INVALID_LMFGA_STRUCTURE)
        throw Dark_Arts_Exception(INPUT, "Invalid LMFGA_structure element value");
        
    }
  }
  
  /// get thermal convergence tolerance
  TiXmlElement* thermal_tolerance_elem = solver_element->FirstChildElement( "Thermal_tolerance");
  if(!thermal_tolerance_elem)
    throw Dark_Arts_Exception(INPUT, "In SOLVER element:Thermal_Tolerance element required");
   
  m_thermal_tolerance = atof( thermal_tolerance_elem->GetText() );
  if(m_thermal_tolerance < 1.E-14)
    throw Dark_Arts_Exception(INPUT, "In SOLVER element:Too small of thermal tolerance.") ; 
    
  if( m_number_groups > 1)
  {
    if( m_bg_tolerance < m_wg_tolerance)
      throw Dark_Arts_Exception(INPUT, "In SOLVER element: Between group tolerance must be greater within group tolerance ");
      
    if( m_thermal_tolerance < m_bg_tolerance)
      throw Dark_Arts_Exception(INPUT, "In SOLVER element: Thermal tolerance must be greater than between group absorption/re-emission tolerance ");
  }
  else
  {
    if(m_thermal_tolerance < m_wg_tolerance)
      throw Dark_Arts_Exception(INPUT, "In SOLVER element:Grey problem.  Within group tolerance must be smaller than thermal iteration tolerance");

  }
    
  return 0;
}

int Input_Reader::load_bc_ic_data(TiXmlElement* bc_ic_element)
{
  /// stream that we may use.  No need to declare everywhere
  std::stringstream err;

  TiXmlElement* temp_ic_type_elem = bc_ic_element->FirstChildElement( "Temperature_ic_type");
  TiXmlElement* rad_ic_type_elem = bc_ic_element->FirstChildElement( "Radiation_ic_type");
  
  TiXmlElement* rad_left_bc_type_elem = bc_ic_element->FirstChildElement( "Left_radiation_bc_type");
  TiXmlElement* rad_right_bc_type_elem = bc_ic_element->FirstChildElement( "Right_radiation_bc_type");
  
  if(!temp_ic_type_elem)
    throw Dark_Arts_Exception(INPUT, "In BC_IC element: Missing Temperature_IC_Type element" );
   
  if(!rad_ic_type_elem)
    throw Dark_Arts_Exception(INPUT, "In BC_IC element: Missing Radiation_IC_Type element" );
  
  
  if(!rad_left_bc_type_elem)
    throw Dark_Arts_Exception(INPUT, "In BC_IC element: Missing Left_Radiation_BC_Type element") ;
  
  if(!rad_right_bc_type_elem)
    throw Dark_Arts_Exception(INPUT, "In BC_IC element: Missing _Right_Radiation_BC_Type element") ;
  
  /// load BC data
  load_bc_data(m_left_bc, rad_left_bc_type_elem);
  load_bc_data(m_right_bc, rad_right_bc_type_elem);
  
  /// Get IC strings and set IC types
  std::string temp_ic_str = temp_ic_type_elem->GetText();
  transform(temp_ic_str.begin() , temp_ic_str.end() , temp_ic_str.begin() , toupper);
  std::string rad_ic_str = rad_ic_type_elem->GetText();
  transform(rad_ic_str.begin() , rad_ic_str.end() , rad_ic_str.begin() , toupper);
  
  if(temp_ic_str == "CONSTANT_TEMPERATURE_IC")
  {
    m_temperature_ic_type = CONSTANT_TEMPERATURE_IC;
  }
  else if(temp_ic_str == "MMS_TEMPERATURE_IC")
  { 
    m_temperature_ic_type = MMS_TEMPERATURE_IC;
  }
  
  if(m_material_temperature_source_type[0] == MMS_SOURCE)
  {
    if( m_temperature_ic_type != MMS_TEMPERATURE_IC)
    {
      err << "Must specify MMS_TEMPERATURE_IC in " << bc_ic_element->Value() << " Since tempature souce is MMS_SOURCE";
      throw Dark_Arts_Exception(INPUT, err );
    }
  }
    
  if(rad_ic_str == "PLANCKIAN_IC")
  {
    m_radiation_ic_type = PLANCKIAN_IC;
  }  
  else if(rad_ic_str == "MMS_RADIATION_IC")
  {
    m_radiation_ic_type = MMS_RADIATION_IC;
  }
  
  if(m_material_radiation_source_type[0] == MMS_SOURCE)
  {
    if( m_radiation_ic_type != MMS_RADIATION_IC)
    {
      err << "Must specify MMS_RADIATION_IC in " << bc_ic_element->Value() << " Since radiation souce is MMS_SOURCE";
      throw Dark_Arts_Exception(INPUT, err);
    }
  }
    
  /// get radiation initial conditions
  switch( m_radiation_ic_type)
  {
    case PLANCKIAN_IC:
    {
      m_region_radiation_temperature.resize(m_number_regions,0.);
      TiXmlElement* rad_ic_reg_elem = rad_ic_type_elem->FirstChildElement("Region");
      for(int reg = 0; reg < m_number_regions; reg++)
      {
        if(!rad_ic_reg_elem)
        {
          err << "In BC_IC element: Missing a region radiation temperature for every region.  Found: " << reg << " Need: " << m_number_regions ;
          throw Dark_Arts_Exception(INPUT,  err );
        }
        
        if( atoi( rad_ic_reg_elem->GetText() ) != reg)
        {
          err << "Expecting Region " << reg << " Radiation Temperature Found: " << atoi( rad_ic_reg_elem->GetText() ) ;
          throw Dark_Arts_Exception(INPUT,  err );
        }
        
        TiXmlElement* rad_temp_value = rad_ic_reg_elem->FirstChildElement("Radiation_temperature");
        if(!rad_temp_value)
        {
          err << "In BC_IC block: Missing required Radiation_Temperature element for region: " << reg ;
          throw Dark_Arts_Exception(INPUT,  err );
        }
        
        m_region_radiation_temperature[reg] = atof( rad_temp_value->GetText() );
        
        if(m_region_radiation_temperature[reg] < 0.)
        {
          err << "Invalid radiation temperature in region: " << reg << " , must be >= 0. ";
          throw Dark_Arts_Exception(INPUT,  err );
        }   
        
        rad_ic_reg_elem = rad_ic_reg_elem->NextSiblingElement("Region");        
      }
      break;
    }
    case MMS_RADIATION_IC:
    {
      
      break;
    }
    case INVALID_RADIATION_IC_TYPE:
    {
      throw Dark_Arts_Exception(INPUT, "In BC_IC element: INVALID_RADIATION_IC_TYPE");
      break;
    }
  }
  
  /// Get temperature initial conditions
  switch(m_temperature_ic_type)
  {
    case INVALID_TEMPERATURE_IC_TYPE:
    {
      throw Dark_Arts_Exception(INPUT,   "In BC_IC block: INVALID_TEMPERATURE_IC_TYPE");
      break;
    }
    case CONSTANT_TEMPERATURE_IC :
    {
      m_region_temperature.resize(m_number_regions,0.);
      TiXmlElement* temp_ic_reg_elem = temp_ic_type_elem->FirstChildElement( "Region");
      for(int reg = 0; reg < m_number_regions; reg++)
      {
        if(!temp_ic_reg_elem)
        {
          err << "Missing a region material temperature for every region.  Found: " << reg << " Need: " << m_number_regions ;
          throw Dark_Arts_Exception(INPUT,  err );
        }
        if( atoi( temp_ic_reg_elem->GetText() ) != reg)
        {
          err << "Expecting Region " << reg << " Found: " << atoi( temp_ic_reg_elem->GetText() ) ;
          throw Dark_Arts_Exception(INPUT,  err );
        }
        
        TiXmlElement* temp_value = temp_ic_reg_elem->FirstChildElement("Material_temperature");
        if(!temp_value)
        {
          err << "Missing required Material_Temperature element for region: " << reg ;
          throw Dark_Arts_Exception(INPUT,  err );
        }
        m_region_temperature[reg] = atof( temp_value->GetText() );
        if(m_region_temperature[reg] < 0.)
        {
          err << "Invalid material temperature in region: " << reg << " , must be >= 0. ";
          throw Dark_Arts_Exception(INPUT,  err );
        }
        
        temp_ic_reg_elem = temp_ic_reg_elem->NextSiblingElement("Region");
      }
      break;
    }  
    case MMS_TEMPERATURE_IC:
    {
      break;
    }
  }
  
  return 0;
}

void Input_Reader::load_mms_poly_constants(TiXmlElement* mms_element, std::vector<double>& poly_constants)
{
  std::stringstream err;
  TiXmlElement* poly_degree_elem = mms_element->FirstChildElement( "MMS_polynomial_degree" );
  if(! poly_degree_elem)
  {
    err << "Missing MMS_polynomial_degree element in " << mms_element->Value() ;
    throw Dark_Arts_Exception(INPUT, err );
  }
  
  int poly_degree = -1;
  poly_degree = atoi(poly_degree_elem->GetText() );
  if(poly_degree < 0 )
  {
    err << "Invalid MMS polynomial degree in " << mms_element->Value() ;
    throw Dark_Arts_Exception( INPUT, err );  
  }
  
  poly_constants.resize(poly_degree + 1, 0.);
  TiXmlElement* poly_val_elem = mms_element->FirstChildElement( "MMS_poly_coefficient" );
  for(int p=0 ; p <= poly_degree ; p++)
  {
    if(!poly_val_elem)
    {
      err << "Missing MMS_poly_coefficent element: " << p << " for MMS element: " << mms_element->Value() ;
      throw Dark_Arts_Exception(INPUT, err );
    }
    
    int degree = atoi(poly_val_elem->GetText() );
    if(degree != p )
    {
      err << "MMS_poly_coefficient out of order or missing in " << mms_element->Value();
      throw Dark_Arts_Exception(INPUT, err );
    }
    
    TiXmlElement* poly_coeff = poly_val_elem->FirstChildElement( "Coefficient" );
    if(!poly_coeff)
    {
      std::stringstream err;
      err << "Missing Coefficient element for " << poly_val_elem->Value() << " in " << mms_element->Value();
      throw Dark_Arts_Exception(INPUT, err );
    }
    poly_constants[p] = atof( poly_coeff->GetText() );
    
    poly_val_elem = poly_val_elem->NextSiblingElement( "MMS_poly_coefficient" );    
  }  
  
  return ;
}

void Input_Reader::load_mms_cos_constants(TiXmlElement* mms_element, std::vector<double>& cos_constants)
{
  TiXmlElement* a_elem = mms_element->FirstChildElement( "MMS_cos_a" );
  TiXmlElement* b_elem = mms_element->FirstChildElement( "MMS_cos_b" );
  TiXmlElement* c_elem = mms_element->FirstChildElement( "MMS_cos_c" );
  TiXmlElement* d_elem = mms_element->FirstChildElement( "MMS_cos_d" );
  
  std::stringstream err;
  if(!a_elem)
    err << "Missing MMS_cos_a element in " << mms_element->Value() ;
    
  if(!b_elem)
    err << "Missing MMS_cos_b element in " << mms_element->Value() ;
    
  if(!c_elem)
    err << "Missing MMS_cos_c element in " << mms_element->Value() ;
  
  if(!d_elem)
    err << "Missing MMS_cos_d element in " << mms_element->Value() ;
  
  if(err.tellp() > 0)
    throw Dark_Arts_Exception(INPUT, err );
  
  cos_constants.resize(4,0.);
  cos_constants[0] = atof( a_elem->GetText() );
  cos_constants[1] = atof( b_elem->GetText() );
  cos_constants[2] = atof( c_elem->GetText() );
  cos_constants[3] = atof( d_elem->GetText() );
  
  return;
}

void Input_Reader::load_bc_data(Radiation_BC_Data& bc_data, TiXmlElement* bc_type_elem)
{
  /// declare this line once rather than at each throw() statement
  std::stringstream err;
  
  /// Get BC string and set BC type
  std::string bc_str = bc_type_elem->GetText();
  transform(bc_str.begin() , bc_str.end() , bc_str.begin() , toupper);
  
  if(bc_str == "VACUUM_BC")
  {
    bc_data.type = VACUUM_BC;
  }
  else if(bc_str == "INCIDENT_BC")
  {
    bc_data.type = INCIDENT_BC;
  }
  else if(bc_str == "REFLECTIVE_BC")
  {
    bc_data.type = REFLECTIVE_BC;
  }
  else if(bc_str == "MMS_BC")
  {
    bc_data.type = MMS_BC;
  }  
  
  /// Make sure if we're specifying MMS problem that we have MMS boundaries (special case of incident flux since we know the solution)
  if( (m_material_radiation_source_type[0] == MMS_SOURCE) && ( bc_data.type != MMS_BC) )
  {
    err <<  "Must specify MMS_BC in " << bc_type_elem->Value() << " Since radiation souce is MMS_SOURCE";
    throw Dark_Arts_Exception(INPUT, err );
  }
  
  /// Handle left radiation boundary condition
  switch( bc_data.type )
  {
    case INVALID_RADIATION_BC_TYPE:
    {
      err << "In: " << bc_type_elem->Value() << " INVALID_RADIATION_BC_TYPE";
      throw Dark_Arts_Exception(INPUT, err );
      break;
    }
    case VACUUM_BC:
    {
      bc_data.value = 0.; 
      break;
    }
    case INCIDENT_BC:
    {
      /// get all required elements
      ///  TiXmlElement* rad_left_bc_energy_elem = rad_left_bc_type_elem->FirstChildElement( "Incident_energy");
      TiXmlElement* bc_angle_incidence_elem = bc_type_elem->FirstChildElement( "BC_angle_dependence");
      TiXmlElement* bc_time_dependence_elem = bc_type_elem->FirstChildElement( "BC_time_dependence");
      TiXmlElement* bc_value_type_elem = bc_type_elem->FirstChildElement("BC_value_type");
      TiXmlElement* bc_value_elem = bc_type_elem->FirstChildElement("BC_value");
      
      /// Check that we have all required elements
      if(!bc_value_type_elem)
      {
        err << "Missing BC_value_type element in " << bc_type_elem->Value() ;
        throw Dark_Arts_Exception(INPUT, err );
      }
      if(!bc_angle_incidence_elem)
      {
        err << "Missing BC_angle_dependence element" << bc_type_elem->Value() ;
        throw Dark_Arts_Exception(INPUT, err );
      }
      if(!bc_time_dependence_elem)
      {
        err << "Missing BC_time_dependence element in " << bc_type_elem->Value() ;
        throw Dark_Arts_Exception(INPUT, err );
      }
      if(!bc_value_elem)
      {
        err << "BC_value element missing in " << bc_type_elem->Value() ;
        throw Dark_Arts_Exception(INPUT , err );
      }
        
      /// get energy/frequency distribution dependence, check for valid input   
      if( m_number_groups > 1 )
      {
        TiXmlElement* bc_energy_dependence_elem = bc_type_elem->FirstChildElement( "BC_energy_dependence");
        if(!bc_energy_dependence_elem)
        {
          err<< "In BC_IC element: Missing BC_Energy_Dependence element in " << bc_type_elem->Value() ;
          throw Dark_Arts_Exception(INPUT, err); 
        }
         
        std::string energy_dep_str = bc_energy_dependence_elem->GetText();
        transform(energy_dep_str.begin() , energy_dep_str.end() , energy_dep_str.begin() , toupper);
        if(energy_dep_str == "PLANCKIAN")
          bc_data.energy_dependence = PLANCKIAN;
          
        if( bc_data.energy_dependence == INVALID_BC_ENERGY_DEPENDENCE)
        {
          err << "In BC_IC element: Invalid BC_ENERGY_DEPENDENCE in: " << bc_type_elem->Value();
          throw Dark_Arts_Exception(INPUT, err );
        }
      }    

      /// get incident value.  Check for basic physcaility. 
      bc_data.value = atof( bc_value_elem->GetText() ); 
      if( bc_data.value < 0.)
      {
        err << "Invalid value for BC_Value in " << bc_type_elem->Value();
        throw Dark_Arts_Exception(INPUT, err );
      }
      
      /// Determine what that numerical value we just read in is describing
      std::string value_type_str = bc_value_type_elem->GetText();
      transform(value_type_str.begin() , value_type_str.end() , value_type_str.begin() , toupper);
      if(value_type_str == "INCIDENT_CURRENT")  
        bc_data.value_type = INCIDENT_CURRENT;
      else if(value_type_str == "INCIDENT_TEMPERATURE")
        bc_data.value_type = INCIDENT_TEMPERATURE;
      
      if(bc_data.value_type == INVALID_INCIDENT_BC_VALUE_TYPE)
      {
        err << "Invalid BC_Value_Type element value in " << bc_type_elem->Value() ;
        throw Dark_Arts_Exception(INPUT, err ) ;
      }
      
      /// get angular dependence, and check for validity
      std::string bc_incidence_str = bc_angle_incidence_elem->GetText();
      transform(bc_incidence_str.begin() , bc_incidence_str.end() , bc_incidence_str.begin() , toupper);
      if(bc_incidence_str == "BC_ISOTROPIC")
      {
        bc_data.angle_dependence = BC_ISOTROPIC;
      }
      else if(bc_incidence_str == "BC_GLANCE")
      {
        bc_data.angle_dependence = BC_GLANCE;
      }
      else if(bc_incidence_str == "BC_NORMAL")
      {
        bc_data.angle_dependence = BC_NORMAL;
      }
      
      if(bc_data.angle_dependence == INVALID_BC_ANGLE_DEPENDENCE)
      {
        err << "Must specify angular dependence in " << bc_type_elem->Value();
        throw Dark_Arts_Exception(INPUT, err);
      }
      
      /// get left BC time dependence
      std::string time_dependence_str = bc_time_dependence_elem->GetText();
      transform(time_dependence_str.begin() , time_dependence_str.end() , time_dependence_str.begin() , toupper);
      if(time_dependence_str == "BC_BURST")
      {
        bc_data.time_dependence = BC_BURST;
      }
      else if(time_dependence_str == "BC_CONSTANT")
      {
        bc_data.time_dependence = BC_CONSTANT;
      }
      
      if(bc_data.time_dependence == INVALID_BC_TIME_DEPENDENCE)
      {
        err << "Invalid left BC time dependence in " << bc_type_elem->Value();
        throw Dark_Arts_Exception(INPUT, err );
      }
      else if(bc_data.time_dependence == BC_CONSTANT)
      {
        /// assume dirichlet conditions last forever
        bc_data.start_time = m_t_start - 0.2*(m_t_end - m_t_start);
        bc_data.end_time = m_t_end + 0.2*(m_t_end - m_t_start);  
      }
      else if(bc_data.time_dependence == BC_BURST)
      {
        TiXmlElement* bc_start_elem = bc_time_dependence_elem->FirstChildElement( "BC_turn_on" );
        TiXmlElement* bc_end_elem = bc_time_dependence_elem->FirstChildElement( "BC_turn_off" );
        
        if(!bc_start_elem || !bc_end_elem)
        {
          err << "BC_Burst requires BC_Turn_On and BC_Turn_Off in " << bc_type_elem->Value();
          throw Dark_Arts_Exception(INPUT, err);
        }
        bc_data.start_time = atof( bc_start_elem->GetText() );
        bc_data.end_time = atof( bc_end_elem->GetText() );
        
        if(bc_data.end_time < bc_data.start_time)
        {
          err << "BC_Turn_Off must be later (in time) than BC_Turn_Off in " << bc_type_elem->Value();
          throw Dark_Arts_Exception(INPUT, err );
        }  
        if( (bc_data.start_time < m_t_start) || (bc_data.end_time > m_t_end) )
        {
          err << "BC Turn_On / Turn_Off time must be within total problem times in " << bc_type_elem->Value();
          throw Dark_Arts_Exception(INPUT, err );
        }
      }
      break;
    }
    case MMS_BC:
    {
      break;
    }
    case REFLECTIVE_BC:
    {
      /// don't need any additional data, but this entry just feels good
      break;
    }
  }
   
  if( bc_data.type == REFLECTIVE_BC) 
  {
    std::string side_str = bc_type_elem->Value();
    /// only allow reflection on the left edge
    if( side_str == "Right_radiation_bc_type") 
      throw Dark_Arts_Exception(INPUT, "In BC_IC element: Cannot use reflective boundary condition on right edge.  Only left edge");
   
    /// do not allow refelective condition and krylov solving
    /// this is not coded yet / requires a not small amount of changes to incorporate
    if( (m_wg_solve_type == KRYLOV_SWEEPS) || (m_wg_solve_type == KRYLOV_DSA) ) 
      throw Dark_Arts_Exception(INPUT, "In BC_IC element:A refelective radiation boundary condition is not allowed with Krylov WGRS");
  }    
    
  return;
}
