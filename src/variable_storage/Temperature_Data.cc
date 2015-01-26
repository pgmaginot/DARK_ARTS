/** @file   Temperature_Data.cc
  *   @author pmaginot
  *   @brief Implement the Temperature_Data class
  *   Store temperature unknowns
*/
#include "Temperature_Data.h"


/// base constructor
Temperature_Data::Temperature_Data(const int n_cells, const Fem_Quadrature& fem_quad)
  /// initilaize range members
  : m_cells{ n_cells } ,     
    m_el_per_cell{fem_quad.get_number_of_interpolation_points() },
    m_t_length{ m_cells*m_el_per_cell} ,
    m_t(m_t_length,0.)    
  {
    fem_quad.get_dfem_interpolation_point_weights(m_dfem_w);
  }
  
/// initial condition constructor
Temperature_Data::Temperature_Data(const Fem_Quadrature& fem_quad, const Input_Reader& input_reader, const Cell_Data& cell_data)
  /// initilaize range members
  : m_cells{ cell_data.get_total_number_of_cells() } ,     
    m_el_per_cell{fem_quad.get_number_of_interpolation_points() },
    m_t_length{ m_cells*m_el_per_cell} ,
    m_t(m_t_length,0.)    
{
  fem_quad.get_dfem_interpolation_point_weights(m_dfem_w);
  if( input_reader.is_restart() )
  {
    std::cout << "This is a restart temperature\n";
    /// load Temperature_Data from a restart file
    std::string temperature_filename;
    input_reader.get_data_dump_path(temperature_filename);
    
    std::string input_file_name_and_path;
    input_reader.get_initial_input_filename(input_file_name_and_path);
    
    /// get only the input file name
    unsigned found = input_file_name_and_path.find_last_of("/");  
    temperature_filename.append( input_file_name_and_path.substr(found+1) );
    /// remove the .xml and add temperature restart file strings
    
    std::string xml_extension = ".xml";
    int time_step_resume = input_reader.get_time_step_restart();
    std::string temperature_extension= "_TemperatureDump_step_";
    temperature_extension += std::to_string(time_step_resume);
    temperature_extension.append(".xml");
    temperature_filename.replace(temperature_filename.find(xml_extension),xml_extension.length() , temperature_extension );
    
    TiXmlDocument doc( temperature_filename.c_str() );  
    bool loaded = doc.LoadFile();  
    if( loaded  )
      std::cout << "Found Temperature restart file: " << temperature_filename << std::endl;
    else
    {
      std::stringstream err;
      err << "Error reading Temperature restart file: " << temperature_filename << std::endl;
      throw Dark_Arts_Exception( INPUT ,  err.str() );
    }
    TiXmlElement* root = doc.FirstChildElement("TEMPERATURE");
    
    if(!root)
      throw Dark_Arts_Exception(INPUT, "Could not find TEMPERATURE block in Temperature restart file");
          
    TiXmlElement* cell_element = root->FirstChildElement("Cell");
    for(int c=0 ; c<m_cells ; c++)
    {
      if(!cell_element)
      {
        std::stringstream err;
        err << "Missing  cell element: " << c << " in temperature restart file";
        throw Dark_Arts_Exception(VARIABLE_STORAGE , err);
      }
        
      int cell_num = atoi(cell_element->GetText() );
      if( cell_num != c)
        throw Dark_Arts_Exception(INPUT , "Cell text not the same as expected cell number");
      
      Eigen::VectorXd local_t = Eigen::VectorXd::Zero(m_el_per_cell);
      
      TiXmlElement* element_element = cell_element->FirstChildElement("Element");    
      for(int el = 0 ; el < m_el_per_cell ; el++)
      {
        if(!element_element)
        {
          std::stringstream err;
          err << "Missing element: " << el << " of cell:  " << c << " in restart temperature loading";
          throw Dark_Arts_Exception(VARIABLE_STORAGE , err);
        }
        
        int el_num = atoi(element_element->GetText() );
        if( el_num != el)
        {
          std::stringstream err;
          err << "Elements out of order in cell: " << c << " of temperature restart file";
          throw Dark_Arts_Exception(VARIABLE_STORAGE, err);
        }
        
        TiXmlElement* val_element = element_element->FirstChildElement("Value");
        if(!val_element)
        {
          std::stringstream err;
          err << "Missing Value element for Element: " << el << " in cell: " << c << " of Temperature restart data file";
          throw Dark_Arts_Exception(VARIABLE_STORAGE, err);
        }
        
        local_t(el) = atof(val_element->GetText() );   
        if( isnan(local_t(el)) )
          throw Dark_Arts_Exception(VARIABLE_STORAGE , "Loading a NAN temperature");
          
        element_element = element_element->NextSiblingElement("Element");
      }
      /// save this data in the Temperature_Data structure
      set_cell_temperature(c,local_t);
      
      cell_element = cell_element->NextSiblingElement("Cell");      
    }
  }
  else{
    /// initialize first time step temperature data
    /// loop over regions.  Each region could have a different initial condition
    TEMPERATURE_IC_TYPE ic_type = input_reader.get_temperature_ic_type();
    switch( ic_type )
    {
      case CONSTANT_TEMPERATURE_IC:
      {
        std::vector<int> cell_per_reg;
        std::vector<double> temp_in_reg;
        input_reader.get_cells_per_region_vector(cell_per_reg);
          
        int cell_cnt = 0;
        for(int reg = 0 ; reg < input_reader.get_n_regions() ; reg++)
        {
          double temp_eval = input_reader.get_region_temperature(reg);
          int n_cell_reg = cell_per_reg[reg];     
          /// assume isotropic planck emission         
          for(int cell = 0; cell < n_cell_reg ; cell++)
          {
            set_cell_temperature( (cell+cell_cnt) , temp_eval );
          }
          cell_cnt += n_cell_reg;          
        }   
        break;
      }
      case MMS_TEMPERATURE_IC:
      {
        std::vector<int> cell_per_reg;      
        input_reader.get_cells_per_region_vector(cell_per_reg);           
        
        Eigen::VectorXd local_t;
        local_t = Eigen::VectorXd::Zero(m_el_per_cell);
        std::vector<double> dfem_loc;
        fem_quad.get_dfem_interpolation_point(dfem_loc);
        
        /// get radiation space angle and temporal components
        MMS_Temperature t_ic(input_reader);
        
        int cell_cnt = 0;
        const double time = input_reader.get_t_start();
        for(int reg = 0 ; reg < input_reader.get_n_regions() ; reg++)
        {
          int n_cell_reg = cell_per_reg[reg];        
          for(int c = 0; c < n_cell_reg ; c++)
          {
            int cell = c + cell_cnt;       
            double xL = cell_data.get_cell_left_edge(cell);
            double dx = cell_data.get_cell_width(cell);        
            for(int el = 0; el < m_el_per_cell ; el++)
            {
              double position = xL + dx/2.*( 1. + dfem_loc[el] );
              local_t(el) = t_ic.get_mms_temperature(position, time);
            }
            set_cell_temperature( cell , local_t);
          }
          cell_cnt += n_cell_reg;
        }    
        break;
      }
      case INVALID_TEMPERATURE_IC_TYPE:
      {
        throw Dark_Arts_Exception( VARIABLE_STORAGE , "Must have a valid Temperature IC initializtion");
        break;
      }
    }  
  }
}
  
/// Public accessor functions
double Temperature_Data::get_temperature(const int el, const int cell) const
{  
  return m_t[temperature_data_locator(el,cell)];
}

void Temperature_Data::set_temperature(const int el, const int cell, const double val)
{
  int loc = temperature_data_locator(el,cell);
  m_t[loc] = val;
  return ;
}

void Temperature_Data::get_cell_temperature(const int cell, Eigen::VectorXd& vec) const
{  
  int base = temperature_data_locator(0,cell);
  for(int i=0; i< m_el_per_cell; i++ )
    vec(i) = m_t[base + i];
  return; 
}

void Temperature_Data::set_cell_temperature(const int cell, const Eigen::VectorXd& vec)
{
  int loc = temperature_data_locator(0,cell);
  for(int i=0; i< m_el_per_cell ; i++)
    m_t[loc+i] = vec(i);
    
  return ;
}

void Temperature_Data::set_cell_temperature(const int cell, const double val)
{
  int loc = temperature_data_locator(0,cell);
  for(int i=0; i< m_el_per_cell ; i++)
    m_t[loc+i] = val;
    
  return ;
}


bool Temperature_Data::temperature_range_check(const int el, const int cell) const
{
  bool is_bad = false;
  
  if( (el >= m_el_per_cell ) || (el < 0) )
    is_bad = true;
       
  if( (cell < 0) || (cell >= m_cells) )
    is_bad = true;
  
  return is_bad;
}
  
int Temperature_Data::temperature_data_locator(const int el, const int cell) const
{
  int loc_val = -1;
  
  if( temperature_range_check(el,cell) )
    throw Dark_Arts_Exception( VARIABLE_STORAGE , "Attempting to access illogical temperature location");
   
  
  /// layout temperature unknowns from left to right
  loc_val = cell*m_el_per_cell + el;
  
  if( temperature_bounds_check(loc_val) )
    throw Dark_Arts_Exception( VARIABLE_STORAGE , "Calculated memory location out of bounds in temperature data");
  
  return loc_val;
}
   
bool Temperature_Data::temperature_bounds_check(const int loc) const
{
  bool is_bad_loc = false;
  if( (loc < 0) || (loc >= m_t_length) )
    is_bad_loc = true;  
  
  return is_bad_loc;
}

Temperature_Data& Temperature_Data::operator= (const Temperature_Data& t_data)
{
  if( (m_cells != t_data.m_cells) ||
      (m_el_per_cell != t_data.m_el_per_cell) ||
      (m_t_length != t_data.m_t_length) )
  {
    throw Dark_Arts_Exception( VARIABLE_STORAGE , "Trying to copy non-identical Temperature_Data objects");
  }
  for(int i=0; i< m_t_length; i++)
    m_t[i] = t_data.m_t[i];
  
  return *this;
}

double Temperature_Data::calculate_average(void)
{
  double val = 0.;
  int cnt = 0;
  for(int c=0; c<m_cells; c++)
  {
    for(int el=0;el<m_el_per_cell;el++)
    {
      val += m_t[cnt]*m_dfem_w[el];
      cnt++;
    }
  }
  val /= (2.* double(m_cells));
  return val;
}

