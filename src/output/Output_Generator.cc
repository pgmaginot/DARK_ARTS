/** @file   Input_Reader.cc
  *   @author pmaginot
  *   @brief Implement the Input_Reader class that reads XML files with TinyXml
 */
         
#include "Output_Generator.h"
// ##########################################################
// Public functions 
// ##########################################################


Output_Generator::Output_Generator(const Angular_Quadrature& angular_quadrature,
    const Fem_Quadrature& fem_quadrature, 
    const Cell_Data& cell_data, 
    const Input_Reader& input_reader)
:
 m_n_dfem( fem_quadrature.get_number_of_interpolation_points() ) , 
 m_n_cells( cell_data.get_total_number_of_cells() ),
 m_n_groups( angular_quadrature.get_number_of_groups()),
 m_n_l_mom( angular_quadrature.get_number_of_leg_moments() ),
 m_n_dir( angular_quadrature.get_number_of_dir() ),
 m_cell_data( cell_data )
{
  input_reader.get_output_directory(m_filename);
  if(input_reader.is_mesh_refinement() )
  {
    std::string base_file_with_path = input_reader.get_initial_input_filename();
    unsigned int found = base_file_with_path.find_last_of("/");
    std::string base_short = base_file_with_path.substr(found+1);  
      
    m_filename.append(base_short);
    std::string xml_ext = ".xml";
    m_filename.replace( m_filename.find(xml_ext) , xml_ext.length() , "_");
    std::string short_input; 
    input_reader.get_short_input_filename(short_input);
    m_filename.append(short_input);
  }  
  else
  {
    std::string short_input; 
    input_reader.get_short_input_filename(short_input);
    /// make sure output 
    m_filename+=short_input;
  }

  /// short_input is the called filename e.g. dark_arts ../inputs/filename.xml
  /// if this is a refinement run, we are going to want to store append the base name too
  
  output_cell_data();
  if(input_reader.get_output_type() == DUMP)
    output_cell_data_text(fem_quadrature);
}

void Output_Generator::write_xml(const bool is_final, const int time_step, const Temperature_Data& temperature_data)
{
  TiXmlDocument doc;
  TiXmlDeclaration* decl = new TiXmlDeclaration( "1.0", "", "" );
  doc.LinkEndChild( decl );
  
  TiXmlElement* root = new TiXmlElement( "TEMPERATURE");  
  doc.LinkEndChild( root ); 
  
  /**
    Want the file to look like this
    <Cell> 0
      <Element> 0
        <Value> val </Value>
      </Element>
      <Element> 1
        <Value> val </Value>
      </Element>
      <Element> 2
        <Value> val </Value>
      </Element>
    </Cell> 
    <Cell> 1
      <Element> 0
        <Value> val </Value>
      </Element>
      <Element> 1
        <Value> val </Value>
      </Element>
      <Element> 2
        <Value> val </Value>
      </Element>
    </Cell> 
  */
  
  TiXmlElement* cell_element;
  TiXmlElement* element_element;
  TiXmlElement* value_element;
  
  TiXmlText* cell_number;
  TiXmlText* element_number;
  TiXmlText* value_number;
  
  Eigen::VectorXd temperature(m_n_dfem);
  
  char buffer [25];
  
  for(int c = 0; c < m_n_cells ; c++)
  {
    temperature_data.get_cell_temperature(c , temperature);
    
    cell_element = new TiXmlElement( "Cell");
    sprintf( buffer , "%i" , c);
    cell_number = new TiXmlText( buffer );
    cell_element->LinkEndChild(cell_number);
    
    for(int el = 0; el < m_n_dfem ; el++)
    {
      element_element = new TiXmlElement("Element");
      sprintf( buffer , "%i" , el);
      element_number = new TiXmlText(buffer);
      
      value_element = new TiXmlElement( "Value" );
      
      if(isnan(temperature(el) ) )
        throw Dark_Arts_Exception(VARIABLE_STORAGE , "Attempting to write a NAN temperature");
      
      sprintf( buffer , "%22.14e" , temperature(el) );
      value_number = new TiXmlText(buffer);
      value_element->LinkEndChild( value_number );
      
      element_element->LinkEndChild( element_number );
      element_element->LinkEndChild( value_element );
      
      cell_element->LinkEndChild( element_element );
    }
    root->LinkEndChild( cell_element );
  }
  
  std::string filename;
  construct_filename(0,is_final,time_step,filename);
  
  doc.SaveFile(filename.data() );
  
  return;
}

void Output_Generator::write_xml( const bool is_final, const int time_step, const Intensity_Data& intensity)
{
  TiXmlDocument doc;
  TiXmlDeclaration* decl = new TiXmlDeclaration( "1.0", "", "" );
  doc.LinkEndChild( decl );
  
  TiXmlElement* root = new TiXmlElement( "INTENSITY");  
  
  doc.LinkEndChild( root ); 
  
  /**
    Want the file to look like this
    <Cell> 0
      <Group> 0 
        <Direction> 0 
          <Element>0
            <Value> </Value>
          </Element>
          <Element>1
            <Value> </Value>
          </Element>
        </Direction>
        <Direction> 1 
          <Element>0
            <Value> </Value>
          </Element>
          <Element>1
            <Value> </Value>
          </Element>
        </Direction>
      </Group>
      <Group> 1 
        <Direction> 0 
          <Element>0
            <Value> </Value>
          </Element>
          <Element>1
            <Value> </Value>
          </Element>
        </Direction>
        <Direction> 1 
          <Element>0
            <Value> </Value>
          </Element>
          <Element>1
            <Value> </Value>
          </Element>
        </Direction>
      </Group>
    </Cell>
  */
  
  TiXmlElement* cell_element;
  TiXmlElement* group_element;
  TiXmlElement* direction_element;
  TiXmlElement* element_element;
  TiXmlElement* value_element;
  
  TiXmlText* cell_number;
  TiXmlText* group_number;
  TiXmlText* direction_number;
  TiXmlText* element_number;
  TiXmlText* value_number;
  
  Eigen::VectorXd i_vec(m_n_dfem);
  
  char buffer [25];
  
  for(int c = 0; c < m_n_cells ; c++)
  {
    cell_element = new TiXmlElement( "Cell");
    sprintf( buffer , "%i" , c);
    cell_number = new TiXmlText( buffer );
    cell_element->LinkEndChild(cell_number);
    for(int g=0; g < m_n_groups ; g++)
    {
      group_element = new TiXmlElement( "Group");
      sprintf( buffer , "%i" , g);
      group_number = new TiXmlText( buffer );
      group_element->LinkEndChild(group_number);
      
      for(int d = 0 ; d < m_n_dir ; d++)
      {
        direction_element = new TiXmlElement("Direction");
        sprintf( buffer , "%i" , d);
        direction_number = new TiXmlText( buffer );
        direction_element->LinkEndChild(direction_number);
        
        intensity.get_cell_intensity(c,g,d,i_vec );
        for(int el = 0; el < m_n_dfem ; el++)
        {
          element_element = new TiXmlElement("Element");
          sprintf( buffer , "%i" , el);
          element_number = new TiXmlText(buffer);
          element_element->LinkEndChild( element_number );
          
          value_element = new TiXmlElement( "Value" );
          if(isnan(i_vec(el) ) )
            throw Dark_Arts_Exception(VARIABLE_STORAGE , "Attempting to write a NAN temperature");
          
          sprintf( buffer , "%22.14e" , i_vec(el) );
          value_number = new TiXmlText(buffer);
          value_element->LinkEndChild( value_number );          
          
          element_element->LinkEndChild( value_element );
          
          direction_element->LinkEndChild(element_element);
        }
        group_element->LinkEndChild(direction_element);
      }
      cell_element->LinkEndChild(group_element);
    }
    root->LinkEndChild(cell_element);
  }
  std::string filename;
  construct_filename(1,is_final,time_step,filename);
  
  doc.SaveFile(filename.data() );
  return;
}

void Output_Generator::write_xml( const bool is_final, const int time_step, const Intensity_Moment_Data& phi)
{
  TiXmlDocument doc;
  TiXmlDeclaration* decl = new TiXmlDeclaration( "1.0", "", "" );
  doc.LinkEndChild( decl );
  
  TiXmlElement* root = new TiXmlElement( "ANGLE_INTEGRATED_INTENSITY");  
  
  doc.LinkEndChild( root ); 
  
  /**
    Want the file to look like this
    <Cell> 0
      <Group> 0 
        <Legendre_Moment> 0 
          <Element>0
            <Value> </Value>
          </Element>
          <Element>1
            <Value> </Value>
          </Element>
        </Legendre_Moment>
        <Legendre_Moment> 1 
          <Element>0
            <Value> </Value>
          </Element>
          <Element>1
            <Value> </Value>
          </Element>
        </Legendre_Moment>
      </Group>
      <Group> 1 
        <Legendre_Moment> 0 
          <Element>0
            <Value> </Value>
          </Element>
          <Element>1
            <Value> </Value>
          </Element>
        </Legendre_Moment>
        <Legendre_Moment> 1 
          <Element>0
            <Value> </Value>
          </Element>
          <Element>1
            <Value> </Value>
          </Element>
        </Legendre_Moment>
      </Group>
    </Cell>
  */
  
  TiXmlElement* cell_element;
  TiXmlElement* group_element;
  TiXmlElement* legendre_element;
  TiXmlElement* element_element;
  TiXmlElement* value_element;
  
  TiXmlText* cell_number;
  TiXmlText* group_number;
  TiXmlText* legendre_number;
  TiXmlText* element_number;
  TiXmlText* value_number;
  
  Eigen::VectorXd phi_vec(m_n_dfem);
  
  char buffer [25];
  
  for(int c = 0; c < m_n_cells ; c++)
  {
    cell_element = new TiXmlElement( "Cell");
    sprintf( buffer , "%i" , c);
    cell_number = new TiXmlText( buffer );
    cell_element->LinkEndChild(cell_number);
    for(int g=0; g < m_n_groups ; g++)
    {
      group_element = new TiXmlElement( "Group");
      sprintf( buffer , "%i" , g);
      group_number = new TiXmlText( buffer );
      group_element->LinkEndChild(group_number);
      
      for(int l_mom = 0 ; l_mom < m_n_l_mom ; l_mom++)
      {
        legendre_element = new TiXmlElement("Legendre_Moment");
        sprintf( buffer , "%i" , l_mom);
        legendre_number = new TiXmlText( buffer );
        legendre_element->LinkEndChild(legendre_number);
        
        phi.get_cell_angle_integrated_intensity(c,g,l_mom,phi_vec );
        for(int el = 0; el < m_n_dfem ; el++)
        {
          element_element = new TiXmlElement("Element");
          sprintf( buffer , "%i" , el);
          element_number = new TiXmlText(buffer);
          element_element->LinkEndChild( element_number );
          
          value_element = new TiXmlElement( "Value" );
          sprintf( buffer , "%22.14e" , phi_vec(el) );
          value_number = new TiXmlText(buffer);
          value_element->LinkEndChild( value_number );          
          
          element_element->LinkEndChild( value_element );
          
          legendre_element->LinkEndChild(element_element);
        }
        group_element->LinkEndChild(legendre_element);
      }
      cell_element->LinkEndChild(group_element);
    }
    root->LinkEndChild(cell_element);
  }
  std::string filename;
  construct_filename(2,is_final,time_step,filename);
  
  doc.SaveFile(filename.data() );
  
  return;
}

void Output_Generator::write_txt(const bool is_final, const int time_step, const Temperature_Data& temperature)
{
  std::string xml_extension = ".xml";
  std::string filename;
  construct_filename(0,is_final,time_step,filename);
  filename.replace(filename.find(xml_extension), xml_extension.length() , ".txt" );
  
  std::ofstream output{ filename , std::ofstream::out};

  std::stringstream err;
  err << "Could not open Temperature_Data output file: " << filename;
  if(!output)
    throw Dark_Arts_Exception(SUPPORT_OBJECT, err.str() ); 
    
  output << "#Cell \n" ;
  for(int el = 0; el < m_n_dfem ; el++)
    output << "# Element_"<< el << std::endl;
    
  Eigen::VectorXd t_vec(m_n_dfem);
  for(int c = 0; c < m_n_cells ; c++)
  {
    output<< c  << std::endl;
    temperature.get_cell_temperature(c,t_vec );
    for(int el = 0; el < m_n_dfem ; el++)
    {
      output << std::scientific << std::setprecision(15) << t_vec(el) << std::endl;
    }
    output << std::endl;
  }    
  output.close();
  return;
}  

void Output_Generator::write_txt(const bool is_final, const int time_step, const Intensity_Data& intensity) 
{
  std::string xml_extension = ".xml";
  std::string filename;
  construct_filename(1,is_final,time_step,filename);
  filename.replace(filename.find(xml_extension), xml_extension.length() , ".txt" );
  
  std::ofstream output{ filename , std::ofstream::out};

  std::stringstream err;
  err << "Could not open Intensity_Data output file: " << filename;
  if(!output)
    throw Dark_Arts_Exception(SUPPORT_OBJECT, err.str() ); 
    
  output << "#Cell  Group   Direction \n" ;
  for(int el = 0; el < m_n_dfem; el++)
    output << "# Element_"<< el << std::endl;
    
  Eigen::VectorXd i_vec(m_n_dfem);
  for(int c = 0; c < m_n_cells ; c++)
  {
    for(int g=0; g < m_n_groups ; g++)
    {
      for(int d = 0 ; d < m_n_dir ; d++)
      {        
        output<< std::setw(5) << c << "  " << g << "  " << d << std::endl;
        intensity.get_cell_intensity(c,g,d,i_vec );
        for(int el = 0; el < m_n_dfem ; el++)
        {
          output << std::scientific << std::setprecision(15) << i_vec(el) << std::endl;
        }
        output << std::endl;
      }
    }
  }
    
  output.close();
  
  return;
}  

void Output_Generator::write_txt(const bool is_final, const int time_step, const Intensity_Moment_Data& phi)
{
  std::string xml_extension = ".xml";
  std::string filename;
  construct_filename(2,is_final,time_step,filename);
  filename.replace(filename.find(xml_extension), xml_extension.length() , ".txt" );
  
  std::ofstream output{ filename , std::ofstream::out};

  std::stringstream err;
  err << "Could not open Intensity_Moment_Data output file: " << filename;
  if(!output)
    throw Dark_Arts_Exception(SUPPORT_OBJECT, err.str() ); 
    
  output << "#Cell    Group    L_mom \n" ;
  for(int el = 0; el < m_n_dfem; el++)
    output << "# Element_"<< el << std::endl;
    
  Eigen::VectorXd phi_vec(m_n_dfem);
  for(int c = 0; c < m_n_cells ; c++)
  {
    for(int g=0; g < m_n_groups ; g++)
    {
      for(int l_mom = 0 ; l_mom < m_n_l_mom ; l_mom++)
      {        
        output<< std::setw(5) << c << "  " << g << "  " << l_mom << std::endl;
        phi.get_cell_angle_integrated_intensity(c,g,l_mom,phi_vec );
        for(int el = 0; el < m_n_dfem ; el++)
        {
          output << std::scientific << std::setprecision(15) << phi_vec(el) << std::endl;
        }
        output << std::endl;
      }
    }
  }
    
  output.close();
  
  return;
}  

void Output_Generator::output_cell_data(void)
{
  /// output cell spacing and DFEM points.  Do this in XML first.
  TiXmlDocument doc;
  TiXmlDeclaration* decl = new TiXmlDeclaration( "1.0", "", "" );
  doc.LinkEndChild( decl );
  
  TiXmlElement* root = new TiXmlElement( "CELL_DATA");    
  doc.LinkEndChild( root ); 
  
  TiXmlElement* cell_element;
  TiXmlElement* material_element;
  TiXmlElement* dx_element;
  TiXmlElement* x_left_element;
  
  TiXmlText* cell_number;
  TiXmlText* material_number;
  TiXmlText* dx_number;
  TiXmlText* x_left_number;
  
  char buffer [25];
  
  for(int c = 0; c < m_n_cells ; c++)
  {
    cell_element = new TiXmlElement( "Cell");
    sprintf( buffer , "%i" , c);
    cell_number = new TiXmlText( buffer );
    cell_element->LinkEndChild(cell_number);
  
    material_element = new TiXmlElement( "Material" );
    dx_element = new TiXmlElement("Cell_width" );
    x_left_element = new TiXmlElement( "X_left");
    
    sprintf( buffer , "%i" , m_cell_data.get_cell_material_number(c) );
    material_number = new TiXmlText( buffer);
    material_element->LinkEndChild(material_number);
    
    sprintf( buffer , "%22.14e" , m_cell_data.get_cell_width(c) );
    dx_number = new TiXmlText( buffer);
    dx_element->LinkEndChild(dx_number);
    
    sprintf( buffer , "%22.14e" , m_cell_data.get_cell_left_edge(c) );
    x_left_number = new TiXmlText( buffer);
    x_left_element->LinkEndChild(x_left_number);
  
    cell_element->LinkEndChild(material_element);
    cell_element->LinkEndChild(dx_element);
    cell_element->LinkEndChild(x_left_element);
    root->LinkEndChild(cell_element);
  }
  std::string filename;
  construct_filename(3,false,1,filename);
  
  doc.SaveFile(filename.data() );
  
  return;
}

void Output_Generator::construct_filename( const int data_type , const bool is_final, const int ts, std::string& output_filename) const
{  
  std::string xml_extension = ".xml";
  std::stringstream new_str;
  switch(data_type)
  {
    case 0:
    {
      new_str << "_TemperatureDump";
      break;
    }
    case 1:
    {
      new_str << "_IntensityDump";
      break;
    }
    case 2:
    {
      new_str << "_PhiDump";
      break;
    }
    case 3:
    {
      new_str << "_CellDataDump";
      break;
    } 
    default:
    {
      throw Dark_Arts_Exception(SUPPORT_OBJECT , "Invalid data type integer in output file generator");
      break;
    }
  }
  
  if(data_type != 3)
  {
    if(is_final)
        new_str << "_Final";
    else
        new_str << "_step_"<< ts;
  }

  new_str << ".xml";
  
  output_filename = m_filename;
  output_filename.replace(output_filename.find(xml_extension),xml_extension.length() , new_str.str() );

  return;
}

void Output_Generator::output_cell_data_text(const Fem_Quadrature& fem_quadrature)
{
  /// Dump text files for MATLAB plotting
  std::string filename;
  construct_filename(3,false,1,filename);
  
  std::string xml_extension;
  xml_extension = ".xml";
  std::string new_str;
  new_str = ".txt";
  filename.replace(filename.find(xml_extension),xml_extension.length() , new_str );
  
  std::ofstream output{ filename , std::ofstream::out};

  std::stringstream err;
  err << "Could not open Cell_Data text file: " << filename;
  if(!output)
    throw Dark_Arts_Exception(SUPPORT_OBJECT, err.str() ); 
    
  output << "#Cell\n" ;
  for(int el = 0; el < m_n_dfem; el++)
    output << "# Element_"<< el << std::endl;
    
  double dx = 0.;
  double xL = 0.;
  
  std::vector<double> fem_pt;
  fem_quadrature.get_dfem_interpolation_point(fem_pt);
  
  for(int c = 0; c < m_n_cells ; c++)
  {
    output<< std::setw(5) << c << std::endl;
    dx = m_cell_data.get_cell_width(c);
    xL = m_cell_data.get_cell_left_edge(c);
    for(int el = 0; el < m_n_dfem ; el++)
    {
      output << std::scientific << std::setprecision(15) << xL + dx/2.*(1. + fem_pt[el]) << std::endl;
    }
    output << std::endl;
  }
    
  output.close();
  
  
  return;
}
