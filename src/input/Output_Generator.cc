/** @file   Input_Reader.cc
  *   @author pmaginot
  *   @brief Implement the Input_Reader class that reads XML files with TinyXml
 */
         
#include "Output_Generator.h"
// ##########################################################
// Public functions 
// ##########################################################


Output_Generator::Output_Generator(const Fem_Quadrature& fem_quadrature, const Cell_Data& cell_data, std::string input_filename)
:
 m_n_dfem( fem_quadrature.get_number_of_interpolation_points() ) , 
 m_n_cells( cell_data.get_total_number_of_cells() ),
 m_filename_base(input_filename)
{
}

void Output_Generator::create_output(const bool is_final, const double time, const int time_step)
{

  
  return;
}  
  
void Output_Generator::write_xml( std::string xmlFilename , const Temperature_Data& temperature_data)
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
  
  Eigen::VectorXd temperature;
  
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
      sprintf( buffer , "%24.14e" , temperature(el) );
      value_number = new TiXmlText(buffer);
      value_element->LinkEndChild( value_number );
      
      element_element->LinkEndChild( element_number );
      element_element->LinkEndChild( value_element );
      
      cell_element->LinkEndChild( element_element );
    }
    root->LinkEndChild( cell_element );
  }
  
  std::cout << "Should be done now \n";
  doc.SaveFile("filename.xml");
  
  return;
}

void Output_Generator::write_xml( std::string xmlFilename , const Intensity_Moment_Data& phi)
{

  return;
}

void Output_Generator::write_xml( std::string xmlFilename , const Intensity_Data& intensity)
{

  return;
}
