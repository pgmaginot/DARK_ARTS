/** @file   Cell_Data.cc
  *   @author pmaginot
  *   @brief Implement the Cell_Data class, holds all cell information (xL,xR,material_number)
*/
#include "Cell_Data.h"

Cell_Data::Cell_Data(Input_Reader&  input_reader)
{
  std::vector<int> cells_per_region;
  input_reader.get_cells_per_region_vector(cells_per_region);
  
  int n_region = input_reader.get_n_regions();
  m_total_cells = 0;
  for(int i=0;i<n_region;i++)
  {
    m_total_cells += cells_per_region[i];
  }
  
  /// initialize vectors of data
  m_x_l.resize(m_total_cells,0.);
  m_dx.resize(m_total_cells,0.);
  m_material_num.resize(m_total_cells,0);
    
  try{
      determine_cell_properties(n_region, cells_per_region, input_reader);
  }
  catch(const Dark_Arts_Exception& da_exception )
  { 
    da_exception.message();
  }
}

double Cell_Data::get_cell_width(int cell_num) const
{
  return m_dx[cell_num];
}

double Cell_Data::get_cell_left_edge(int cell_num) const
{
  return m_x_l[cell_num];
}

int Cell_Data::get_cell_material_number(int cell_num) const
{
  return m_material_num[cell_num];
}

int Cell_Data::get_total_number_of_cells(void) const
{
  return m_total_cells;
}

// void Cell_Data::read_mesh_and_refine(const Input_Reader& input_reader)
// {
  // std::cout << "In this file";
  // /** determine CellDataDump filename
   // trim .xml from initial filename, add _CellDataDump.xml
   // */
  // std::string celldata_filename;
  // input_reader.get_data_dump_path(celldata_filename);
  
  // std::string input_file_name_and_path;
  // input_reader.get_initial_input_filename(input_file_name_and_path);
  
  // /// get only the input file name
  // unsigned found = input_file_name_and_path.find_last_of("/");  
  // celldata_filename.append( input_file_name_and_path.substr(found+1) );
  // /// remove the .xml and add _CellDataDump.xml
  
  // std::string xml_extension = ".xml";
  // celldata_filename.replace(celldata_filename.find(xml_extension),xml_extension.length() , "_CellDataDump.xml" );
  
  // TiXmlDocument doc( celldata_filename.c_str() );  
  // bool loaded = doc.LoadFile();  
  // if( loaded  )
    // std::cout << "Found CellDataDump file: " << celldata_filename << std::endl;
  // else
  // {
    // std::stringstream err;
    // err << "Error reading CellDataDump file: " << celldata_filename << std::endl;
    // throw Dark_Arts_Exception( INPUT ,  err.str() );
  // }
  // TiXmlElement* root = doc.FirstChildElement("CELL_DATA");
  
  // if(!root)
    // throw Dark_Arts_Exception(INPUT, "Could not find CELL_DATA block in CellDataDump file");
    
  // const int incr_factor = input_reader.get_refinement_factor();
  // const int n_original_cells = m_total_cells/incr_factor;
    
  // TiXmlElement* cell_element = root->FirstChildElement("Cell");
  // int new_cell_cnt = 0;
  // for(int c=0 ; c<n_original_cells ; c++)
  // {
    // int old_cell_num = atoi(cell_element->GetText() );
    // if( old_cell_num != c)
      // throw Dark_Arts_Exception(INPUT , "Cell text not the same as expected cell number");
    
    // TiXmlElement* material_element = cell_element->FirstChildElement("Material");    
    // TiXmlElement* dx_old_element = cell_element->FirstChildElement("Cell_width");
    // TiXmlElement* x_left_old_element = cell_element->FirstChildElement("X_left");
      
    // double dx_old = atof(dx_old_element->GetText() );
    // int material_old = atoi(material_element->GetText() );
    // double x_left_old = atof(x_left_old_element->GetText() );
    
    // if(dx_old < 0.)
    // {
      // std::stringstream err;
      // err << "Found dx < 0 in CellDataDump that is being read for cell: " << c ;
      // throw Dark_Arts_Exception(INPUT, err.str() );
    // }
    // double dx_new = dx_old/(double (incr_factor));
    // const double dbl_incr_factor = double(incr_factor);
    // for(double n = 0. ; n < dbl_incr_factor ; n+=1. )
    // {
      // m_x_l[new_cell_cnt] = x_left_old + dx_new*n;
      // m_dx[new_cell_cnt] = dx_new;
      // m_material_num[new_cell_cnt] = material_old;
      
      // new_cell_cnt++;
    // }  
    
    // cell_element = cell_element->NextSiblingElement("Cell");
  // }
   
  // return;
// }

void Cell_Data::determine_cell_properties(const int n_reg, const std::vector<int>& cell_reg,
  const Input_Reader& input_reader)
{
  /// looop over the number of regions
  int cell_cnt = 0;
  for(int i=0;i<n_reg;i++)
  {
    /// Get the relevant data for this region
    int n_cell = cell_reg[i];
    int mat_num = input_reader.get_region_material_number(i);
    double x_l = input_reader.get_region_left_bound(i);
    double x_r = input_reader.get_region_right_bound(i);
    double x_l_cell = x_l;
    GRID_SPACING spacing = input_reader.get_region_spacing(i);
    
    if(spacing == EQUAL)
    {
      /// loop over the number of cells in the region
      double dx = (x_r - x_l)/( double(n_cell) );
      
      for(int c=0;c<n_cell;c++)
      {
        m_x_l[cell_cnt] = x_l_cell;
        m_dx[cell_cnt] = dx;
        x_l_cell += dx;
        m_material_num[cell_cnt] = mat_num;
        cell_cnt++;
      }
    }
    else if(spacing == LOG)
    {
      /// calculate leftmost cell width, assuming purely geometric series
      double r = input_reader.get_r_factor(i);
      double min_size = input_reader.get_min_cell_size(i);
      
      double dx_leftmost = (x_r - x_l)*(1. - r)/(1. - pow(r,n_cell)) ;
      double dx_rightmost = dx_leftmost* pow(r , n_cell-1);
      
      double x_l_cell = x_l;
      
      if(dx_leftmost < 0.)
        throw Dark_Arts_Exception(SUPPORT_OBJECT , "Calculated a negative cell width in log spacing");
      
      if(dx_rightmost < 0.)
        throw Dark_Arts_Exception(SUPPORT_OBJECT , "Calculated a negative cell width in log spacing");
      
      /// enforce minimum cell width condition
      if( (dx_leftmost < min_size) && (dx_rightmost < min_size) )
      {
        std::cout << "dx_left: " << dx_leftmost << "\n dx_right: " << dx_rightmost << std::endl;
        throw Dark_Arts_Exception(SUPPORT_OBJECT ,  "Minimum cell size is too large for the desired number of cells");        
      }
      else if( dx_leftmost < min_size) 
      {
        /// cells on the left of the region are too small
        
        double dx = dx_leftmost;
        int c=0;
        for( ; c<n_cell ; c++)
        {
          if(dx < min_size)
          {
            m_x_l[cell_cnt] = x_l_cell;
            m_dx[cell_cnt] = min_size;
            m_material_num[cell_cnt] = mat_num;
            x_l_cell += min_size;
            cell_cnt++;
            dx *=r;
          }
          else{
            /// c cells were adjusted, fill the space
            break;
          }          
        }
        
        dx = (x_r - x_l_cell)*(1.-r)/(1. - pow(r,n_cell - c));
        for( ; c< n_cell ; c++)
        {
          m_x_l[cell_cnt] = x_l_cell;
          m_dx[cell_cnt] = dx;
          m_material_num[cell_cnt] = mat_num;
          
          /// get ready for the next cell
          x_l_cell += dx;
          dx *= r;          
          cell_cnt++;
        }        
      }
      else if( dx_rightmost < min_size)
      {
        /// cells on the right of the region are too small
        int c=0;
        /// count the number of cells that would be too small using pure log spacing
        double dx = dx_rightmost;
        for( ; c<n_cell; c++)
        {
          if( dx > min_size)
            break;
            
          dx /=r;
        }
        /// c now holds the number of cells that were adjusted (enlarged)
        double x_r_new = x_r - double(c) * min_size;
        dx = (x_r_new - x_l)*(1.-r)/(1.-pow(r,n_cell - c));
        
        /// fill in the cells that will actually be log spaced
        int cell = 0;
        for( ; cell < n_cell - c; cell++)
        {
          m_x_l[cell_cnt] = x_l_cell;
          m_dx[cell_cnt] = dx;
          m_material_num[cell_cnt] = mat_num;
          
          x_l_cell += dx;
          dx *=r;
          cell_cnt++;
        }
        /// fill in the cells that will only be the minimum size
        dx = min_size;
        for( ; cell <n_cell ; cell++)
        {
          m_x_l[cell_cnt] = x_l_cell;
          m_dx[cell_cnt] = dx;
          m_material_num[cell_cnt] = mat_num;
          
          x_l_cell += dx;
          cell_cnt++;
        }
        
      }
      else
      {      
        /// log spacing works as defined
        double dx = dx_leftmost;
        for(int c=0; c<n_cell; c++)
        {
          /// save this cell's data
          m_x_l[cell_cnt] = x_l_cell;
          m_dx[cell_cnt] = dx;
          m_material_num[cell_cnt] = mat_num;
          
          /// get ready for the next cell
          x_l_cell += dx;
          dx *= r;          
          cell_cnt++;
        }
      }
      
    }
    

  }
  
  return;
}


