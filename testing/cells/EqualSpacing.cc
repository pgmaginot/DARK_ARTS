#include <iostream>
#include <vector>
#include <math.h>
#include <string>

#include "Input_Reader.h"
#include "Cell_Data.h"
#include "Dark_Arts_Exception.h"

int main(int argc, char** argv)
{
  int val = 0;
  
  /// expected spacing we should be able to see
  const int expected_n_cells = 5;
  const double dx = (22.-1.)/double(expected_n_cells);
  std::vector<double> expected_left_edges(expected_n_cells,0.);
  std::vector<double> expected_cell_widths(expected_n_cells,0.);
  
  expected_left_edges[0] = 1.;
  expected_left_edges[1] = 1. + dx;
  expected_left_edges[2] = 1. + 2.*dx;
  expected_left_edges[3] = 22. - 2.*dx;
  expected_left_edges[4] = 22. - dx;
  
  expected_cell_widths[0] = dx;
  expected_cell_widths[1] = dx;
  expected_cell_widths[2] = dx;
  expected_cell_widths[3] = dx;
  expected_cell_widths[4] = dx;
  
  Input_Reader input_reader;    
    
  try{
    input_reader.read_xml(argv[1]);
  }
  catch( const Dark_Arts_Exception& da_exception )
  {
    da_exception.testing_message();
    val = -1;
  }
  
  Cell_Data cell_data( input_reader ); 
  try{
    int n_cells = cell_data.get_total_number_of_cells();
    if( n_cells != expected_n_cells)
    {
      std::stringstream err;
      err << "Did not return expected number of cells in LogSpacing.cc Test";
      throw Dark_Arts_Exception( SUPPORT_OBJECT , err.str() ) ;
    }
    
    for(int i=0; i < expected_n_cells ; i++)
    {
      double cell_width = cell_data.get_cell_width( i ) ;
      double left_edge = cell_data.get_cell_left_edge( i ) ;
      
      if( fabs(cell_width - expected_cell_widths[i] ) > 0.01 )
      {
        std::stringstream err;
        err << "Cell " << i << " width expected to be: " << expected_cell_widths[i] << " actually is: " << cell_width ;
        throw Dark_Arts_Exception( SUPPORT_OBJECT , err.str() ) ; 
      }
      
      if( fabs(left_edge - expected_left_edges[i] ) > 0.01 )
      {
        std::stringstream err;
        err << "Cell " << i << " left edge expected to be: " << expected_left_edges[i] << " actually is: " << left_edge ;
        throw Dark_Arts_Exception( SUPPORT_OBJECT , err.str() ) ; 
      }   
    }
  }
  catch(const Dark_Arts_Exception& da_exception )
  {
    da_exception.testing_message();
    val = -1;
  }
  
  // Return 0 if tests passed, somethnig else if failing
  return val;
}
