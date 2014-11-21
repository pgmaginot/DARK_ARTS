#include <iostream>
#include <vector>
#include <math.h>
#include <string>

#include "Input_Reader.h"
#include "Cell_Data.h"
#include "Dark_Arts_Exception.h"

int main(int argc, char** argv)
{
  const double tol = 1.0E-6;
  int val = 0;
  
  /// expected spacing we should be able to see
  const int factor = 4;
  const double incr_factor = double(factor);
  const int expected_n_cells = 8*incr_factor;
  std::vector<double> expected_left_edges(expected_n_cells,0.);
  std::vector<double> expected_cell_widths(expected_n_cells,0.);
  std::vector<int> expected_mat_number(expected_n_cells,-1);
  
  expected_cell_widths[0] = 1./incr_factor;
  expected_cell_widths[1] = 1./incr_factor;
  expected_cell_widths[2] = 1./incr_factor;
  expected_cell_widths[3] = 1./incr_factor;
  
  expected_cell_widths[4] = 2./incr_factor;
  expected_cell_widths[5] = 2./incr_factor;
  expected_cell_widths[6] = 2./incr_factor;
  expected_cell_widths[7] = 2./incr_factor;
  
  expected_cell_widths[8] = 4./incr_factor;
  expected_cell_widths[9] = 4./incr_factor;
  expected_cell_widths[10] = 4./incr_factor;
  expected_cell_widths[11] = 4./incr_factor;
  
  expected_cell_widths[12] = 8./incr_factor;
  expected_cell_widths[13] = 8./incr_factor;
  expected_cell_widths[14] = 8./incr_factor;
  expected_cell_widths[15] = 8./incr_factor;
  
  expected_cell_widths[16] = 16./incr_factor;
  expected_cell_widths[17] = 16./incr_factor;
  expected_cell_widths[18] = 16./incr_factor;
  expected_cell_widths[19] = 16./incr_factor;
  
  expected_cell_widths[20] = 4./incr_factor;
  expected_cell_widths[21] = 4./incr_factor;
  expected_cell_widths[22] = 4./incr_factor;
  expected_cell_widths[23] = 4./incr_factor;
  
  expected_cell_widths[24] = 1./incr_factor;
  expected_cell_widths[25] = 1./incr_factor;
  expected_cell_widths[26] = 1./incr_factor;
  expected_cell_widths[27] = 1./incr_factor;
  
  expected_cell_widths[28] = 0.25/incr_factor;  
  expected_cell_widths[29] = 0.25/incr_factor;  
  expected_cell_widths[30] = 0.25/incr_factor;  
  expected_cell_widths[31] = 0.25/incr_factor;  
  
  double dx_total = 0.;
  for(int i = 0 ; i < expected_n_cells ; i++)
  {
    expected_left_edges[i] = 0. + dx_total;
    dx_total += expected_cell_widths[i];
  }
  
  expected_mat_number[0] = 0;
  expected_mat_number[1] = 0;
  expected_mat_number[2] = 0;
  expected_mat_number[3] = 0;
  expected_mat_number[4] = 0;
  expected_mat_number[5] = 0;
  expected_mat_number[6] = 0;
  expected_mat_number[7] = 0;
  expected_mat_number[8] = 0;
  expected_mat_number[9] = 0;
  expected_mat_number[10] = 0;
  expected_mat_number[11] = 0;
  expected_mat_number[12] = 0;
  expected_mat_number[13] = 0;
  expected_mat_number[14] = 0;
  expected_mat_number[15] = 0;
  expected_mat_number[16] = 0;
  expected_mat_number[17] = 0;
  expected_mat_number[18] = 0;
  expected_mat_number[19] = 0;
  
  
  expected_mat_number[20] = 1;
  expected_mat_number[21] = 1;
  expected_mat_number[22] = 1;
  expected_mat_number[23] = 1;
  expected_mat_number[24] = 1;
  expected_mat_number[25] = 1;
  expected_mat_number[26] = 1;
  expected_mat_number[27] = 1;
  expected_mat_number[28] = 1;
  expected_mat_number[29] = 1;
  expected_mat_number[30] = 1;
  expected_mat_number[31] = 1;
  
  Input_Reader input_reader;    
    
  try{
    input_reader.read_xml(argv[1]);
  }
  catch( const Dark_Arts_Exception& da_exception )
  {
    da_exception.testing_message();
    val = -1;
  }
  std::cout << "Completed reading input file\n";
  
  Cell_Data cell_data(input_reader);
  
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
      int mat_num = cell_data.get_cell_material_number( i ) ;
      
      if( fabs(cell_width - expected_cell_widths[i] ) > tol )
      {
        std::stringstream err;
        err << "Cell " << i << " width expected to be: " << expected_cell_widths[i] << " actually is: " << cell_width ;
        throw Dark_Arts_Exception( SUPPORT_OBJECT , err.str() ) ; 
      }
      
      if( fabs(left_edge - expected_left_edges[i] ) > tol )
      {
        std::stringstream err;
        err << "Cell " << i << " left edge expected to be: " << expected_left_edges[i] << " actually is: " << left_edge ;
        throw Dark_Arts_Exception( SUPPORT_OBJECT , err.str() ) ; 
      }   
      if( (mat_num - expected_mat_number[i]) != 0 )
      {
        std::stringstream err;
        err << "Cell " << i << " Material number supposed to be: " << expected_mat_number[i] << " actually is: " << mat_num ;
        throw Dark_Arts_Exception( SUPPORT_OBJECT , err.str() ) ; 
      }   
    }
  }
  catch(const Dark_Arts_Exception& da_exception )
  {
    da_exception.testing_message();
  }

  
  // Return 0 if tests passed, somethnig else if failing
  return val;
}
