#include <iostream>
#include <vector>
#include <math.h>
#include <string>

#include "Input_Reader.h"
#include "Cell_Data.h"
#include "Output_Generator.h"
#include "Angular_Quadrature.h"
#include "Fem_Quadrature.h"
#include "Quadrule_New.h"
#include "Dark_Arts_Exception.h"

int main(int argc, char** argv)
{
  int val = 0;
  
  /// expected spacing we should be able to see
  const int expected_n_cells = 8;
  std::vector<double> expected_left_edges(expected_n_cells,0.);
  std::vector<double> expected_cell_widths(expected_n_cells,0.);
  std::vector<int> expected_mat_number(expected_n_cells,-1);
  
  expected_left_edges[0] = 0.;
  expected_left_edges[1] = 1.;
  expected_left_edges[2] = 3.;
  expected_left_edges[3] = 7.;
  expected_left_edges[4] = 15.;
  expected_left_edges[5] = 31.;
  expected_left_edges[6] = 35.;
  expected_left_edges[7] = 36.;
  
  expected_cell_widths[0] = 1.;
  expected_cell_widths[1] = 2.;
  expected_cell_widths[2] = 4.;
  expected_cell_widths[3] = 8.;
  expected_cell_widths[4] = 16.;
  expected_cell_widths[5] = 4.;
  expected_cell_widths[6] = 1.;
  expected_cell_widths[7] = 0.25;  
  
  expected_mat_number[0] = 0;
  expected_mat_number[1] = 0;
  expected_mat_number[2] = 0;
  expected_mat_number[3] = 0;
  expected_mat_number[4] = 0;
  expected_mat_number[5] = 1;
  expected_mat_number[6] = 1;
  expected_mat_number[7] = 1;
  
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
  
  /// Do all this to dump an output file that I want to use for LogSpacingRefinementTest
  Cell_Data cell_data( input_reader ); 
  Quadrule_New quad_fun;  
  Fem_Quadrature fem_quadrature( input_reader , quad_fun);  
  Angular_Quadrature angular_quadrature( input_reader , quad_fun );  
  Output_Generator output( angular_quadrature, fem_quadrature, cell_data, argv[1] );
  
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
