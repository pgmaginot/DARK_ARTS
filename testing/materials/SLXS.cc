#include <iostream>
#include <vector>
#include <math.h>
#include <string>

#include "Input_Reader.h"
#include "Materials.h"
#include "Fem_Quadrature.h"
#include "Quadrule_New.h"
#include "Cell_Data.h"
#include "Angular_Quadrature.h"
#include "Dark_Arts_Exception.h"

int main(int argc, char** argv)
{
  int val = 0;
  const double tol = 1.E-6;
  
  Input_Reader input_reader;    
  try
  {
    input_reader.read_xml(argv[1]);
  }
  catch(const Dark_Arts_Exception& da_exception )
  {
    da_exception.message() ;
    val = -1;
  }       
  /// Initialize a Quadrule object to be able to get all of the quadrature we need
  Quadrule_New quad_fun;  
  Fem_Quadrature fem_quadrature( input_reader , quad_fun);  
  Cell_Data cell_data( input_reader );  
  Angular_Quadrature angular_quadrature( input_reader , quad_fun );  
  /// Create a Materials object that contains all opacity, heat capacity, and source objects
  try
  {
    Materials materials( input_reader, fem_quadrature , cell_data, angular_quadrature );  
    
    /// verify that we get expected (constant values for material properties evaluations)
    const double cv_1 = 77.1;
    const double cv_2 = 3.0;
    
    const double sig_a_1 = 0.9; 
    const double sig_a_2 = 0.01; 
    
    const double sig_s_1 = 1.3; 
    const double sig_s_2 = 2.7; 
    
    const int n_p = fem_quadrature.get_number_of_interpolation_points();
    const int n_dfem_int_pt = fem_quadrature.get_number_of_integration_points();
    const int n_xs_eval_pt = fem_quadrature.get_number_of_xs_point();
    
    std::vector<double> xs_quad(n_p,0.);
    std::vector<double> dfem_interp_points(n_p,0.);
    std::vector<double> integration_points(n_p,0.);
  
  // try{
    if( n_p  != n_dfem_int_pt)
      throw Dark_Arts_Exception(FEM , "For SLXS, number of integration points must equal the number of DFEM interpolation points");
      
    if( n_p  != n_xs_eval_pt)
      throw Dark_Arts_Exception(FEM , "For SLXS, number of xs evaluation points must equal the number of DFEM interpolation points");
      
    fem_quadrature.get_xs_eval_points(xs_quad);
    fem_quadrature.get_dfem_interpolation_point(dfem_interp_points);
    fem_quadrature.get_dfem_integration_points(integration_points);
    
    for(int i=0; i< n_p ; i++)
      if( fabs( dfem_interp_points[i] - xs_quad[i] ) > tol)
        throw Dark_Arts_Exception(FEM, "For SLXS, dfem interpolation points must be the same as the xs evaluation point quadrature");
        
    for(int i=0; i< n_p ; i++)
      if( fabs( integration_points[i] - xs_quad[i] ) > tol)
        throw Dark_Arts_Exception(FEM, "For SLXS, integration quadrature must be the same as the xs evaluation point quadrature");
  // }
  // catch(const Dark_Arts_Exception& da)
  // {
    // val = -1;
    // da.testing_message();
  // }
  
  
    double cv = 0., sig_a = 0. , sig_s = 0.;
    
    Eigen::VectorXd t_vec = Eigen::VectorXd::Zero(n_p);
    for(int i = 0; i < n_p ; i++)
      t_vec(i) = 3. + double(i) ; 
    
    const int total_cells = cell_data.get_total_number_of_cells();
    
    /// variables that we will test in every cell
    std::vector<double> calculated_cv_edge(2,0.);
    std::vector<double> calculated_sig_a_edge(2,0.);
    std::vector<double> calculated_sig_s_edge(2,0.);
    std::vector<double> calculated_cv_at_dfem_int_pts(n_dfem_int_pt,0.);
    std::vector<double> calcualted_sig_a_at_dfem_int_pts(n_dfem_int_pt, 0.);
    std::vector<double> calcualted_sig_s_at_dfem_int_pts(n_dfem_int_pt, 0.);

  // try
  // {
    for(int c = 0; c<total_cells ; c++)
    {
      materials.calculate_local_temp_and_position(c , t_vec);
      materials.get_cv_boundary(calculated_cv_edge);
      materials.get_sigma_a_boundary(0, calculated_sig_a_edge);
      materials.get_sigma_s_boundary(0, 0, calculated_sig_s_edge);
      materials.get_cv(calculated_cv_at_dfem_int_pts);
      materials.get_sigma_s(0,0,calcualted_sig_s_at_dfem_int_pts);  
      materials.get_sigma_a(0,calcualted_sig_a_at_dfem_int_pts);
      
      int mat_num = cell_data.get_cell_material_number(c);
      
      if(mat_num == 0)
      {
        cv = cv_1;
        sig_a = sig_a_1;
        sig_s = sig_s_1;
      }
      else if(mat_num == 1)
      {
        cv = cv_2;
        sig_a = sig_a_2;
        sig_s = sig_s_2;
      }
      
      for(int i=0 ; i < 2 ; i++)
      {
        if( fabs( calculated_cv_edge[i] - cv ) > tol )
          throw Dark_Arts_Exception( SUPPORT_OBJECT, "Incorrect calculation of cv at edges" );

        if(fabs(sig_s - calculated_sig_s_edge[i]) > tol )
          throw Dark_Arts_Exception(SUPPORT_OBJECT , "Difference in calculated edge sigma_s");

        if(fabs(calculated_sig_a_edge[i] - sig_a ) > tol )
          throw Dark_Arts_Exception(SUPPORT_OBJECT , "Difference in calculated edge sigma_a");
      }
      
      for(int i = 0; i < n_dfem_int_pt ; i++)
      {
        // std::cout <<"Cell_num: " << c << "sig_s expected: " << sig_s << " sig_s calculated: "  << calcualted_sig_s_at_dfem_int_pts[i] << std::endl;
        // if( fabs( sig_s - calcualted_sig_s_at_dfem_int_pts[i] ) > tol)
          // throw Dark_Arts_Exception(SUPPORT_OBJECT, "Discrepancy in calculating scattering cross section in space");
        
        // std::cout <<"Cell_num: " << c << "sig_a expected: " << sig_a << " sig_a calculated: "  << calcualted_sig_a_at_dfem_int_pts[i] << std::endl;               
        // if( fabs( sig_a - calcualted_sig_a_at_dfem_int_pts[i] ) > tol)
          // throw Dark_Arts_Exception(SUPPORT_OBJECT, "Discrepancy in calculating absorption cross section in space");
          
        if( fabs( cv - calculated_cv_at_dfem_int_pts[i] ) > tol)
          throw Dark_Arts_Exception(SUPPORT_OBJECT, "Discrepancy in calculating Cv in space");          
      }      
    }    
  }
  catch(const Dark_Arts_Exception& da)
  {
    val = -1;
    da.testing_message();
  }
  
  return val;
}
