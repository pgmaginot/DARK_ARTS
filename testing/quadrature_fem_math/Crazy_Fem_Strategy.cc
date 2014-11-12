#include <iostream>
#include <vector>
#include <math.h>
#include <string>

#include "Input_Reader.h"
#include "Fem_Quadrature.h"
#include "Quadrule_New.h"
#include "Dark_Arts_Exception.h"

int main(int argc, char** argv)
{
  int val = 0;
  const double tol = 1.0E-6;
  
  /// get all the objects necessary to create an Angular_Quadrature object
  Input_Reader input_reader;  
  try{
    input_reader.read_xml(argv[1]);  
  }
  catch(const Dark_Arts_Exception& da_exception)
  {
    da_exception.testing_message();
    val = -1;
  }
  Quadrule_New quad_fun;   
  Fem_Quadrature fem_quadrature( input_reader , quad_fun);   

  const int n_dfem_pt = 4;
  const int n_xs_pt = 5;
  const int n_int_pt = 8;
  
  std::vector<double> dfem_pts(n_dfem_pt,0.);
  std::vector<double> xs_pts(n_xs_pt,0.);
  std::vector<double> integration_pts(n_int_pt,0.);
  std::vector<double> dfem_at_integration_pts(n_dfem_pt*n_int_pt,0.);
  std::vector<double> xs_at_integration_pts(n_xs_pt*n_int_pt,0.);
  std::vector<double> dfem_deriv_at_integration_pts(n_dfem_pt*n_int_pt,0.);
  /// DFEM interpolation pts, XS interpolation pts , integration quadrature pts
  dfem_pts[0] =   -1.00000000 ;
  dfem_pts[1] =   -0.44721360 ;
  dfem_pts[2] =    0.44721360 ;
  dfem_pts[3] =    1.00000000 ;
  xs_pts[0] =   -1.00000000 ;
  xs_pts[1] =   -0.50000000 ;
  xs_pts[2] =    0.00000000 ;
  xs_pts[3] =    0.50000000 ;
  xs_pts[4] =    1.00000000 ;
  integration_pts[0] =   -0.96028986 ;
  integration_pts[1] =   -0.79666648 ;
  integration_pts[2] =   -0.52553241 ;
  integration_pts[3] =   -0.18343464 ;
  integration_pts[4] =    0.18343464 ;
  integration_pts[5] =    0.52553241 ;
  integration_pts[6] =    0.79666648 ;
  integration_pts[7] =    0.96028986 ;
  
  /// Matlab DFEM functions at integration pts
  dfem_at_integration_pts[0] =   0.88477267 ;
  dfem_at_integration_pts[1] =   0.48810653 ;
  dfem_at_integration_pts[2] =   0.07263852 ;
  dfem_at_integration_pts[3] =  -0.12304150 ;
  dfem_at_integration_pts[4] =  -0.08489816 ;
  dfem_at_integration_pts[5] =   0.02259187 ;
  dfem_at_integration_pts[6] =   0.05524031 ;
  dfem_at_integration_pts[7] =   0.01792309 ;
  dfem_at_integration_pts[8] =   0.15312152 ;
  dfem_at_integration_pts[9] =   0.63506763 ;
  dfem_at_integration_pts[10] =   0.98399404 ;
  dfem_at_integration_pts[11] =   0.85170155 ;
  dfem_at_integration_pts[12] =   0.35623812 ;
  dfem_at_integration_pts[13] =  -0.07922443 ;
  dfem_at_integration_pts[14] =  -0.17841448 ;
  dfem_at_integration_pts[15] =  -0.05581728 ;
  dfem_at_integration_pts[16] =  -0.05581728 ;
  dfem_at_integration_pts[17] =  -0.17841448 ;
  dfem_at_integration_pts[18] =  -0.07922443 ;
  dfem_at_integration_pts[19] =   0.35623812 ;
  dfem_at_integration_pts[20] =   0.85170155 ;
  dfem_at_integration_pts[21] =   0.98399404 ;
  dfem_at_integration_pts[22] =   0.63506763 ;
  dfem_at_integration_pts[23] =   0.15312152 ;
  dfem_at_integration_pts[24] =   0.01792309 ;
  dfem_at_integration_pts[25] =   0.05524031 ;
  dfem_at_integration_pts[26] =   0.02259187 ;
  dfem_at_integration_pts[27] =  -0.08489816 ;
  dfem_at_integration_pts[28] =  -0.12304150 ;
  dfem_at_integration_pts[29] =   0.07263852 ;
  dfem_at_integration_pts[30] =   0.48810653 ;
  dfem_at_integration_pts[31] =   0.88477267 ;
  
  dfem_deriv_at_integration_pts[0] =  -2.80440596 ;
  dfem_deriv_at_integration_pts[1] =  -2.06085336 ;
  dfem_deriv_at_integration_pts[2] =  -1.04976110 ;
  dfem_deriv_at_integration_pts[3] =  -0.16738381 ;
  dfem_deriv_at_integration_pts[4] =   0.29120280 ;
  dfem_deriv_at_integration_pts[5] =   0.26406992 ;
  dfem_deriv_at_integration_pts[6] =  -0.06918717 ;
  dfem_deriv_at_integration_pts[7] =  -0.40368132 ;
  dfem_deriv_at_integration_pts[8] =   3.66907895 ;
  dfem_deriv_at_integration_pts[9] =   2.25925682 ;
  dfem_deriv_at_integration_pts[10] =   0.41731096 ;
  dfem_deriv_at_integration_pts[11] =  -1.02717453 ;
  dfem_deriv_at_integration_pts[12] =  -1.48576114 ;
  dfem_deriv_at_integration_pts[13] =  -0.89652006 ;
  dfem_deriv_at_integration_pts[14] =   0.26759063 ;
  dfem_deriv_at_integration_pts[15] =   1.26835431 ;
  dfem_deriv_at_integration_pts[16] =  -1.26835431 ;
  dfem_deriv_at_integration_pts[17] =  -0.26759063 ;
  dfem_deriv_at_integration_pts[18] =   0.89652006 ;
  dfem_deriv_at_integration_pts[19] =   1.48576114 ;
  dfem_deriv_at_integration_pts[20] =   1.02717453 ;
  dfem_deriv_at_integration_pts[21] =  -0.41731096 ;
  dfem_deriv_at_integration_pts[22] =  -2.25925682 ;
  dfem_deriv_at_integration_pts[23] =  -3.66907895 ;
  dfem_deriv_at_integration_pts[24] =   0.40368132 ;
  dfem_deriv_at_integration_pts[25] =   0.06918717 ;
  dfem_deriv_at_integration_pts[26] =  -0.26406992 ;
  dfem_deriv_at_integration_pts[27] =  -0.29120280 ;
  dfem_deriv_at_integration_pts[28] =   0.16738381 ;
  dfem_deriv_at_integration_pts[29] =   1.04976110 ;
  dfem_deriv_at_integration_pts[30] =   2.06085336 ;
  dfem_deriv_at_integration_pts[31] =   2.80440596 ;
  
  xs_at_integration_pts[0] =   0.84353255 ;
  xs_at_integration_pts[1] =   0.36707052 ;
  xs_at_integration_pts[2] =   0.01399493 ;
  xs_at_integration_pts[3] =  -0.03131084 ;
  xs_at_integration_pts[4] =   0.02160436 ;
  xs_at_integration_pts[5] =  -0.00435267 ;
  xs_at_integration_pts[6] =  -0.04154235 ;
  xs_at_integration_pts[7] =  -0.01708768 ;
  xs_at_integration_pts[8] =   0.29109309 ;
  xs_at_integration_pts[9] =   1.00635222 ;
  xs_at_integration_pts[10] =   1.04026891 ;
  xs_at_integration_pts[11] =   0.32305934 ;
  xs_at_integration_pts[12] =  -0.14964035 ;
  xs_at_integration_pts[13] =   0.02589930 ;
  xs_at_integration_pts[14] =   0.23024500 ;
  xs_at_integration_pts[15] =   0.09175384 ;
  xs_at_integration_pts[16] =  -0.20929180 ;
  xs_at_integration_pts[17] =  -0.56212539 ;
  xs_at_integration_pts[18] =  -0.07581047 ;
  xs_at_integration_pts[19] =   0.83628748 ;
  xs_at_integration_pts[20] =   0.83628748 ;
  xs_at_integration_pts[21] =  -0.07581047 ;
  xs_at_integration_pts[22] =  -0.56212539 ;
  xs_at_integration_pts[23] =  -0.20929180 ;
  xs_at_integration_pts[24] =   0.09175384 ;
  xs_at_integration_pts[25] =   0.23024500 ;
  xs_at_integration_pts[26] =   0.02589930 ;
  xs_at_integration_pts[27] =  -0.14964035 ;
  xs_at_integration_pts[28] =   0.32305934 ;
  xs_at_integration_pts[29] =   1.04026891 ;
  xs_at_integration_pts[30] =   1.00635222 ;
  xs_at_integration_pts[31] =   0.29109309 ;
  xs_at_integration_pts[32] =  -0.01708768 ;
  xs_at_integration_pts[33] =  -0.04154235 ;
  xs_at_integration_pts[34] =  -0.00435267 ;
  xs_at_integration_pts[35] =   0.02160436 ;
  xs_at_integration_pts[36] =  -0.03131084 ;
  xs_at_integration_pts[37] =   0.01399493 ;
  xs_at_integration_pts[38] =   0.36707052 ;
  xs_at_integration_pts[39] =   0.84353255 ;

  /// make sure Fem_Quadrature has the correct number of points
  try{
    if( (n_xs_pt - fem_quadrature.get_number_of_xs_point() ) != 0 )
      throw Dark_Arts_Exception(FEM , "Expected a different number of xs interpolation points");
  }
  catch(const Dark_Arts_Exception& da)
  {
    da.testing_message();
    val = -1;
  }
  try{
    if( (n_int_pt - fem_quadrature.get_number_of_integration_points() ) != 0 )
      throw Dark_Arts_Exception(FEM , "Expected a different number of xs interpolation points");
  }
  catch(const Dark_Arts_Exception& da)
  {
    da.testing_message();
    val = -1;
  }
  try{
    if( (n_dfem_pt - fem_quadrature.get_number_of_interpolation_points() ) != 0 )
      throw Dark_Arts_Exception(FEM , "Expected a different number of xs interpolation points");
  }
  catch(const Dark_Arts_Exception& da)
  {
    da.testing_message();
    val = -1;
  }
  
  std::vector<double> da_dfem_pt;
  std::vector<double> da_xs_pt;
  std::vector<double> da_int_pt;
  std::vector<double> da_dfem_at_int_pt;
  std::vector<double> da_xs_at_int_pt;
  std::vector<double> da_dfem_deriv_at_int_pt;
  
  fem_quadrature.get_dfem_interpolation_point(da_dfem_pt);
  fem_quadrature.get_xs_eval_points(da_xs_pt);
  fem_quadrature.get_dfem_integration_points(da_int_pt);
  
  fem_quadrature.get_dfem_at_integration_points(da_dfem_at_int_pt);
  fem_quadrature.get_xs_at_dfem_integration_points(da_xs_at_int_pt);
  fem_quadrature.get_dfem_derivatives_at_integration_points(da_dfem_deriv_at_int_pt);
  
  try{
    for(int i=0; i < n_dfem_pt ; i++)
    {
      if(fabs(da_dfem_pt[i] - dfem_pts[i] ) > tol)
        throw Dark_Arts_Exception(FEM , "DFEM interpolation points not correct" );
    }
  }
  catch(const Dark_Arts_Exception& da)
  {
    da.testing_message();
    val = -1;
  }
  
  try{
    for(int i=0; i < n_xs_pt ; i++)
    {
      if(fabs(da_xs_pt[i] - xs_pts[i] ) > tol)
        throw Dark_Arts_Exception(FEM , "XS interpolation points not correct" );
    }
  }
  catch(const Dark_Arts_Exception& da)
  {
    da.testing_message();
    val = -1;
  }
  
  try{
    for(int i=0; i < n_int_pt ; i++)
    {
      if(fabs(da_int_pt[i] - integration_pts[i] ) > tol)
        throw Dark_Arts_Exception(FEM , "Integration points not correct" );
    }
  }
  catch(const Dark_Arts_Exception& da)
  {
    da.testing_message();
    val = -1;
  }
  
  try{
    if( da_dfem_at_int_pt.size() != n_int_pt*n_dfem_pt)
    {
      std::stringstream err;
      err << "DFEM function evaluations length is: " << da_dfem_deriv_at_int_pt.size() << " Excepted to be of length: " <<n_int_pt*n_dfem_pt;
      throw Dark_Arts_Exception(FEM , err.str() );
    }
    for(int i=0; i< (n_int_pt*n_dfem_pt) ; i++)
    {
      if(fabs( da_dfem_at_int_pt[i] - dfem_at_integration_pts[i] ) > tol )
        throw Dark_Arts_Exception(FEM , "Basis functions at integration points not correct" );
    }      
  }
  catch(const Dark_Arts_Exception& da)
  {
    da.testing_message();
    val = -1;
  }
  
  try{
    if( da_dfem_deriv_at_int_pt.size() != n_int_pt*n_dfem_pt)
    {
      std::stringstream err;
      err << "Derivative evaluations length is: " << da_dfem_deriv_at_int_pt.size() << " Expected to be of length: " <<n_int_pt*n_dfem_pt;
      throw Dark_Arts_Exception(FEM , err.str() );
    }
    for(int i=0; i< (n_int_pt*n_dfem_pt) ; i++)
    {
      std::cout << "Expected derivative: " << dfem_deriv_at_integration_pts[i] << " calcualted derivative: " << da_dfem_deriv_at_int_pt[i] << std::endl;
      if(fabs( da_dfem_deriv_at_int_pt[i] - dfem_deriv_at_integration_pts[i] ) > tol )
        throw Dark_Arts_Exception(FEM , "Basis derivatives at integration points not correct" );
    }      
  }
  catch(const Dark_Arts_Exception& da)
  {
    da.testing_message();
    val = -1;
  }
  
  try{
    if( da_xs_at_int_pt.size() != n_int_pt*n_xs_pt)
    {
      std::stringstream err;
      err << "XS function evaluations length is: " << da_xs_at_int_pt.size() << " Expected to be of length: " << n_int_pt*n_xs_pt;
      throw Dark_Arts_Exception(FEM , err.str() );
    }
    for(int i=0; i < (n_int_pt*n_xs_pt) ; i++)
    {
      if(fabs( da_xs_at_int_pt[i] - xs_at_integration_pts[i] ) > tol )
        throw Dark_Arts_Exception(FEM , "XS functions at integration points not correct" );
    }      
  }
  catch(const Dark_Arts_Exception& da)
  {
    da.testing_message();
    val = -1;
  }
  
  
  // Return 0 if tests passed, something else if failing
  return val;
}


