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
  Materials materials( input_reader, fem_quadrature , cell_data, angular_quadrature.get_number_of_groups() , angular_quadrature.get_sum_w() );  
  
  // cubic dfem trial space, lobatto DFEM points, Interpolatory opacity treatment with P4 Equally_Spaced points
  // Exact matrix integration (via Gauss quadrature using (3 + 1 + 4) points),
  // Let T(x) = 1 + 2x + 3x^2 + 4x^3
  // sig_a(x) = 3.5  + 2.1 x^2 
  // sig_s(T) = 10/(T^3 + 0.1)
  // Cv(T) = 11/(T + 2)
  
  const double dx = 1.0;
  const double x_left = 0.;
  
  std::vector<double> sig_a_coeff(3,0.);  
  sig_a_coeff[0] = 3.5;
  sig_a_coeff[1] = 0.;
  sig_a_coeff[2] = 2.1;
  
  const int cell_num = 0;  
  const int n_p = fem_quadrature.get_number_of_interpolation_points();  
  Eigen::VectorXd t_vec(n_p);
  
  std::vector<double> dfem_interp_points(n_p,0.);
  fem_quadrature.get_dfem_interpolation_point(dfem_interp_points);
  
  std::vector<double> x_interp(n_p,0.);
  for(int i = 0 ; i < n_p ; i++)
  {
    x_interp[i] = x_left + dx/2.*(1. + dfem_interp_points[i]);
    t_vec(i) = 1. + 2. * x_interp[i] + 3.*x_interp[i]*x_interp[i] + 4.*pow(x_interp[i],3); 
  }
  
  /// calculate x and T at xs interpolation points
  const int n_xs = fem_quadrature.get_number_of_xs_point();
  std::vector<double> xs_quad(n_xs,0.);
  fem_quadrature.get_xs_eval_points(xs_quad);  
  std::vector<double> x_xs(n_xs,0.); 
  std::vector<double> t_xs(n_xs,0.);
  x_xs[0] = x_left;
  x_xs[1] = x_left + .25*dx;
  x_xs[2] = x_left + dx/2.;
  x_xs[3] = x_left + .75*dx;
  x_xs[4] = x_left + dx;
  /// trial space contains analytic T, T(x_xs) = T(x)
  for(int i=0; i < n_xs ; i++)
    t_xs[i] = 1. + 2. * x_xs[i] + 3.*x_xs[i]*x_xs[i] + 4.*pow(x_xs[i],3); 
    
  std::vector<double> sig_a(n_xs,0.);
  std::vector<double> sig_s(n_xs,0.);
  std::vector<double> cv(n_xs,0.);
  for(int i=0; i < n_xs ; i++)
  {
    sig_a[i] = sig_a_coeff[0] + sig_a_coeff[1] * x_xs[i] + sig_a_coeff[2]*x_xs[i]*x_xs[i];
    /// 10/(T^3 + 0.1)
    sig_s[i] = 10./(pow(t_xs[i] , 3) + 0.1);
    cv[i] = 11./(t_xs[i] + 2.);
  }
  
  /// get xs polynomials at dfem integration points
  std::vector<double> xs_at_dfem;
  const int n_int = fem_quadrature.get_number_of_integration_points();
  fem_quadrature.get_xs_at_dfem_integration_points(xs_at_dfem);
  std::vector<double> expected_sig_a(n_int,0.);
  std::vector<double> expected_sig_s(n_int,0.);
  std::vector<double> expected_cv(n_int,0.);
  int cnt = 0;
  for(int xs=0; xs< n_xs ; xs++)
  {
    for(int int_pt = 0 ; int_pt < n_int ; int_pt++)
    {  
      expected_sig_a[int_pt] += sig_a[xs]*xs_at_dfem[cnt];
      expected_sig_s[int_pt] += sig_s[xs]*xs_at_dfem[cnt];
      expected_cv[int_pt] += cv[xs]*xs_at_dfem[cnt];
      cnt++;
    }
  }
  
  materials.calculate_local_temp_and_position(cell_num , t_vec);
  
  // /// test edge values
  // std::vector<double> calculated_cv_edge(2,0.);
  // materials.get_cv_boundary(calculated_cv_edge);
  // try{
    // for(int i=0 ; i < 2 ; i++)
    // {
      // std::cout << "Calculated cv: " << calculated_cv_edge[i] << " Expected Cv: " << cv << std::endl;
      // if( fabs( calculated_cv_edge[i] - cv ) > tol )
        // throw Dark_Arts_Exception( SUPPORT_OBJECT, "Incorrect calculation of cv at edges" );
    // }
  // }
  // catch(const Dark_Arts_Exception& da)
  // {
    // val = -1;
    // da.testing_message();
  // }
  
  // std::vector<double> expected_sig_s_edge(2,0.);
  // std::vector<double> expected_sig_a_edge(2,0.);
  // std::vector<double> calculated_sig_a_edge(2,0.);
  // std::vector<double> calculated_sig_s_edge(2,0.);
  
  // expected_sig_s_edge[0] = sig_s_coeff[0] + x_left*sig_s_coeff[1];
  // expected_sig_s_edge[1] = sig_s_coeff[0] + x_right*sig_s_coeff[1];
  
  // expected_sig_a_edge[0] = sig_a/(t_left*t_left);
  // expected_sig_a_edge[1] = sig_a/(t_right*t_right);
  // materials.get_sigma_a_boundary(0, calculated_sig_a_edge);
  // materials.get_sigma_s_boundary(0, 0, calculated_sig_s_edge);
  
  // try{
    // for(int i=0; i<2 ; i++)
      // if(fabs(expected_sig_s_edge[i] - calculated_sig_s_edge[i]) > tol )
        // throw Dark_Arts_Exception(SUPPORT_OBJECT , "Difference in calculated edge sigma_s");
        
    // for(int i=0; i<2 ; i++)
      // if(fabs(expected_sig_a_edge[i] - calculated_sig_a_edge[i]) > tol )
        // throw Dark_Arts_Exception(SUPPORT_OBJECT , "Difference in calculated edge sigma_a");
  // }
  // catch(const Dark_Arts_Exception& da)
  // {
    // val = -1;
    // da.testing_message();
  // }
  
  
  /// get the values from Materials
  std::vector<double> calculated_cv_at_dfem_int_pts(n_int,0.);
  std::vector<double> calcualted_sig_a_at_dfem_int_pts(n_int, 0.);
  std::vector<double> calcualted_sig_s_at_dfem_int_pts(n_int, 0.);
  try{
    materials.get_cv(calculated_cv_at_dfem_int_pts);
    materials.get_sigma_s(0,0,calcualted_sig_s_at_dfem_int_pts);  
    materials.get_sigma_a(0,calcualted_sig_a_at_dfem_int_pts);
    
    for(int i = 0; i < n_int ; i++)
    {
      std::cout << " Expected sig_s: " << expected_sig_s[i] << " Calculated sig_s: " << calcualted_sig_s_at_dfem_int_pts[i] << std::endl;
      if( fabs( expected_sig_s[i] - calcualted_sig_s_at_dfem_int_pts[i] ) > tol)
        throw Dark_Arts_Exception(SUPPORT_OBJECT, "Discrepancy in calculating scattering cross section in space");
        
      if( fabs( expected_sig_a[i] - calcualted_sig_a_at_dfem_int_pts[i] ) > tol)
        throw Dark_Arts_Exception(SUPPORT_OBJECT, "Discrepancy in calculating absorption cross section in space");
        
      if( fabs( expected_cv[i] - calculated_cv_at_dfem_int_pts[i] ) > tol)
        throw Dark_Arts_Exception(SUPPORT_OBJECT, "Discrepancy in calculating Cv in space");
        
    }
  }
  catch(const Dark_Arts_Exception& da)
  {
    val = -1;
    da.testing_message();
  }
  
  // Return 0 if tests passed, something else if failing
  return val;
}
