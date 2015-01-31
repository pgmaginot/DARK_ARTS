#include <iostream>
#include <vector>
#include <math.h>
#include <string>
#include <iomanip>

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
  const double tol = 1.0E-8;
  
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
  Materials materials( input_reader, fem_quadrature , cell_data, angular_quadrature);  
  
  /// verify that we get expected (constant values for material properties evaluations);
  /// constant in space
  const double cv = 1.5;
  /// make temperature linear in space
  const double t_left = 1.;
  const double t_right = 3.;
  const double x_left = t_left;
  const double x_right = t_right;
  
  const double sig_a = 2.2; /// sig_{a} = sig_a/(T^2 + 0.)

  std::vector<double> sig_s_coeff(2,0.); ///sig_s = sig_s_coeff[0] + x * sig_s_coeff[1] 
  sig_s_coeff[0] = 1.0;
  sig_s_coeff[1] = 1.5;
  
  const int cell_num = 0;
  
  const int n_p = fem_quadrature.get_number_of_interpolation_points();
  
  Eigen::VectorXd t_vec(n_p);
  t_vec(0) = t_left;
  t_vec(1) = (t_left + t_right)/2.;
  t_vec(2) = t_right;
  
  /// calculate local temperature and position in Materials object
  try{
    materials.calculate_local_temp_and_position(cell_num , t_vec);
  }
  catch(const Dark_Arts_Exception& da)
  {
    val = -1;
    da.testing_message();
  }
  
  /// Generate the data we expect to get from Materials automatically
  const int n_xs = fem_quadrature.get_number_of_xs_point();
  std::vector<double> xs_quad;
  std::vector<double> xs_wt;
  fem_quadrature.get_xs_eval_points(xs_quad);
  fem_quadrature.get_xs_eval_weights(xs_wt);
  
  std::vector<double> x_pt(n_xs,0.); /// same as t_pt
  for(int i=0; i < n_xs ; i++)
  {
    x_pt[i] = x_left + (x_right-x_left)/2.*( 1. + xs_quad[i] );    
  }
  
  /// test edge values
  std::vector<double> calculated_cv_edge(2,0.);
  materials.get_cv_boundary(calculated_cv_edge);
  try{
    for(int i=0 ; i < 2 ; i++)
    {
      std::cout << "Calculated cv: " << calculated_cv_edge[i] << " Expected Cv: " << cv << std::endl;
      if( fabs( (calculated_cv_edge[i] - cv )/cv) > tol )
        throw Dark_Arts_Exception( SUPPORT_OBJECT, "Incorrect calculation of cv at edges" );
    }
  }
  catch(const Dark_Arts_Exception& da)
  {
    val = -1;
    da.testing_message();
  }
  
  std::vector<double> expected_sig_s_edge(2,0.);
  std::vector<double> expected_sig_a_edge(2,0.);
  std::vector<double> calculated_sig_a_edge(2,0.);
  std::vector<double> calculated_sig_s_edge(2,0.);
  
  expected_sig_s_edge[0] = sig_s_coeff[0] + x_left*sig_s_coeff[1];
  expected_sig_s_edge[1] = sig_s_coeff[0] + x_right*sig_s_coeff[1];
  
  expected_sig_a_edge[0] = sig_a/(t_left*t_left);
  expected_sig_a_edge[1] = sig_a/(t_right*t_right);
  materials.get_sigma_a_boundary(0, calculated_sig_a_edge);
  materials.get_sigma_s_boundary(0, 0, calculated_sig_s_edge);
  
  try{
    for(int i=0; i<2 ; i++)
      if(fabs((expected_sig_s_edge[i] - calculated_sig_s_edge[i])/expected_sig_s_edge[i] ) > tol )
        throw Dark_Arts_Exception(SUPPORT_OBJECT , "Difference in calculated edge sigma_s");
        
    for(int i=0; i<2 ; i++)
      if(fabs((expected_sig_a_edge[i] - calculated_sig_a_edge[i])/expected_sig_a_edge[i]) > tol )
        throw Dark_Arts_Exception(SUPPORT_OBJECT , "Difference in calculated edge sigma_a");
  }
  catch(const Dark_Arts_Exception& da)
  {
    val = -1;
    da.testing_message();
  }
  
  
  const int n_dfem_int_pt = fem_quadrature.get_number_of_integration_points();
  std::vector<double> evaluated_sig_a_at_dfem_int_pts(n_dfem_int_pt,0.);
  std::vector<double> evaluated_sig_s_at_dfem_int_pts(n_dfem_int_pt,0.);
  std::vector<double> evaluated_cv_at_dfem_int_pts(n_dfem_int_pt,cv);
  std::vector<double> dfem_int_pts;
  fem_quadrature.get_dfem_integration_points(dfem_int_pts);
  std::vector<double> x_int(n_dfem_int_pt,0.);
  
  /// quadratic moment preservation will get exact (linear in space) scattering cross section
  for(int i=0; i< n_dfem_int_pt ; i++)
  {
    x_int[i] = x_left + (x_right - x_left)/2.*( 1. + dfem_int_pts[i] );
    evaluated_sig_s_at_dfem_int_pts[i] = sig_s_coeff[0] + sig_s_coeff[1]*x_int[i];
  }
  
  /// have to do some pretty serious work to calculate Legendre moments of sig_a
  std::vector<double> raw_sig_a_evals(n_xs,0.);
  for(int i=0; i < n_xs ; i++)
  {
    raw_sig_a_evals[i] = sig_a/(x_pt[i]*x_pt[i]);  
  }
  double sig_a_0 = 0.;
  double sig_a_1 = 0.;
  double sig_a_2 = 0.;
  for(int i=0; i< n_xs ; i++)
  {
    /// 0-th legendre polynomial is 1
    sig_a_0 += raw_sig_a_evals[i]*xs_wt[i]*(1./2.);
    /// 1st legendre polynomial is x
    sig_a_1 += raw_sig_a_evals[i]*xs_wt[i]*(3./2.)*xs_quad[i];
    /// 2nd legendre polynomial is (3x^2 - 1)/2
    sig_a_2 += raw_sig_a_evals[i]*xs_wt[i]*(5./2.)*(3.*xs_quad[i]*xs_quad[i]  - 1.)/2.;
  }
  for(int i=0; i < n_dfem_int_pt ; i++)
  {
    evaluated_sig_a_at_dfem_int_pts[i] = sig_a_0 + dfem_int_pts[i]*sig_a_1 + (3.*dfem_int_pts[i]*dfem_int_pts[i]  - 1.)/2.*sig_a_2;
  }
  
  
  /** these Materials functions return cross section evaluated at the DFEM integration points
    void get_sigma_a(const int grp, std::vector<double>& sig_a);
    void get_sigma_s(const int grp, const int l_mom, std::vector<double>& sig_s);
    void get_cv(std::vector<double>& cv);
  */  
  /// get the values from Materials
  std::vector<double> calculated_cv_at_dfem_int_pts(n_dfem_int_pt,0.);
  std::vector<double> calcualted_sig_a_at_dfem_int_pts(n_dfem_int_pt, 0.);
  std::vector<double> calcualted_sig_s_at_dfem_int_pts(n_dfem_int_pt, 0.);
  try{
    materials.get_cv(calculated_cv_at_dfem_int_pts);
    materials.get_sigma_s(0,0,calcualted_sig_s_at_dfem_int_pts);  
    materials.get_sigma_a(0,calcualted_sig_a_at_dfem_int_pts);
    
    for(int i = 0; i < n_dfem_int_pt ; i++)
    {
      std::cout << "x: " << x_int[i] << " Expected sig_s: " << evaluated_sig_s_at_dfem_int_pts[i] << " Calculated sig_s: " << calcualted_sig_s_at_dfem_int_pts[i] << std::endl;
      if( fabs( (evaluated_sig_s_at_dfem_int_pts[i] - calcualted_sig_s_at_dfem_int_pts[i] )/evaluated_sig_s_at_dfem_int_pts[i]) > tol)
        throw Dark_Arts_Exception(SUPPORT_OBJECT, "Discrepancy in calculating scattering cross section in space");
        
      if( fabs( ( evaluated_sig_a_at_dfem_int_pts[i] - calcualted_sig_a_at_dfem_int_pts[i] )/evaluated_sig_a_at_dfem_int_pts[i]) > tol)
        throw Dark_Arts_Exception(SUPPORT_OBJECT, "Discrepancy in calculating absorption cross section in space");
        
      if( fabs( (evaluated_cv_at_dfem_int_pts[i] - evaluated_cv_at_dfem_int_pts[i] )/evaluated_cv_at_dfem_int_pts[i]) > tol)
        throw Dark_Arts_Exception(SUPPORT_OBJECT, "Discrepancy in calculating Cv in space");
        
    }
  }
  catch(const Dark_Arts_Exception& da)
  {
    val = -1;
    da.testing_message();
  }
  
  Eigen::VectorXd planck = Eigen::VectorXd::Zero(n_p);
  Eigen::MatrixXd d_planck = Eigen::MatrixXd::Zero(n_p,n_p);
  
  materials.get_grey_planck(t_vec , planck);
  materials.get_grey_planck_derivative(t_vec, d_planck);
  
  try{
    for(int i=0; i < n_p ; i++)
    {
      if( fabs( planck(i) - pow(t_vec(i) , 4)/angular_quadrature.get_sum_w() ) > tol)
        throw Dark_Arts_Exception(SUPPORT_OBJECT, "Problem calculating vector planck");
      
      if( fabs( materials.get_grey_planck( t_vec(i) ) - pow(t_vec(i) , 4)/angular_quadrature.get_sum_w() ) > tol)
        throw Dark_Arts_Exception(SUPPORT_OBJECT, "Problem calculating pointwise planck");
      
      if( fabs( d_planck(i,i) - 4.*pow(t_vec(i) , 3)/angular_quadrature.get_sum_w() ) > tol)
        throw Dark_Arts_Exception(SUPPORT_OBJECT, "Problem calculating matrix planck derivative");
              
      for(int j=0; j< n_p ; j++)
      {
        if( i==j)
          continue;
        else
          if( fabs( d_planck(i,j) ) > tol )
          {
            std::cout << std::setprecision(15) << std::scientific << d_planck << std::endl;
            throw Dark_Arts_Exception( SUPPORT_OBJECT , "Off-diagonal element of D matrix expected to be zero");            
          }
      }        
    }
  }
  catch(const Dark_Arts_Exception& da)
  {
    val = -1;
    da.testing_message();
  }
  
  try{
    if( fabs( materials.get_c() - 1. ) > tol)
      throw Dark_Arts_Exception( SUPPORT_OBJECT, "Expected c to be unity");
      
    if( fabs( materials.get_cell_width() - (x_right - x_left) ) > tol)
      throw Dark_Arts_Exception(SUPPORT_OBJECT , "Materials returning incorrect cell width");
  }
  catch(const Dark_Arts_Exception& da)
  {
    val = -1;
    da.testing_message();
  }
  
  // Return 0 if tests passed, something else if failing
  return val;
}
