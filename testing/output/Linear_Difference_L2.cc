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
#include "Time_Data.h"
#include "Intensity_Moment_Data.h"
#include "Temperature_Data.h"
#include "MMS_Intensity.h"
#include "MMS_Temperature.h"
#include "L2_Error_Calculator.h"

#include "Dark_Arts_Exception.h"

/**
  Goal of this unit test is to verify that L2 error estimate is zero when the analytic solution is within the DFEM trial space
*/ 

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
  
  Quadrule_New quad_fun;  
  Fem_Quadrature fem_quadrature( input_reader , quad_fun);  
  Cell_Data cell_data( input_reader );  
  Angular_Quadrature angular_quadrature( input_reader , quad_fun );    
  
  /// Create a Materials object that contains all opacity, heat capacity, and source objects
  Materials materials( input_reader, fem_quadrature , cell_data, angular_quadrature );  
  Time_Data time_data(input_reader);
  
  Temperature_Data t_old(fem_quadrature, input_reader, cell_data);  
  Intensity_Data i_old(cell_data, angular_quadrature, fem_quadrature, materials,input_reader);    
  Intensity_Moment_Data phi(cell_data,angular_quadrature,fem_quadrature,i_old);  
  
  MMS_Intensity mms_intensity(input_reader,angular_quadrature);
  MMS_Temperature mms_temperature(input_reader);
  
  try{ 
    L2_Error_Calculator l2_error_calculator(angular_quadrature,fem_quadrature, cell_data, input_reader);
    const double time_eval = time_data.get_t_start();
    
    double temperature_err = l2_error_calculator.calculate_l2_error(time_eval , t_old);
    double phi_err = l2_error_calculator.calculate_l2_error(time_eval , phi);
    
    double temperature_avg_err = l2_error_calculator.calculate_cell_avg_error(time_eval,t_old);
    double phi_avg_err = l2_error_calculator.calculate_cell_avg_error(time_eval , phi);
    
    const int n_int_pts = fem_quadrature.get_number_of_integration_points();
    const int n_xs_pts = fem_quadrature.get_number_of_xs_point();
    /// integrating a 4th order polynomial, need at least 3 gauss points
    std::cout << "N xs pt: " << n_xs_pts << "\n n integration_pts: " << n_int_pts << std::endl;
    
    if(n_int_pts < 3)
      throw Dark_Arts_Exception(SUPPORT_OBJECT , "Not enough gauss points to get exact integration");
      
    if(n_xs_pts < 3)
      throw Dark_Arts_Exception(SUPPORT_OBJECT, "Not enough xs evaluation points to get exact integration");
  
    /// use the xs quadrature to verify what we have calculated with the integration quadrature
    const int n_cell = cell_data.get_total_number_of_cells();
    std::vector<double> xs_quad_weight;
    std::vector<double> xs_quad_pt;
    std::vector<double> dfem_at_xs_quad;
    fem_quadrature.get_xs_eval_weights(xs_quad_weight);
    fem_quadrature.get_xs_eval_points(xs_quad_pt);
    fem_quadrature.get_dfem_at_xs_eval_points(dfem_at_xs_quad);
    
    
    Eigen::VectorXd num_soln = Eigen::VectorXd::Zero(fem_quadrature.get_number_of_interpolation_points() );
    std::vector<double> evaluated_numerical(n_xs_pts,0.);
    double err = 0.;
    double xL = 0.; 
    double dx = 0.;
      
    /// Angle integrated intensity (phi) errors    
    double expected_phi_err = 0.;
    double expected_phi_avg_err = 0.;
    double numeric_phi_A = 0.;   
    double analytic_phi_A = 0.;
    for(int c= 0 ; c < n_cell ; c++)
    {
      phi.get_cell_angle_integrated_intensity(c,0,0,num_soln);
      fem_quadrature.evaluate_variable_at_quadrature_pts(num_soln, dfem_at_xs_quad , evaluated_numerical );
      xL = cell_data.get_cell_left_edge(c);
      dx = cell_data.get_cell_width(c);
      numeric_phi_A = 0.;
      analytic_phi_A = 0.;
      err = 0.;
      for(int q = 0; q < n_xs_pts ; q++)
      {
        double x_eval  = xL + dx/2.*(1. + xs_quad_pt[q]);
        double phi_analytic = mms_intensity.get_mms_phi(x_eval , time_eval);
        numeric_phi_A += xs_quad_weight[q]*evaluated_numerical[q];
        analytic_phi_A += xs_quad_weight[q]*phi_analytic;
        
        err += xs_quad_weight[q]*(phi_analytic - evaluated_numerical[q] )*(phi_analytic - evaluated_numerical[q] );
      }
      expected_phi_err += dx/2.*err;
      
      numeric_phi_A /= 2.;
      analytic_phi_A /= 2.;
      
      expected_phi_avg_err += dx*(numeric_phi_A - analytic_phi_A)*(numeric_phi_A - analytic_phi_A);
    }
    expected_phi_err = sqrt(expected_phi_err);
    expected_phi_avg_err = sqrt(expected_phi_avg_err);
    
    std::cout << "Expected_phi_err:  " << expected_phi_err << std::endl;
    std::cout << "Expected phi_avg err: " << expected_phi_avg_err << std::endl;    
    
    std::cout << "Calculated Phi L2 err: " << phi_err << std::endl;
    std::cout << "Calculated Phi avg err: " << phi_avg_err << std::endl;
    
    /// Temperature errors
    
    double expected_t_err = 0.;
    double expected_t_avg_err = 0.;
    for(int c= 0 ; c < n_cell ; c++)
    {
      t_old.get_cell_temperature(c,num_soln);
      fem_quadrature.evaluate_variable_at_quadrature_pts(num_soln, dfem_at_xs_quad , evaluated_numerical );
      xL = cell_data.get_cell_left_edge(c);
      dx = cell_data.get_cell_width(c);
      
      err = 0.;
      double numeric_t_avg = 0.;
      double analytic_t_avg = 0.;    
      for(int q = 0; q < n_xs_pts ; q++)
      {
        double x_eval  = xL + dx/2.*(1. + xs_quad_pt[q]);
        double t_analytic = mms_temperature.get_mms_temperature(x_eval , time_eval);
        numeric_t_avg += xs_quad_weight[q]*evaluated_numerical[q];
        analytic_t_avg += xs_quad_weight[q]*t_analytic;
        
        err += xs_quad_weight[q]*(t_analytic - evaluated_numerical[q] )*(t_analytic - evaluated_numerical[q] );
      }
      expected_t_err += dx/2.*err;
      
      numeric_t_avg /= 2.;
      analytic_t_avg /= 2.;
      
      expected_t_avg_err += dx*(numeric_t_avg - analytic_t_avg)*(numeric_t_avg - analytic_t_avg);
    }
    expected_t_err = sqrt(expected_t_err);
    expected_t_avg_err = sqrt(expected_t_avg_err);
    
    std::cout << "Expected_t_err:  " << expected_t_err << std::endl;
    std::cout << "Expected t_avg err: " << expected_t_avg_err << std::endl;
    
    std::cout << "Temperature L2 err: " << temperature_err << std::endl;    
    std::cout << "Temperature avg err: " << temperature_avg_err << std::endl;
    
    /// check L2 errors
    if(fabs(phi_err - expected_phi_err) > tol)
      throw Dark_Arts_Exception(SUPPORT_OBJECT , "Not calculating phi err identically");
      
    if(fabs(temperature_err - expected_t_err) > tol)
      throw Dark_Arts_Exception(SUPPORT_OBJECT , "Not calculating temperature err identically");
 
    /// check L2 like average errors
    if(fabs(temperature_avg_err - expected_t_avg_err) > tol)
      throw Dark_Arts_Exception(SUPPORT_OBJECT, "Not calculating t avg error identically");
      
    if(fabs(phi_avg_err - expected_phi_avg_err) > tol)
      throw Dark_Arts_Exception(SUPPORT_OBJECT, "Not calculating phi avg error identically");
    
  }
  catch(const Dark_Arts_Exception& da)
  {
    val = -1;
    da.testing_message();
  }
  
  return val;
}
