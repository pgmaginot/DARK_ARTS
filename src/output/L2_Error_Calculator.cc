/** @file   L2_Error_Calculator.cc
  *   @author pmaginot
  *   @brief A class that numerically approximates the L2 spatial errors of Temperature, Intensity, and Angle Integrated Intensity for MMS problems
 */
         
#include "L2_Error_Calculator.h"
// ##########################################################
// Public functions 
// ##########################################################


L2_Error_Calculator::L2_Error_Calculator(const Angular_Quadrature& angular_quadrature,
    const Fem_Quadrature& fem_quadrature, 
    const Cell_Data& cell_data, 
    const Input_Reader& input_reader)
:
 m_n_dfem( fem_quadrature.get_number_of_interpolation_points() ) , 
 m_n_cells( cell_data.get_total_number_of_cells() ),
 /// m_n_quad_pts(fem_quadrature.get_number_of_integration_points() ),
 m_n_quad_pts(fem_quadrature.get_number_of_source_points() ),
 m_cell_data( cell_data ),
 m_fem_quadrature( fem_quadrature ),
 m_analytic_temperature(input_reader) , 
 m_analytic_intensity(input_reader, angular_quadrature)
{
  /// old way
  // fem_quadrature.get_dfem_integration_points( m_integration_quadrature_pts );
  // fem_quadrature.get_integration_weights( m_integration_quadrature_wts );
  // fem_quadrature.get_dfem_at_integration_points( m_dfem_at_integration_pts );
  
  /// more correct way  
  fem_quadrature.get_source_points(m_integration_quadrature_pts);  
  fem_quadrature.get_source_weights(m_integration_quadrature_wts);
  fem_quadrature.get_dfem_at_source_points(m_dfem_at_integration_pts); 
}


double L2_Error_Calculator::calculate_l2_error(const double time_eval, const Intensity_Moment_Data& phi) const
{
  /// only going to calculate the error in the scalar flux
  double err = 0.;
  
  Eigen::VectorXd num_soln = Eigen::VectorXd::Zero(m_n_dfem);
  
  for(int cell = 0 ; cell < m_n_cells ; cell++)
  {
    phi.get_cell_angle_integrated_intensity(cell,0,0 , num_soln);
    err += calculate_local_l2_error(num_soln , m_cell_data.get_cell_left_edge(cell) , m_cell_data.get_cell_width(cell) , time_eval , PHI );
  }
  
  return sqrt(err);
}

double L2_Error_Calculator::calculate_l2_error(const double time_eval, const Temperature_Data& temperature) const
{
  double err = 0.;
  Eigen::VectorXd num_soln = Eigen::VectorXd::Zero(m_n_dfem);

  for(int cell = 0 ; cell < m_n_cells ; cell++)  {
    temperature.get_cell_temperature(cell,num_soln);
    err += calculate_local_l2_error(num_soln , m_cell_data.get_cell_left_edge(cell) , m_cell_data.get_cell_width(cell) ,  time_eval, TEMPERATURE );
  }
  
  return sqrt(err);
}

double L2_Error_Calculator::calculate_cell_avg_error(const double time_eval, const Intensity_Moment_Data& phi) const
{
  /// only going to calculate the error in the scalar flux
  double err = 0.;
  
  Eigen::VectorXd num_soln = Eigen::VectorXd::Zero(m_n_dfem);
  
  for(int cell = 0 ; cell < m_n_cells ; cell++)
  {
    phi.get_cell_angle_integrated_intensity(cell,0,0 , num_soln);
    err += calculate_local_cell_avg_error(num_soln , m_cell_data.get_cell_left_edge(cell) , m_cell_data.get_cell_width(cell) , time_eval, PHI );
  }  
  
  return sqrt(err);
}

double L2_Error_Calculator::calculate_cell_avg_error(const double time_eval, const Temperature_Data& temperature) const
{
  double err = 0.;
  
  Eigen::VectorXd num_soln = Eigen::VectorXd::Zero(m_n_dfem);

  for(int cell = 0 ; cell < m_n_cells ; cell++)  {
    temperature.get_cell_temperature(cell,num_soln);
    err += calculate_local_cell_avg_error(num_soln , m_cell_data.get_cell_left_edge(cell) , m_cell_data.get_cell_width(cell) , time_eval , TEMPERATURE );
  }
  
  return sqrt(err);
}

double L2_Error_Calculator::calculate_local_l2_error(
  const Eigen::VectorXd& numeric_at_dfem, const double xL , const double dx , const double time , const Data_Type data ) const
{
  double err = 0.;
  
  /// add in dx/2 contribution
  std::vector<double> numeric_soln;
  m_fem_quadrature.evaluate_variable_at_quadrature_pts(numeric_at_dfem , m_dfem_at_integration_pts , numeric_soln);
  double x_eval = 0.;
  double analytic_soln = 0.;
  for(int q = 0 ; q < m_n_quad_pts ; q++)
  {
    x_eval = xL + dx/2.*(1. + m_integration_quadrature_pts[q] );
    switch(data)
    {
      case PHI:
      {
        analytic_soln = m_analytic_intensity.get_mms_phi(x_eval, time);
        break;
      }
      case TEMPERATURE:
      {
        analytic_soln = m_analytic_temperature.get_mms_temperature(x_eval, time);
        break;
      }
      default:
      {
        throw Dark_Arts_Exception(SUPPORT_OBJECT , "Trying to find the L2 error of a weird data type");
        break;
      }
    }
    // err += m_integration_quadrature_wts[q]*( analytic_soln*analytic_soln - 2.*numeric_soln[q]*analytic_soln +  numeric_soln[q]*numeric_soln[q] );
    err += m_integration_quadrature_wts[q]*( analytic_soln - numeric_soln[q])*( analytic_soln - numeric_soln[q]) ;
  }
    err *= dx/2.;
  return err;
}
 
double L2_Error_Calculator::calculate_local_cell_avg_error(
  const Eigen::VectorXd& numeric_at_dfem , const double xL , const double dx , const double time, const Data_Type data ) const
{
  double err = 0.;
  
  std::vector<double> numeric_soln;
  m_fem_quadrature.evaluate_variable_at_quadrature_pts(numeric_at_dfem , m_dfem_at_integration_pts , numeric_soln);
  double x_eval = 0.;
  double analytic_soln = 0.;
  double analytic_average = 0.;
  double numeric_average = 0.;
  for(int q = 0 ; q < m_n_quad_pts ; q++)
  {
    x_eval = xL + dx/2.*(1. + m_integration_quadrature_pts[q] );
    switch(data)
    {
      case PHI:
      {
        analytic_soln = m_analytic_intensity.get_mms_phi(x_eval, time);
        break;
      }
      case TEMPERATURE:
      {
        analytic_soln = m_analytic_temperature.get_mms_temperature(x_eval, time);
        break;
      }
      default:
      {
        throw Dark_Arts_Exception(SUPPORT_OBJECT , "Trying to find the L2 error of a weird data type");
        break;
      }
    }
    analytic_average += m_integration_quadrature_wts[q]* analytic_soln ; 
    numeric_average +=  m_integration_quadrature_wts[q]*numeric_soln[q];
  }
  analytic_average /= 2.;
  numeric_average /= 2.;
  // err = dx*(analytic_average*analytic_average - 2.*analytic_average*numeric_average + numeric_average*numeric_average);
  
  err = dx*(analytic_average-numeric_average)*(analytic_average-numeric_average) ;
  return err;
}

  
