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
  Materials materials( input_reader, fem_quadrature , cell_data, angular_quadrature );  
  
  try{
    /**
      First we will test that the Materials::get_intensity_source() and Materials::get_temperature_source() functions which evalaute
      the MMS sources at the source points return the correct values
    */
    std::vector<double> source_quad_pts;
    fem_quadrature.get_source_points(source_quad_pts);
    const int n_source_pts = source_quad_pts.size();
    // /*
      // <N_cells> 5 </N_cells>
      // <Left_bound> 1.0 </Left_bound>
      // <Right_bound> 22.0 </Right_bound>
      // <Spacing> Equal   </Spacing>
      
      // <Scattering_opacity_type> Constant_xs 
        // <Constant_value> 1.3 </Constant_value>
      // </Scattering_opacity_type>      
      // <Absorption_opacity_type> Constant_xs 
        // <Constant_value> 0.9 </Constant_value>
      // </Absorption_opacity_type>
      
      // <Cv_type> Constant_cv 
        // <Cv_constant> 2.0 </Cv_constant>
      // </Cv_type>      
      
      // R(\mu) = 1 + mu + mu^2
      // R(x) = 1 + 2.1x + 4.5x^2 + 5x^3
      // T(x) = 1*cos(pi*x/2 + 0) + 1
      // f(t) = 0.5 + 2t
    // */
    
    std::vector<double> ex_temp_space_coeff(4,0.);
    std::vector<double> ex_rad_space_coeff(5,0.);
    std::vector<double> ex_rad_angle_coeff(3,0.);
    std::vector<double> ex_time_coeff(2,0.);
    
    ex_rad_space_coeff[0] = 0.6;
    ex_rad_space_coeff[1] = 2.1;
    ex_rad_space_coeff[2] = 4.5;
    ex_rad_space_coeff[3] = 5.;
    ex_rad_space_coeff[4] = 3.2;
    
    ex_rad_angle_coeff[0] = 1.1;
    ex_rad_angle_coeff[1] = 0.1;
    ex_rad_angle_coeff[2] = 1.3;
    
    ex_temp_space_coeff[0] = 1.;
    ex_temp_space_coeff[1] = 2.;
    ex_temp_space_coeff[2] = 0.;
    ex_temp_space_coeff[3] = 3.;
    
    ex_time_coeff[0] = 0.5;
    ex_time_coeff[1] = 1.2;
    
    std::vector<double> calc_temp_space_coeff;
    std::vector<double> calc_time_coeff;
    std::vector<double> calc_rad_space_coeff;
    std::vector<double> calc_angle_coeff;
    
    input_reader.get_mms_angle_coeff(calc_angle_coeff);
    input_reader.get_mms_time_coeff(calc_time_coeff);
    input_reader.get_mms_radiation_space_coeff(calc_rad_space_coeff);
    input_reader.get_mms_temperature_space_coeff(calc_temp_space_coeff);
    
    /// verify that correct input are being read
    
    for(unsigned int i=0;i < ex_temp_space_coeff.size() ; i++)
      if(fabs( (calc_temp_space_coeff[i] - ex_temp_space_coeff[i] )/ex_temp_space_coeff[i]) > tol )
        throw Dark_Arts_Exception(INPUT, "Difference in MMS temperature space coeff");
        
    for(unsigned int i=0;i < ex_rad_space_coeff.size() ; i++)
      if(fabs( (calc_rad_space_coeff[i] - ex_rad_space_coeff[i] )/ex_rad_space_coeff[i] ) > tol )
        throw Dark_Arts_Exception(INPUT, "Difference in MMS radiation space coeff");
        
    for(unsigned int i=0;i < ex_time_coeff.size() ; i++)
      if(fabs( (calc_time_coeff[i] - ex_time_coeff[i] )/calc_time_coeff[i]) > tol )
        throw Dark_Arts_Exception(INPUT, "Difference in MMS time coeff");
        
    for(unsigned int i=0;i < ex_rad_angle_coeff.size() ; i++)
      if(fabs( (calc_angle_coeff[i] - ex_rad_angle_coeff[i] )/calc_angle_coeff[i]) > tol )
        throw Dark_Arts_Exception(INPUT, "Difference in MMS radiation angle coeff");
  
    const int n_dir = angular_quadrature.get_number_of_dir();
    std::vector<double> dfem_interp_pts;
    fem_quadrature.get_dfem_interpolation_point(dfem_interp_pts);
    const int n_el = dfem_interp_pts.size();
    Eigen::VectorXd temp_dfem;
    temp_dfem = Eigen::VectorXd::Zero(n_el);
    
    const double cv = 2.0;
    const double sig_a = 0.9;
    const double sig_s = 1.3;
    const double c_speed = 1.0;
    const int n_cells = cell_data.get_total_number_of_cells();
    const double sn_w = angular_quadrature.get_sum_w();    
    
    MMS_Intensity calc_intensity( input_reader,angular_quadrature);
    MMS_Temperature calc_temperature( input_reader );
    
    /// verify composite manufactured  solutions for a single point in space/time
    std::shared_ptr<V_MMS_Time> time_component = std::shared_ptr<V_MMS_Time> (new MMS_Time_Poly(ex_time_coeff) );
    std::shared_ptr<V_MMS_Space> rad_space_component = std::shared_ptr<V_MMS_Space> (new MMS_Space_Poly(ex_rad_space_coeff) );
    std::shared_ptr<V_MMS_Space> temp_space_component = std::shared_ptr<V_MMS_Space> (new MMS_Space_Cos(ex_temp_space_coeff) );
    std::shared_ptr<V_MMS_Angle> rad_angle_component = std::shared_ptr<V_MMS_Angle> (new MMS_Angle_Poly(calc_angle_coeff, angular_quadrature) );
   
    const double time_eval = 1.3;
    const double space_eval = 1.4;
    
    const double ex_temp = time_component->get_time_component(time_eval) * temp_space_component->get_position_component(space_eval);
    if( fabs( (ex_temp - calc_temperature.get_mms_temperature(space_eval,time_eval) )/ex_temp) > tol)
      throw Dark_Arts_Exception(SUPPORT_OBJECT , "MMS_tmperature not equal to expected components product");
      
    const double ex_temp_d = time_component->get_time_derivative(time_eval) * temp_space_component->get_position_component(space_eval);
    if( fabs( (ex_temp_d - calc_temperature.get_mms_temperature_time_derivative(space_eval,time_eval) )/ex_temp_d) > tol)
      throw Dark_Arts_Exception(SUPPORT_OBJECT , "MMS_tmperature tiem derivative not equal to expected components product");
      
    for(int d = 0 ; d < n_dir ; d++)
    {
      double ex_intensity = time_component->get_time_component(time_eval) * rad_space_component->get_position_component(space_eval) * rad_angle_component->get_angle_component(d);
       if( fabs((ex_intensity - calc_intensity.get_mms_intensity(space_eval,time_eval,d) )/ex_intensity) > tol)
        throw Dark_Arts_Exception(SUPPORT_OBJECT , "MMS_Intensity not equal to expected components product");
        
      double ex_intensity_dt = time_component->get_time_derivative(time_eval) * rad_space_component->get_position_component(space_eval) * rad_angle_component->get_angle_component(d);
       if( fabs( (ex_intensity_dt - calc_intensity.get_mms_intensity_time_derivative(space_eval,time_eval,d) )/ex_intensity_dt) > tol)
        throw Dark_Arts_Exception(SUPPORT_OBJECT , "MMS_Intensity d wrt time not equal to expected components product");
        
      double ex_intensity_dx = time_component->get_time_component(time_eval) * rad_space_component->get_position_derivative(space_eval) * rad_angle_component->get_angle_component(d);
       if( fabs( (ex_intensity_dx - calc_intensity.get_mms_intensity_space_derivative(space_eval,time_eval,d) )/ex_intensity_dx) > tol)
        throw Dark_Arts_Exception(SUPPORT_OBJECT , "MMS_Intensity d wrt x not equal to expected components product");
    }
    
    double ex_phi = 0.;
    double sum_phi = 0.;
    for(int d = 0; d < n_dir ; d++)
    {
      double ex_intensity = time_component->get_time_component(time_eval) * rad_space_component->get_position_component(space_eval) * rad_angle_component->get_angle_component(d);
      double calc_i = calc_intensity.get_mms_intensity(space_eval,time_eval,d) ;
      if( fabs( (ex_intensity - calc_i )/ex_intensity) > tol)
        throw Dark_Arts_Exception(SUPPORT_OBJECT , "MMS_Intensity not equal to expected components product");
      
      sum_phi += angular_quadrature.get_w(d)* calc_i;
      
      ex_phi += angular_quadrature.get_w(d)* time_component->get_time_component(time_eval) * rad_space_component->get_position_component(space_eval) * rad_angle_component->get_angle_component(d);
    }
    double calc_phi = calc_intensity.get_mms_phi(space_eval,time_eval);
    std::cout << "Expected phi: " << ex_phi << " Calculated_phi: " << calc_phi << std::endl;
 
    if(fabs( (ex_phi - calc_phi )/ex_phi) > tol)
      throw Dark_Arts_Exception(SUPPORT_OBJECT, "MMS_phi does not agree with expected phi");
    
    /// verify temperature source
    double xL , dx, x_dfem;
    for(int c = 0; c < n_cells ; c++)
    {
      xL = cell_data.get_cell_left_edge(c);
      dx = cell_data.get_cell_width(c);
      
      for(int el = 0; el < n_el; el++)
      {        
        x_dfem = xL + dx/2.*(1. + dfem_interp_pts[el] );
        temp_dfem(el) = calc_temperature.get_mms_temperature(x_dfem,time_eval);      
      }
      std::vector<double> calculated_temperature_source(n_source_pts,0.);
      materials.calculate_local_temp_and_position(c,temp_dfem);
      materials.get_temperature_source(time_eval, calculated_temperature_source);
      std::vector<double> expected_t_sources(n_source_pts,0.);
      for(int p=0; p < n_source_pts; p++)
      {
        double x = xL + dx/2.*(1. + source_quad_pts[p]);
        /**
          \f[
            C_v \frac{\partial T}{\partial t} - \sigma_a \left( \phi - \text{m_sn_w} B(T) \right) = S_T
          \f]
        */
        expected_t_sources[p] = cv*calc_temperature.get_mms_temperature_time_derivative(x,time_eval) ;
        expected_t_sources[p] -= sig_a*(calc_intensity.get_mms_phi(x,time_eval) - pow(calc_temperature.get_mms_temperature(x,time_eval),4) );
      
        std::cout << "Calculated t_source: " << calculated_temperature_source[p] << " Expected t_source: " << expected_t_sources[p] << std::endl;
        if(fabs( (expected_t_sources[p] - calculated_temperature_source[p])/expected_t_sources[p])>tol )
          throw Dark_Arts_Exception(SUPPORT_OBJECT, "Difference in MMS temperature evaluations");
      }
    
      for(int d=0; d < n_dir ; d++)
      {
        std::vector<double> ex_i_src(n_source_pts,0.);
        std::vector<double> calc_i_source(n_source_pts,0.);
        materials.get_intensity_source( time_eval, 0 , d, calc_i_source);
        
        for(int p=0; p<n_source_pts; p++)
        {
          double x = xL + dx/2.*(1. + source_quad_pts[p]);
          double mu = angular_quadrature.get_mu(d);
 
          ex_i_src[p] = calc_intensity.get_mms_intensity_time_derivative(x,time_eval,d)/c_speed;
          ex_i_src[p] += mu*calc_intensity.get_mms_intensity_space_derivative(x,time_eval,d);
          ex_i_src[p] += (sig_a+sig_s)*calc_intensity.get_mms_intensity(x,time_eval,d);
          ex_i_src[p] -=  sig_s/sn_w*calc_intensity.get_mms_phi(x,time_eval);
          ex_i_src[p] -= sig_a/sn_w*pow(calc_temperature.get_mms_temperature(x,time_eval),4);
        
          std::cout << "Calculated I source: " << calc_i_source[p] << " Expected i src: " << ex_i_src[p] << std::endl;
          if(fabs( (ex_i_src[p] - calc_i_source[p] )/ex_i_src[p]) > tol )
            throw Dark_Arts_Exception(SUPPORT_OBJECT, "Miscalculating intensity mms sources");
        
        }        
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
