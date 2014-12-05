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
  Materials materials( input_reader, fem_quadrature , cell_data, angular_quadrature );  
  
  try{
    /**
      First we will test that the Materials::get_intensity_source() and Materials::get_temperature_source() functions which evalaute
      the MMS sources at the source points return the correct values
    */
    std::vector<double> source_quad_pts;
    fem_quadrature.get_source_points(source_quad_pts);
    const int n_source_pts = source_quad_pts.size();
    
    std::vector<double> calculated_temperature_source(n_source_pts,0.);    
    std::vector<double> calculated_intensity_source(n_source_pts,0.);
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
    
    const double time = 1.1;
    const double fT = 0.5 + 2.*time;
    const double dt_fT = 2.;
    const double pi = 3.1415926;
    const int n_dir = angular_quadrature.get_number_of_dir();
    std::vector<double> dfem_interp_pts;
    fem_quadrature.get_dfem_interpolation_point(dfem_interp_pts);
    const int n_el = dfem_interp_pts.size();
    Eigen::VectorXd temp_dfem;
    temp_dfem = Eigen::VectorXd::Zero(n_el);
    const double cv = 2.0;
    const double sig_a = 0.9;
    const double sig_s = 1.3;
    const int n_cells = cell_data.get_total_number_of_cells();
    const double sn_w = 2.;
    
    double xL , dx, x_dfem;
    for(int c = 0; c < n_cells ; c++)
    {
      xL = cell_data.get_cell_left_edge(c);
      dx = cell_data.get_cell_width(c);
      
      for(int el = 0; el < n_el; el++)
      {
        x_dfem = xL + dx/2.*(1. + dfem_interp_pts[el] );
        temp_dfem(el) = fT*(1.*cos(pi*x_dfem/2.) + 1.);
      }
      
      materials.calculate_local_temp_and_position(c,temp_dfem);
      materials.get_temperature_source(time, calculated_temperature_source);
      std::vector<double> expected_t_sources(n_source_pts,0.);
      for(int p=0; p < n_source_pts; p++)
      {
        double x = xL + dx/2.*(1. + source_quad_pts[p]);
        
        double Tx = 1.*cos(x*pi/2.) + 1.;
        
        double Rx = 1. + 2.1*x + 4.5*x*x + 5.*pow(x,3);
        double dx_Rx = 2.1 + 2.*4.5*x + 5.*3.*x*x;
        double temp = fT*Tx;
        /**
          \f[
            C_v \frac{\partial T}{\partial t} - \sigma_a \left( \phi - \text{m_sn_w} B(T) \right) = S_T
          \f]
        */
        double phi_angle = 0.;
        double mu , Rmu;
        for(int d=0; d < n_dir ; d++)
        {
          mu = angular_quadrature.get_mu(d);
          Rmu = 1. + mu + mu*mu;
          phi_angle += Rmu*angular_quadrature.get_w(d);
        }
                  
        expected_t_sources[p] = cv*(dt_fT*Tx) - sig_a*(phi_angle*fT*Rx - sn_w*pow(temp,4) );
        
       
        
        // double i_src = 0.;
        // for(int d=0; d < n_dir ; d++)
        // {
          // mu = angular_quadrature.get_mu(d);
          // Rmu = 1. + mu + mu*mu;
          // i_src = dt_fT + mu*(Rmu*dx_Rx*fT) + (sig_a+sig_s)*(Rmu*Rx*fT) - sig_s*(phi_angle*Rx*fT) - sig_a*pow(temp,4);
          // materials.get_intensity_source( time, 0 , d, calculated_intensity_source);
          // // if( fabs(i_src - calculated_intensity_source[p]) > tol )
            // // throw Dark_Arts_Exception(SUPPORT_OBJECT, "Not getting I source moments correct");
        // }

      }
      for(int p=0; p < n_source_pts; p++)
      {
        std::cout << "Expected source: " << expected_t_sources[p] << " Calculated_Source: " << calculated_temperature_source[p] << std::endl;
        if( fabs(expected_t_sources[p] - calculated_temperature_source[p]) > tol )
        {         
          throw Dark_Arts_Exception(SUPPORT_OBJECT, "Not getting T source moments correct");
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
