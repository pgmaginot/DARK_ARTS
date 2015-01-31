#include <iostream>
#include <vector>
#include <math.h>
#include <string>

#include "Input_Reader.h"
#include "Angular_Quadrature.h"
#include "Quadrule_New.h"
#include "Dark_Arts_Exception.h"

#include "MMS_Angle_Poly.h"
#include "MMS_Angle_Isotropic.h"
#include "MMS_Space_Poly.h"
#include "MMS_Space_Cos.h"
#include "MMS_Time_Cos.h"
#include "MMS_Time_Poly.h"

int main(int argc, char** argv)
{
  int val = 0;
  const double tol = 1.E-8;
  
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
  Angular_Quadrature angular_quadrature( input_reader , quad_fun );  
  
  std::vector<double> t_trial(3,0.);
  std::vector<double> x_trial(3,0.);
  
  std::vector<double> space_poly_coeff(3,0.);
  std::vector<double> space_cos_coeff(4,0.);
  std::vector<double> time_poly_coeff(2,0.);
  std::vector<double> time_cos_coeff(4,0.);
  std::vector<double> angle_poly_coeff(4,0.);
  
  t_trial[0] = -0.1;
  t_trial[1] = 0.;
  t_trial[2] = 4.5;
  
  x_trial[0] = -0.4;
  x_trial[1] = 0.0;
  x_trial[2] = 2.2;
  
  space_poly_coeff[0] = 1.;
  space_poly_coeff[1] = 0.5;
  space_poly_coeff[2] = 3.;
  
  space_cos_coeff[0] = 1.;
  space_cos_coeff[1] = 2.;
  space_cos_coeff[2] = 3.1;
  space_cos_coeff[3] = 4.;
  
  time_poly_coeff[0] = 3.;
  time_poly_coeff[1] = 5.;
  
  time_cos_coeff[0] = 3.1;
  time_cos_coeff[1] = 1.2;
  time_cos_coeff[2] = 0.1;
  time_cos_coeff[3] = 5.;  
  
  angle_poly_coeff[0] = 0.1;
  angle_poly_coeff[1] = 2.3;
  angle_poly_coeff[2] = 3.;
  angle_poly_coeff[3] = 4.2;

  const double pi = 3.14159265358979323846;
  
  std::shared_ptr<V_MMS_Space> space_ptr;
  std::shared_ptr<V_MMS_Angle> angle_ptr;
  std::shared_ptr<V_MMS_Time> time_ptr;  
   
  double ex_val = 0.; 
  double ex_deriv = 0.;
  
  double calc_val;
  double calc_deriv;
  
  try{
    space_ptr = std::shared_ptr<V_MMS_Space> (new MMS_Space_Poly(space_poly_coeff) );
    for(int i=0; i<3; i++)
    {
      ex_val = 1. + x_trial[i]*0.5 + x_trial[i]*x_trial[i]*3.;
      ex_deriv = 0.5 + 3.*2.*x_trial[i];
      
      if(fabs((ex_val - space_ptr->get_position_component(x_trial[i]) )/ex_val) > tol)
        throw Dark_Arts_Exception(SUPPORT_OBJECT, "Polynomial spatial component wrong");          
        
      if(fabs( (ex_deriv - space_ptr->get_position_derivative(x_trial[i]) )/ex_deriv) > tol)
        throw Dark_Arts_Exception(SUPPORT_OBJECT, "Polynomial spatial deriviative wrong");
    }
    
    space_ptr = std::shared_ptr<V_MMS_Space> (new MMS_Space_Cos(space_cos_coeff) );
    for(int i=0; i<3; i++)
    {
      ex_val = space_cos_coeff[0]*cos( x_trial[i]*pi/space_cos_coeff[1] + space_cos_coeff[2]) + space_cos_coeff[3];        
      ex_deriv = -space_cos_coeff[0]*sin( x_trial[i]*pi/space_cos_coeff[1] + space_cos_coeff[2]) * pi/space_cos_coeff[1];
      
      calc_val = space_ptr->get_position_component(x_trial[i]);
      calc_deriv = space_ptr->get_position_derivative(x_trial[i]);
      
      std::cout << "Calculated: " << calc_val << " Expected:  " << ex_val << std::endl;
      if(fabs( (ex_val - calc_val)/ex_val) > tol)
      {
        throw Dark_Arts_Exception(SUPPORT_OBJECT, "Cos spatial component wrong");          
      }
      if(fabs( (ex_deriv - calc_deriv )/ex_deriv) > tol)
        throw Dark_Arts_Exception(SUPPORT_OBJECT, "Cos spatial deriviative wrong");
    }
    
    time_ptr = std::shared_ptr<V_MMS_Time> (new MMS_Time_Poly(time_poly_coeff) );
    for(int i=0; i<3; i++)
    {
      ex_val = time_poly_coeff[0] + t_trial[i]*time_poly_coeff[1];
      ex_deriv = time_poly_coeff[1];
      
      calc_val = time_ptr->get_time_component(t_trial[i]);
      calc_deriv = time_ptr->get_time_derivative(t_trial[i]);
      
      if(fabs( (ex_val - calc_val )/ex_val) > tol)
        throw Dark_Arts_Exception(SUPPORT_OBJECT, "Polynomial time component wrong");          
        
      if(fabs( (ex_deriv - calc_deriv)/ex_deriv) > tol)
        throw Dark_Arts_Exception(SUPPORT_OBJECT, "Polynomial time deriviative wrong");
    }
    
    
    time_ptr = std::shared_ptr<V_MMS_Time> (new MMS_Time_Cos(time_cos_coeff) );
    for(int i=0; i<3; i++)
    {        
      ex_val = time_cos_coeff[0]*cos( t_trial[i]*pi/time_cos_coeff[1] + time_cos_coeff[2]) + time_cos_coeff[3];        
      ex_deriv = -time_cos_coeff[0]*sin( t_trial[i]*pi/time_cos_coeff[1] + time_cos_coeff[2]) * pi/time_cos_coeff[1];
      
      calc_val = time_ptr->get_time_component(t_trial[i]);
      calc_deriv = time_ptr->get_time_derivative(t_trial[i]);
      
      if(fabs((ex_val - calc_val )/ex_val) > tol)
        throw Dark_Arts_Exception(SUPPORT_OBJECT, "Cos time component wrong");          
        
      if(fabs((ex_deriv - calc_deriv )/ex_deriv) > tol)
        throw Dark_Arts_Exception(SUPPORT_OBJECT, "Cos time derivative wrong");
    }
    
    const int n_dir = angular_quadrature.get_number_of_dir();
    
    angle_ptr = std::shared_ptr<V_MMS_Angle> (new MMS_Angle_Isotropic(2.) );
    double calc_val_2;
    for(int i=0; i<n_dir;i++)
    {
      ex_val = 1./2.;
      
      calc_val = angle_ptr->get_angle_component(i  );
      calc_val_2 = angle_ptr->get_angle_component( angular_quadrature.get_mu(i) );
      if( fabs(calc_val - calc_val_2) > tol)
        throw Dark_Arts_Exception(SUPPORT_OBJECT , "Methods of mu retrieval not equal for MMS components");
        
      if( fabs(calc_val - ex_val ) > tol )
        throw Dark_Arts_Exception(SUPPORT_OBJECT, "Mu component not isotropic for MMS");
    }
    
    angle_ptr = std::shared_ptr<V_MMS_Angle> (new MMS_Angle_Poly(angle_poly_coeff,angular_quadrature) );
    double mu;
    for(int i=0; i<n_dir;i++)
    {
      mu = angular_quadrature.get_mu(i);
      ex_val = angle_poly_coeff[0] + mu*angle_poly_coeff[1] + mu*mu*angle_poly_coeff[2] + mu*mu*mu*angle_poly_coeff[3];
      
      calc_val = angle_ptr->get_angle_component( i );
      calc_val_2 = angle_ptr->get_angle_component( mu );
      if( fabs( (calc_val - calc_val_2)/calc_val) > tol)
        throw Dark_Arts_Exception(SUPPORT_OBJECT , "Methods of mu retrieval not equal for MMS components");
        
      if( fabs( (calc_val - ex_val )/calc_val) > tol )
        throw Dark_Arts_Exception(SUPPORT_OBJECT, "Mu component not isotropic for MMS");
    }
      
  }
  catch(const Dark_Arts_Exception& da)
  {
    val = -1;
    da.testing_message();
  }
  
  return val;
}
