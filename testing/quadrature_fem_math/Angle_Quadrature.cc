#include <iostream>
#include <vector>
#include <math.h>
#include <string>

#include "Input_Reader.h"
#include "Angular_Quadrature.h"
#include "Quadrule_New.h"
#include "Dark_Arts_Exception.h"

int main(int argc, char** argv)
{
  int val = 0;
  
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
  Angular_Quadrature angular_quadrature( input_reader , quad_fun ); 

  
  /// Testing an S8 Gauss-Legendre quadrature  
  const int n_dir = 8;
  const int n_grp = 3;
  const int n_l_mom = 6;
  std::vector<double> expected_grp_low_bnd(n_grp,0.);
  std::vector<double> expected_grp_high_bnd(n_grp,0.);
  expected_grp_low_bnd[0] = 0.; expected_grp_high_bnd[0] = 1.;
  expected_grp_low_bnd[1] = 1.; expected_grp_high_bnd[1] = 2.;
  expected_grp_low_bnd[2] = 2.; expected_grp_high_bnd[2] = 3.;
  
  std::vector<double> expected_mu(n_dir,0.);  
  std::vector<double> expected_w(n_dir, 0.);
  expected_mu[0] =  -0.96028986 ; expected_w[0] =   0.10122854 ; 
  expected_mu[1] =  -0.79666648 ; expected_w[1] =   0.22238103 ; 
  expected_mu[2] =  -0.52553241 ; expected_w[2] =   0.31370665 ; 
  expected_mu[3] =  -0.18343464 ; expected_w[3] =   0.36268378 ; 
  expected_mu[4] =   0.18343464 ; expected_w[4] =   0.36268378 ; 
  expected_mu[5] =   0.52553241 ; expected_w[5] =   0.31370665 ; 
  expected_mu[6] =   0.79666648 ; expected_w[6] =   0.22238103 ; 
  expected_mu[7] =   0.96028986 ; expected_w[7] =   0.10122854 ; 

  /**
    Test all of the following functions:
    int get_number_of_dir(void) const;
    int get_number_of_groups(void) const;
    int get_number_of_leg_moments(void) const;
    double get_leg_poly(const int dir, const int mom) const;
    double get_leg_moment_coeff(const int dir, const int mom) const;
    double get_mu(const int dir) const;
    double get_w(const int dir) const;
    double get_sum_w(void) const;
    
    double get_group_low_bound(const int grp) const;
    double get_group_upper_bound(const int grp) const;
    
    bool has_left_reflection(void) const;  
    
    double most_glance_mu(void) const;
    double most_normal_mu(void) const;  
  */
  try
  {
    if( (angular_quadrature.get_number_of_dir() - n_dir) != 0)
      throw Dark_Arts_Exception(SUPPORT_OBJECT , "Not getting the correct number of angles");
  }
  catch(const Dark_Arts_Exception da_exception)
  {
    da_exception.testing_message() ;
    val = -1;
  }  
  
  try
  {
    if( (angular_quadrature.get_number_of_groups() - n_grp) != 0)
      throw Dark_Arts_Exception(SUPPORT_OBJECT , "Not getting the correct number of groups");
  }
  catch(const Dark_Arts_Exception da_exception)
  {
    da_exception.testing_message() ;
    val = -1;
  }  
  
  try
  {
    if( (angular_quadrature.get_number_of_leg_moments() - n_l_mom) != 0)
      throw Dark_Arts_Exception(SUPPORT_OBJECT , "Not getting the correct number of legendre moments");
  }
  catch(const Dark_Arts_Exception da_exception)
  {
    da_exception.testing_message() ;
    val = -1;
  }  
  
  try
  {
    if( (angular_quadrature.get_sum_w() - 2.0) > 1.0E-6)
      throw Dark_Arts_Exception(SUPPORT_OBJECT , "Unexpected sum of quadrature weights");
  }
  catch(const Dark_Arts_Exception da_exception)
  {
    da_exception.testing_message() ;
    val = -1;
  }  
  
  try
  {
    if( !angular_quadrature.has_left_reflection() )
      throw Dark_Arts_Exception(SUPPORT_OBJECT , "Should have detected left reflective boundary condition");
  }
  catch(const Dark_Arts_Exception da_exception)
  {
    da_exception.testing_message() ;
    val = -1;
  }  
  
  try{
    for(int d = 0; d < n_dir ; d++)
    {
      std::cout << "Expected mu: " << expected_mu[d] << " Got mu= " << angular_quadrature.get_mu(d) << std::endl;
      std::cout << "Expected w: " << expected_w[d] << " Got w= " << angular_quadrature.get_w(d) << std::endl;
      if( fabs( expected_mu[d] - angular_quadrature.get_mu(d) ) > 1.E-6 )
      {
        throw Dark_Arts_Exception( SUPPORT_OBJECT , "mu not matching up");
        break;
      }
      if( fabs( expected_w[d] - angular_quadrature.get_w(d) ) > 1.E-6 )
      {
        throw Dark_Arts_Exception( SUPPORT_OBJECT , "w not matching up");   
        break;        
      }
    }
  }
  catch(const Dark_Arts_Exception da_exception)
  {
    da_exception.testing_message() ;
    val = -1;
  }  
  
  
 
    
  // Return 0 if tests passed, something else if failing
  return val;
}


