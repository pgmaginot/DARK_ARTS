#include <iostream>
#include <vector>
#include <math.h>
#include <string>

#include "Input_Reader.h"
#include "Cell_Data.h"
#include "Intensity_Data.h"
#include "Intensity_Moment_Data.h"
#include "Angular_Quadrature.h"
#include "Materials.h"
#include "Fem_Quadrature.h"
#include "Quadrule_New.h"
#include "Dark_Arts_Exception.h"

int main(int argc, char** argv)
{
  int val = 0;
  const double tol = 1.0E-6;
  
  /// expected spacing we should be able to see
  std::vector<int> expected_n_cells(2,0);
  expected_n_cells[0] = 5;
  expected_n_cells[1] = 3;
  const double phi_0 = pow(0.4, 4 );
  const double phi_1 = pow(0.98,4 );
  
  Input_Reader input_reader;    
    
  try{
    input_reader.read_xml(argv[1]);
  }
  catch( const Dark_Arts_Exception& da_exception )
  {
    da_exception.testing_message();
    val = -1;
  }
  
  Cell_Data cell_data( input_reader );   
  Quadrule_New quad_fun;   
  Fem_Quadrature fem_quadrature( input_reader , quad_fun);    
  Angular_Quadrature angular_quadrature( input_reader , quad_fun );  
  Materials materials( input_reader, fem_quadrature , cell_data, angular_quadrature);  
  
  Intensity_Data intensity_ic(cell_data, angular_quadrature, fem_quadrature, materials,input_reader);
  
  Intensity_Moment_Data phi_ic(cell_data,angular_quadrature, fem_quadrature, intensity_ic);
  
  
  
  int cell_cnt = 0;
  const int n_l_mom = angular_quadrature.get_number_of_leg_moments();
  const int n_grp = 1;
  const int n_p = fem_quadrature.get_number_of_interpolation_points();
  
  
  Eigen::VectorXd local_phi_vec(n_p);
  
  std::vector<Eigen::VectorXd> all_flux_moment(n_l_mom , Eigen::VectorXd(n_p) );
  
  try
  {
    for(int reg=0; reg < 2 ; reg++)
    {
      double phi;
      if( reg==0)
        phi = phi_0;
      else
        phi = phi_1;
        
      for(int cell = 0; cell < expected_n_cells[reg] ; cell++)
      {
        phi_ic.get_all_moments(all_flux_moment,cell_cnt , 0 );
        for(int leg = 0; leg < n_l_mom ; leg++)
        {       
          phi_ic.get_cell_angle_integrated_intensity(cell_cnt,n_grp-1,leg,local_phi_vec); 
          for(int el = 0; el < n_p ; el++)
          {
            if(leg ==0)
            {
              if( fabs(local_phi_vec(el) - phi) > tol )
                throw Dark_Arts_Exception( VARIABLE_STORAGE , "Incorrect intensity_moment_data initial condition");
                
              if( fabs(all_flux_moment[leg](el) - phi) > tol )
                throw Dark_Arts_Exception( VARIABLE_STORAGE , "Incorrect intensity_moment_data when retireving all flux moments");
            }
            else
            {
              if(fabs( local_phi_vec(el) ) > tol )
                throw Dark_Arts_Exception( VARIABLE_STORAGE , "Incorrect intensity_moment_data should be isotropic");
                
              if( fabs(all_flux_moment[leg](el) ) > tol )
                throw Dark_Arts_Exception( VARIABLE_STORAGE , "Incorrect zero intensity_moment_data when retireving all flux moments");
            }            
          }          
        }
        cell_cnt++;
      }
    }
    
    /// test the average
    std::vector<double> ref_norm(n_grp,0.);
    phi_ic.get_phi_norm(ref_norm);
    
    double expected_norm = (5.*phi_0 + 3.*phi_1)/8. *1.0E-6;
    std::cout << "Calculated norm for difference calculations: " << ref_norm[0] << " Expected norm for difference calculations: " << expected_norm << std::endl;
    if(fabs( expected_norm - ref_norm[0] )/expected_norm > tol )
    {
      throw Dark_Arts_Exception( VARIABLE_STORAGE , "Not calculating expected phi_norm" );
    }
    
    Intensity_Moment_Data blank_phi(cell_data, angular_quadrature, fem_quadrature, ref_norm);
    for(int cell = 0; cell < 8 ; cell++)
    {
      blank_phi.get_all_moments(all_flux_moment,cell , 0 );
      for(int leg = 0; leg < n_l_mom ; leg++)
      {       
        blank_phi.get_cell_angle_integrated_intensity(cell,n_grp-1,leg,local_phi_vec); 
        for(int el = 0; el < n_p ; el++)
        {
          if(fabs( local_phi_vec(el) ) > tol )
            throw Dark_Arts_Exception( VARIABLE_STORAGE , "Incorrect intensity_moment_data.  Should be isotropic");
                
          if( fabs(all_flux_moment[leg](el) ) > tol )
            throw Dark_Arts_Exception( VARIABLE_STORAGE , "Incorrect zero intensity_moment_data when retireving all flux moments");
        }
      }
    }
    
    Err_Phi err_phi;
    /// difference should be the maximum of phi_0 and phi_1
    phi_ic.normalized_difference(blank_phi, err_phi);
    
    std::cout << "Worst error:  " << err_phi.get_worst_err() << " phi_1: " << phi_1 << std::endl;
    if( fabs( err_phi.get_worst_err() - 1.) > tol)
      throw Dark_Arts_Exception(VARIABLE_STORAGE, "Did not find correct error magnitude");
      
    if( err_phi.get_legendre_moment_with_worst_err() != 0)
      throw Dark_Arts_Exception(VARIABLE_STORAGE, "Largest error found in not scalar flux moment for isotropic intensity");
      
    if(err_phi.get_group_with_worst_err() != 0 )
      throw Dark_Arts_Exception(VARIABLE_STORAGE, "Largest error not in group 0 for a grey problem");
      
    blank_phi = phi_ic;
    
    err_phi.clear();
    phi_ic.normalized_difference(blank_phi, err_phi);
    if( fabs(err_phi.get_worst_err() ) > tol )
      throw Dark_Arts_Exception( VARIABLE_STORAGE, "Should calculate zero error if assignment operator works correctly for Intensity_Moment_Data");
    
    phi_ic.clear_angle_integrated_intensity();
    for(int cell = 0; cell < 8 ; cell++)
    {
      for(int leg = 0; leg < n_l_mom ; leg++)
      {       
        phi_ic.get_cell_angle_integrated_intensity(cell,n_grp-1,leg,local_phi_vec); 
        for(int el = 0; el < n_p ; el++)
        {
          if(fabs( local_phi_vec(el) ) > tol )
            throw Dark_Arts_Exception( VARIABLE_STORAGE , "clear_angle_integrated_intensity did not work correctly");
        }
      }
    }
    
    std::vector<double> phi_moment(n_l_mom,0.);    
    phi_moment[0] = 1.;
    phi_moment[1] = 0.1;
    phi_moment[2] = 2.2;
    phi_moment[3] = 0.5;
    phi_moment[4] = 0.;
    
    for(int l_mom = 0; l_mom < n_l_mom ; l_mom++)
    {
      for(int el = 0; el<n_p ; el++)
        local_phi_vec(el) = phi_moment[l_mom];

      phi_ic.set_cell_angle_integrated_intensity(0,0,l_mom ,local_phi_vec);        
    }
    
    Intensity_Moment_Data i_copy_constructor(phi_ic);
    for(int l_mom = 0; l_mom < n_l_mom ; l_mom++)
    {
      phi_ic.get_cell_angle_integrated_intensity(0,0,l_mom ,local_phi_vec);  
      if(fabs( local_phi_vec(0) - phi_moment[l_mom] ) > tol)
        throw Dark_Arts_Exception(VARIABLE_STORAGE , "Error with set/get or copy constructor");
    }    
  } 
  catch(const Dark_Arts_Exception& da_exception )
  {
    da_exception.testing_message();
    val = -1;
  }
  
  // Return 0 if tests passed, somethnig else if failing
  return val;
}
