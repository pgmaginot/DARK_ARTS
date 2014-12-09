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
#include "Intensity_Data.h"
#include "Temperature_Data.h"
#include "K_Intensity.h"
#include "K_Temperature.h"
#include "Sweep_Matrix_Creator_Grey.h"

#include "Dark_Arts_Exception.h"

/**
  Goal of this unit test is to check source moment formation, and reaction matrix for SLXS Lobatto scheme
  -This is a MMS problem with a spatially varying sigma_a, constant cv , zero sig_s
*/ 

int main(int argc, char** argv)
{
  int val = 0;
  const double tol = 1.0E-4;
  
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
  const int n_cell = cell_data.get_total_number_of_cells();
  const int n_dir = angular_quadrature.get_number_of_dir();
  
  /// Initialize a Quadrule object to be able to get all of the quadrature we need
  Quadrule_New quad_fun;  
  Fem_Quadrature fem_quadrature( input_reader , quad_fun);  
  Cell_Data cell_data( input_reader );  
  Angular_Quadrature angular_quadrature( input_reader , quad_fun );  
  /// Create a Materials object that contains all opacity, heat capacity, and source objects
  Materials materials( input_reader, fem_quadrature , cell_data, angular_quadrature );  
  Time_Data time_data(input_reader);
  Temperature_Data t_old(fem_quadrature, input_reader, cell_data);  
  Intensity_Data i_old(cell_data, angular_quadrature, fem_quadrature, materials,input_reader);
  Temperature_Data t_star(n_cell,fem_quadrature);
  t_star = t_old;
  
  
  
  try{
  
    K_Temperature kt(n_cell , time_data.get_number_of_stages(), fem_quadrature);
    K_Intensity ki(n_cell ,  time_data.get_number_of_stages(),fem_quadrature, angular_quadrature);
    
    std::shared_ptr<V_Sweep_Matrix_Creator> matrix_creator;
    matrix_creator = std::shared_ptr<V_Sweep_Matrix_Creator> (new Sweep_Matrix_Creator_Grey(
      fem_quadrature, materials, time_data.get_number_of_stages(),
      angular_quadrature.get_sum_w() , angular_quadrature.get_number_of_leg_moments(),
      t_old, i_old, kt, ki, t_star) );
    
    
      

    /**
      functions to test out:
        void update_cell_dependencies(const int cell) override;
        void update_group_dependencies(const int grp) override;           
        void update_direction_dependencies(const int dir) override;
        
        void construct_l_matrix(const double mu, Eigen::MatrixXd& l_matrix);
        void construct_f_vector(const double mu, Eigen::VectorXd& f_vector);
        void set_time_data( const double dt, const double time_stage, const std::vector<double>& rk_a_of_stage_i, const int stage );

        /// retrieve the already constructed \f$ \bar{\bar{\mathbf{R}}}_{\sigma_{t,g}} \f$.  Constructed in \fn update_group_dependencies
        void get_r_sig_t(Eigen::MatrixXd& r_sig_t) const;
        void get_r_sig_s(Eigen::MatrixXd& r_sig_s, const int l_mom) const;

        void get_s_i(Eigen::VectorXd& s_i) const;

        void get_mass_inverse(Eigen::MatrixXd& m_inv) const;
    
        void set_ard_phi_ptr(Intensity_Moment_Data* ard_phi_ptr);

        /// used only by Solution_Saver_K_I
        void calculate_k_i_quantities(void);
        
        /// these assume that the linearization matrices have already been created!!!!
        void k_i_get_r_sig_a(Eigen::MatrixXd& r_sig_a) const;
        void k_i_get_r_sig_s_zero(Eigen::MatrixXd& r_sig_s_zero) const;
        void k_i_get_r_sig_t(Eigen::MatrixXd& r_sig_t) const;
        void k_i_get_s_i(Eigen::VectorXd& s_i) const;
        void k_i_get_planck_vec(Eigen::VectorXd& planck) const;
    */
    const double dt = time_data.get_dt_min();
    const double time_stage = time_data.get_t_start + dt;
    std::vector<double> rk_a;    
    const int stage_num = 0;
    
    matrix_creator->set_time_data(dt,time_stage,rk_a,stage_num);
    
    for(int cell = 0 ; cell < n_cell ; cell++)
    {
      /// update cell dependencies  
      matrix_creator->update_cell_dependencies(cell);
      /// update group (grey) dependencies  
      matrix_creator->update_group_dependencies(0);
      for(int dir = 0; dir < n_dir ; dir++)
      {
        /// update direction dependencies        
        matrix_creator->update_direction_dependencies(dir);
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
