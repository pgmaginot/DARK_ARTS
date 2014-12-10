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
#include "MMS_Temperature.h"

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
  
  /// Initialize a Quadrule object to be able to get all of the quadrature we need
  Quadrule_New quad_fun;  
  Fem_Quadrature fem_quadrature( input_reader , quad_fun);  
  Cell_Data cell_data( input_reader );  
  Angular_Quadrature angular_quadrature( input_reader , quad_fun );    
  
  const int n_cell = cell_data.get_total_number_of_cells();
  const int n_dir = angular_quadrature.get_number_of_dir();
  const int n_dfem_p = fem_quadrature.get_number_of_interpolation_points();
  
  /// Create a Materials object that contains all opacity, heat capacity, and source objects
  Materials materials( input_reader, fem_quadrature , cell_data, angular_quadrature );  
  Time_Data time_data(input_reader);
  Temperature_Data t_old(fem_quadrature, input_reader, cell_data);  
  Intensity_Data i_old(cell_data, angular_quadrature, fem_quadrature, materials,input_reader);
  Temperature_Data t_star(n_cell,fem_quadrature);
  t_star = t_old;
  
  std::vector<std::vector<double>> expected_dimensionless_mass( n_dfem_p , std::vector<double> (n_dfem_p , 0.) );   
  std::vector<std::vector<double>> expected_L_mu_positive( n_dfem_p , std::vector<double> (n_dfem_p , 0.) );    
  std::vector<std::vector<double>> expected_L_mu_negative( n_dfem_p , std::vector<double> (n_dfem_p , 0.) ); 
  std::vector<double> expected_f_mu_positive(n_dfem_p,0.);
  std::vector<double> expected_f_mu_negative(n_dfem_p,0.);
    
  
  expected_dimensionless_mass[0][0] = 0.16666667 ; 
  expected_dimensionless_mass[0][1] = 0.00000000 ; 
  expected_dimensionless_mass[0][2] = 0.00000000 ; 
  expected_dimensionless_mass[0][3] = 0.00000000 ; 
  expected_dimensionless_mass[1][0] = 0.00000000 ; 
  expected_dimensionless_mass[1][1] = 0.83333333 ; 
  expected_dimensionless_mass[1][2] = 0.00000000 ; 
  expected_dimensionless_mass[1][3] = 0.00000000 ; 
  expected_dimensionless_mass[2][0] = 0.00000000 ; 
  expected_dimensionless_mass[2][1] = 0.00000000 ; 
  expected_dimensionless_mass[2][2] = 0.83333333 ; 
  expected_dimensionless_mass[2][3] = 0.00000000 ; 
  expected_dimensionless_mass[3][0] = 0.00000000 ; 
  expected_dimensionless_mass[3][1] = 0.00000000 ; 
  expected_dimensionless_mass[3][2] = 0.00000000 ; 
  expected_dimensionless_mass[3][3] = 0.16666667 ;   
  
  expected_L_mu_positive[0][0] = 0.50000000 ; 
  expected_L_mu_positive[0][1] = 0.67418083 ; 
  expected_L_mu_positive[0][2] = -0.25751416 ; 
  expected_L_mu_positive[0][3] = 0.08333333 ; 
  expected_L_mu_positive[1][0] = -0.67418083 ; 
  expected_L_mu_positive[1][1] = 0.00000000 ; 
  expected_L_mu_positive[1][2] = 0.93169499 ; 
  expected_L_mu_positive[1][3] = -0.25751416 ; 
  expected_L_mu_positive[2][0] = 0.25751416 ; 
  expected_L_mu_positive[2][1] = -0.93169499 ; 
  expected_L_mu_positive[2][2] = -0.00000000 ; 
  expected_L_mu_positive[2][3] = 0.67418083 ; 
  expected_L_mu_positive[3][0] = -0.08333333 ; 
  expected_L_mu_positive[3][1] = 0.25751416 ; 
  expected_L_mu_positive[3][2] = -0.67418083 ; 
  expected_L_mu_positive[3][3] = 0.50000000 ; 
  expected_L_mu_negative[0][0] = -0.50000000 ; 
  expected_L_mu_negative[0][1] = 0.67418083 ; 
  expected_L_mu_negative[0][2] = -0.25751416 ; 
  expected_L_mu_negative[0][3] = 0.08333333 ; 
  expected_L_mu_negative[1][0] = -0.67418083 ; 
  expected_L_mu_negative[1][1] = 0.00000000 ; 
  expected_L_mu_negative[1][2] = 0.93169499 ; 
  expected_L_mu_negative[1][3] = -0.25751416 ; 
  expected_L_mu_negative[2][0] = 0.25751416 ; 
  expected_L_mu_negative[2][1] = -0.93169499 ; 
  expected_L_mu_negative[2][2] = -0.00000000 ; 
  expected_L_mu_negative[2][3] = 0.67418083 ; 
  expected_L_mu_negative[3][0] = -0.08333333 ; 
  expected_L_mu_negative[3][1] = 0.25751416 ; 
  expected_L_mu_negative[3][2] = -0.67418083 ; 
  expected_L_mu_negative[3][3] = -0.50000000 ; 
  expected_f_mu_positive[0] = 1.00000000 ; 
  expected_f_mu_positive[1] = 0.00000000 ; 
  expected_f_mu_positive[2] = -0.00000000 ; 
  expected_f_mu_positive[3] = 0.00000000 ; 
  expected_f_mu_negative[0] = 0.00000000 ; 
  expected_f_mu_negative[1] = -0.00000000 ; 
  expected_f_mu_negative[2] = 0.00000000 ; 
  expected_f_mu_negative[3] = -1.00000000 ; 
  
  std::vector<double> source_quad_pt;
  std::vector<double> source_quad_wt;
  std::vector<double> dfem_at_source;
  const int n_source_q = fem_quadrature.get_number_of_source_points();
  fem_quadrature.get_source_points(source_quad_pt);
  fem_quadrature.get_source_weights(source_quad_wt);
  fem_quadrature.get_dfem_at_source_points(dfem_at_source);
  
  try{
  
    K_Temperature kt(n_cell , time_data.get_number_of_stages(), fem_quadrature);
    K_Intensity ki(n_cell ,  time_data.get_number_of_stages(),fem_quadrature, angular_quadrature);
    
    std::shared_ptr<V_Sweep_Matrix_Creator> matrix_creator;
    matrix_creator = std::shared_ptr<V_Sweep_Matrix_Creator> (new Sweep_Matrix_Creator_Grey(
      fem_quadrature, materials, time_data.get_number_of_stages(),
      angular_quadrature.get_sum_w() , angular_quadrature.get_number_of_leg_moments(),
      t_old, i_old, kt, ki, t_star) );
       
    MMS_Temperature temp_mms(input_reader);
    MMS_Intensity i_mms(input_reader,angular_quadrature);
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
    const double c_speed = materials.get_c();
    const double dt = time_data.get_dt_min();
    const double time_stage = time_data.get_t_start() + dt;
    std::vector<double> rk_a(1,0.);  

    rk_a[0] = 1.;
    
    const int stage_num = 0;
    
    matrix_creator->set_time_data(dt,time_stage,rk_a,stage_num);
    
    std::vector<double> dfem_pt;
    fem_quadrature.get_dfem_interpolation_point(dfem_pt);
    
    for(int cell = 0 ; cell < n_cell ; cell++)
    {
      /// change t_star to be the temperature at t = t_start + dt_min = time_stage;
      double dx = cell_data.get_cell_width(cell);
      double xL = cell_data.get_cell_left_edge(cell);
      Eigen::VectorXd t_star_vec = Eigen::VectorXd::Zero(n_dfem_p);
      for(int el = 0;el < n_dfem_p ; el++)
      {
        double x_dfem = xL + dx/2.*(1. + dfem_pt[el]); 
        t_star_vec(el) = temp_mms.get_mms_temperature(x_dfem,time_stage);
      }
      t_star.set_cell_temperature(cell,t_star_vec);
      
      /// update cell dependencies  
      matrix_creator->update_cell_dependencies(cell);
      /// update group (grey) dependencies  
      matrix_creator->update_group_dependencies(0);
      
      Eigen::MatrixXd Minv = Eigen::MatrixXd(n_dfem_p,n_dfem_p);
      matrix_creator->get_mass_inverse(Minv);
      
      Eigen::MatrixXd r_sig_t = Eigen::MatrixXd(n_dfem_p,n_dfem_p);
      matrix_creator->get_r_sig_t(r_sig_t);
      /// check mass matrix inverse
      for(int i = 0; i< n_dfem_p ; i++)
      {
        for(int j=0 ; j<n_dfem_p ; j++)
        {          
          if(j==i)
          {
            if(fabs(1./(expected_dimensionless_mass[i][j]*dx/2. ) - Minv(i,j) ) > tol )              
              throw Dark_Arts_Exception(FEM , "Miscalculating mass inverse");
          }
          else
          {
            if(fabs(expected_dimensionless_mass[i][j] - Minv(i,j) ) > tol )
              throw Dark_Arts_Exception(FEM , "Miscalculating mass inverse");
          }
        }
      }
      /// check \mathbf{R}_{\sigma_{\tau}} 
      /** 
        \f[ \mathbf{R}_{\sigma_{\tau}} = \frac{1}{c\Delta t}\mathbf{M} + \mathbf{R}_{\sigma_t} \f]
      */
      for(int i = 0; i< n_dfem_p ; i++)
      {
        for(int j=0 ; j<n_dfem_p ; j++)
        {
          if(j==i)
          {             
            double sig_t = 1./pow(t_star_vec(j),3);
            double ex_val = expected_dimensionless_mass[i][j]*dx/(2.*c_speed*dt) + 
              expected_dimensionless_mass[i][j]*dx/2.*sig_t;
            if(fabs( ex_val - r_sig_t(i,j) ) > tol )          
            {              
              throw Dark_Arts_Exception(FEM , "Miscalculating r_sig_t");
            }
          }
          else
          {
            if(fabs( r_sig_t(i,j) ) > tol )
              throw Dark_Arts_Exception(FEM , "Off diagonals of r_sig_t are wrong");
          }
        }
      }
      
      for(int dir = 0; dir < n_dir ; dir++)
      {
        /// update direction dependencies  
        Eigen::MatrixXd L = Eigen::MatrixXd::Zero(n_dfem_p,n_dfem_p);
        Eigen::VectorXd f = Eigen::VectorXd::Zero(n_dfem_p);
        
        double mu = angular_quadrature.get_mu(dir);
        matrix_creator->update_direction_dependencies(dir);
        matrix_creator->construct_l_matrix(mu,L);
        matrix_creator->construct_f_vector(mu,f);
        
        if(mu < 0.)
        {
          for(int i=0; i < n_dfem_p; i++)
          {
            for(int j=0; j<n_dfem_p ; j++)
            {
              if(fabs(L(i,j) - mu*expected_L_mu_negative[i][j]) > tol )
                throw Dark_Arts_Exception(FEM, "L for mu < 0 is wrong");
            }
            if(fabs(f(i) - mu*expected_f_mu_negative[i]) > tol )
              throw Dark_Arts_Exception(FEM, "f for mu < 0 is wrong");
          }
        }
        else
        {
          for(int i=0; i < n_dfem_p ; i++)
          {
            for(int j=0; j<n_dfem_p ; j++)
            {
              if(fabs(L(i,j) - mu*expected_L_mu_positive[i][j]) > tol )
                throw Dark_Arts_Exception(FEM, "L for mu > 0 is wrong");
            }
            if(fabs(f(i) - mu*expected_f_mu_positive[i]) > tol )
              throw Dark_Arts_Exception(FEM, "f for mu > 0 is wrong");
          }
        }
        
        Eigen::VectorXd s_i = Eigen::VectorXd::Zero(n_dfem_p);        
        matrix_creator->k_i_get_s_i(s_i);
        Eigen::VectorXd my_s_i = Eigen::VectorXd::Zero(n_dfem_p);
        int pnt = 0;
        for(int el = 0; el < n_dfem_p ; el++)
        {
          /// evaluate S_I and form to verify matrix_creator return value;
          
          // source_quad_pt; source_quad_wt;
          for(int q = 0; q < n_source_q ; q++)
          {
            double x_q = xL + dx/2.*(1.+source_quad_pt[q]);
            double t = temp_mms.get_mms_temperature(x_q,time_stage);
            double sig_a = 1. / pow(t,3);
            double i_d = i_mms.get_mms_intensity(x_q,time_stage,dir);
            double i_dx = i_mms.get_mms_intensity_space_derivative(x_q , time_stage,dir);
            double i_dt = i_mms.get_mms_intensity_time_derivative(x_q,time_stage,dir);
            
            /// take adavantage of this being a pure absorber problem
            double source_q = i_dt/c_speed + mu*i_dx + sig_a*i_d - sig_a*pow(t,4)/2.;
            my_s_i(el) += source_quad_wt[q]*source_q*dfem_at_source[pnt];
            pnt++;
          }
          my_s_i(el) *= dx/2.;
          
          if( fabs(my_s_i(el) - s_i(el) ) > tol)
            throw Dark_Arts_Exception(FEM, "Not calculating driving intensity source moments correctly");
        }
        
      }
      /** test these simple features that we would need while finding k_I
        void calculate_k_i_quantities(void);
        
        /// these assume that the linearization matrices have already been created!!!!
        void k_i_get_r_sig_a(Eigen::MatrixXd& r_sig_a) const;
        void k_i_get_r_sig_s_zero(Eigen::MatrixXd& r_sig_s_zero) const;
        void k_i_get_r_sig_t(Eigen::MatrixXd& r_sig_t) const;
        void k_i_get_s_i(Eigen::VectorXd& s_i) const;
        void k_i_get_planck_vec(Eigen::VectorXd& planck) const;
      */
      matrix_creator->calculate_k_i_quantities();
      Eigen::MatrixXd k_i_mat = Eigen::MatrixXd::Zero(n_dfem_p , n_dfem_p);
      matrix_creator->k_i_get_r_sig_a(k_i_mat);
      for(int i = 0; i< n_dfem_p ; i++)
      {
        for(int j=0 ; j<n_dfem_p ; j++)
        {
          if(j==i)
          {             
            double sig_t = 1./pow(t_star_vec(j),3);
            double ex_val = expected_dimensionless_mass[i][j]*dx/2.*sig_t;
            if(fabs( ex_val - k_i_mat(i,j) ) > tol )          
            {              
              throw Dark_Arts_Exception(FEM , "Miscalculating r_sig_a");
            }
          }
          else
          {
            if(fabs( k_i_mat(i,j) ) > tol )
              throw Dark_Arts_Exception(FEM , "Off diagonals of r_sig_a are nonzero");
          }
        }
      }
      
      matrix_creator->k_i_get_r_sig_s_zero(k_i_mat);
      for(int i = 0; i< n_dfem_p ; i++)
      {
        for(int j=0 ; j<n_dfem_p ; j++)
        {
          if(fabs( k_i_mat(i,j) ) > tol )
              throw Dark_Arts_Exception(FEM , "r_sig_s is nonzero");
        }
      }
      
      matrix_creator->k_i_get_r_sig_t(k_i_mat);
      for(int i = 0; i< n_dfem_p ; i++)
      {
        for(int j=0 ; j<n_dfem_p ; j++)
        {
          if(j==i)
          {             
            double sig_t = 1./pow(t_star_vec(j),3);
            double ex_val = expected_dimensionless_mass[i][j]*dx/2.*sig_t;
            if(fabs( ex_val - k_i_mat(i,j) ) > tol )          
            {              
              throw Dark_Arts_Exception(FEM , "Miscalculating r_sig_t");
            }
          }
          else
          {
            if(fabs( k_i_mat(i,j) ) > tol )
              throw Dark_Arts_Exception(FEM , "Off diagonals of r_sig_t are nonzero");
          }
        }
      }
      
      Eigen::VectorXd k_i_vec = Eigen::VectorXd::Zero(n_dfem_p);
      matrix_creator->k_i_get_planck_vec(k_i_vec);
      for(int el = 0; el < n_dfem_p ; el++)
      {
        if(fabs(k_i_vec(el) - pow(t_star_vec(el),4)/2.) > tol)
          throw Dark_Arts_Exception(FEM, "Incorrect Planck vector for planck");
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
