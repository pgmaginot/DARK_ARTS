/** @file   Diffusion_Operator.cc
  *   @author pmaginot
  *   @brief Implement a MIP DSA operator
 */

#include "Diffusion_Operator.h"

Diffusion_Operator::Diffusion_Operator(const Input_Reader& input_reader, const Fem_Quadrature& fem_quadrature, 
  const Cell_Data& cell_data, Materials& materials, const Angular_Quadrature& angular_quadrature,
  const int n_groups, const Temperature_Data& t_eval, const bool is_wg_solve)
:
  m_sn_w( angular_quadrature.get_sum_w() ),
  m_np(fem_quadrature.get_number_of_interpolation_points()),
  m_n_cell(cell_data.get_total_number_of_cells() ),
  m_n_mip_loops( is_wg_solve ? ( angular_quadrature.get_number_of_groups() ) : (1) ),
  m_n_mip_sys_size( (is_wg_solve) ? (m_n_cell*m_np) 
      : ( (input_reader.get_lmfga_structure() == NO_COLLAPSE) ? (m_n_cell*m_np*angular_quadrature.get_number_of_groups() ) : (m_n_cell*m_np) ) ),
  m_r_sig_s(Eigen::MatrixXd::Zero(m_np,m_np) ),
  m_r_sig_a(Eigen::MatrixXd::Zero(m_np,m_np) ), 
  
  
  m_no_mu_pos_l_matrix(Eigen::MatrixXd::Zero(m_np,m_np) ),
  m_no_mu_neg_l_matrix(Eigen::MatrixXd::Zero(m_np,m_np) ),
  
  m_sdirk_a_stage(-1.),
  m_dt(-1.),  
  m_time_stage(-1.),
  m_local_assembler(fem_quadrature,input_reader)
{    
  if(is_wg_solve)
  {
    if(n_groups ==1)
    {
      m_diffusion_matrix_creator = std::make_shared<Diffusion_Matrix_Creator_Grey>
        (fem_quadrature,materials,angular_quadrature,t_eval,m_n_cell,input_reader);
        
      m_diffusion_ordering = std::make_shared<Grey_Diffusion_Ordering>(cell_data,angular_quadrature);
    }
    else
    {
      throw Dark_Arts_Exception(TIME_MARCHER , "MF within group scattering MIP not coded");
      
      m_diffusion_ordering = std::make_shared<MG_WG_Diffusion_Ordering>(cell_data,angular_quadrature);
    }  
  }
  else
  {
    LMFGA_STRUCTURE lmfga_type = input_reader.get_lmfga_structure();
    if( lmfga_type == GROUP_COLLAPSE)
    {
      throw Dark_Arts_Exception(TIME_MARCHER , "MF LMFGA with Group Collapse Not Coded");
      
      m_diffusion_ordering = std::make_shared<MG_LMFGA_Group_Collapse_Diffusion_Ordering>(cell_data,angular_quadrature);
    }
    else if( lmfga_type == NO_COLLAPSE)
    {
      throw Dark_Arts_Exception(TIME_MARCHER , "MF LMFGA without Group Collapse Not Coded");
      
      
      m_diffusion_ordering = std::make_shared<MG_LMFGA_No_Collapse_Group_Outer_Diffusion_Ordering>(cell_data,angular_quadrature);
      
      m_diffusion_ordering = std::make_shared<MG_LMFGA_No_Collapse_Cell_Outer_Diffusion_Ordering>(cell_data,angular_quadrature);
    }  
  }
  
  if( !m_diffusion_matrix_creator )
    throw Dark_Arts_Exception(TRANSPORT, "Unable to allocate a Diffusion_Matrix_Creator in Diffusion_Operator");
    
  /// initialize PETSC data
  m_petsc_err = VecCreate(PETSC_COMM_WORLD,&m_mip_solution);   CHKERRV(m_petsc_err);  
  m_petsc_err = VecCreate(PETSC_COMM_WORLD,&m_mip_rhs);   CHKERRV(m_petsc_err); 

  m_petsc_err = VecSetSizes(m_mip_rhs,PETSC_DECIDE,m_n_mip_sys_size); CHKERRV(m_petsc_err);
  m_petsc_err = VecSetSizes(m_mip_solution,PETSC_DECIDE,m_n_mip_sys_size); CHKERRV(m_petsc_err);
  
  m_petsc_err = VecSetFromOptions(m_mip_rhs); CHKERRV(m_petsc_err) ;
  m_petsc_err = VecSetFromOptions(m_mip_solution); CHKERRV(m_petsc_err) ;
  
  /**
    To change values:
    Call this each time we are changing entries:
      VecSetValues(Vec x, int n_vals_to_change, int *array_global_indices, PetscScalar *arr_values, INSERT_VALUES);
    After all enties have been changed:
     VecAssemblyBegin(Vec x);
     VecAssemblyEnd(Vec x);
     
     can define the following with Eigen:
     Eigen::VectorXd eigen_variable_name(n_dim);
     double *dbl_array_vector = &eigen_variable_name(0) ; 
     
     E.g. this will print identical numbers:
     for(int i = 0 ; i < n_dim ; i++)
      std::cout << "Array representation: " << dbl_array_vector[i] << " Eigen representation: " << eigen_varaible_name(i) << std::endl;
     
     --or--
     VecSetValue( Vec v, PetscInt row, PetscScalar value, INSERT_VALUES);
     After all enties have been changed:
     VecAssemblyBegin(Vec v);
     VecAssemblyEnd(Vec v);
     
     Eigen by default stores its matrices in COLUMN-MAJOR-ORDERING
     
     When setting PETSc Matrix values via:
     
    
    std::cout << "Eigen representation: \n" << dum_mat ;
    std::cout << "\n\n Dumb array loop\n" ;
    PetscReal *arr_ptr = &dum_mat(0,0);
    for(int i=0;i<4*4 ; i++)
      std::cout << "arr_ptr[" << i << "]: " << arr_ptr[i] << std::endl;
  */
  
  /// set MIP matrix structure one time for faster assembly
  
  
}

// void Diffusion_Operator::build_rhs(const int mip_loop_number)
// {
  // int cell = 0; 
  // int grp = 0;
  // for(int i=0; i < m_n_mip_blocks ; i++)
  // {
    // m_cell_group_order->get_cell_and_group(i,mip_loop_number,cell,grp);
  // }
  // return;
// } 

// void Diffusion_Operator::build_matrix(const int mip_loop_number)
// {

  // return;
// }

// void Diffusion_Operator::build_matrix_and_rhs(const int mip_loop_number)
// {

// }

Diffusion_Operator::~Diffusion_Operator()
{
  VecDestroy( &m_mip_solution);
  VecDestroy( &m_mip_rhs);
}


bool Diffusion_Operator::check_all_eigen_variables_for_finite(void)
{
  bool bad_eigen_variables = false;
  std::stringstream err;
  err << std::endl;
  
  // if(!m_r_sig_t.allFinite())
  // {
    // bad_eigen_variables = true;
    // err << "m_r_sig_t has non finite values!\n";
  // }
  
  std::cout << err.str();
    
  return bad_eigen_variables;
}
