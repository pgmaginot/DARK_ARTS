/** @file   Diffusion_Operator.cc
  *   @author pmaginot
  *   @brief Implement a MIP DSA operator
 */

#include "Diffusion_Operator.h"


Diffusion_Operator::Diffusion_Operator(const Input_Reader& input_reader, const Fem_Quadrature& fem_quadrature, 
  const Cell_Data& cell_data, Materials& materials, const Angular_Quadrature& angular_quadrature,
  const int n_groups, const Temperature_Data& t_eval, const bool is_wg_solve,
  const double abs_tol, const double rel_tol , const double div_tol)
:
  m_sn_w( angular_quadrature.get_sum_w() ),
  m_np(fem_quadrature.get_number_of_interpolation_points()),
  m_n_cell(cell_data.get_total_number_of_cells() ),
  m_n_mip_blocks(   (is_wg_solve) ? (m_n_cell*angular_quadrature.get_number_of_groups()) : 
  ( (input_reader.get_lmfga_structure() == NO_COLLAPSE) ? (m_n_cell*angular_quadrature.get_number_of_groups() ) : m_n_cell)    ),
  m_n_mip_sys_size( m_n_mip_blocks*m_np) ,
      
  m_r_sig_s(Eigen::MatrixXd::Zero(m_np,m_np) ),
  m_r_sig_a(Eigen::MatrixXd::Zero(m_np,m_np) ), 
  
  m_cell_cm1(Eigen::MatrixXd::Zero(m_np,m_np) ), 
  m_cell_c(Eigen::MatrixXd::Zero(m_np,m_np) ), 
  m_cell_cp1(Eigen::MatrixXd::Zero(m_np,m_np) ), 
  
  m_rhs_local(Eigen::VectorXd::Zero(m_np) ),
  m_phi_new_vec(Eigen::VectorXd::Zero(m_np) ),
  m_phi_old_vec(Eigen::VectorXd::Zero(m_np) ),
  
  m_d_r_cm1(-1.) ,  
  m_d_l_c(-1.) , 
  m_d_r_c(-1.) ,  
  m_d_l_cp1(-1.) , 
  m_kappa_cm12(-1.) ,
  m_kappa_cp12(-1.) , 
  m_dx_cm1(0.) , 
  m_dx_c(0.),
  m_dx_cp1(0.),
  m_cell_data(cell_data),
  m_kappa_calculator( m_np-1, (is_wg_solve ? (input_reader.get_wg_z_mip() ) : (-1.) ) ),
  m_local_assembler(fem_quadrature,input_reader),
  m_pointer_to_eigen_m_rhs( &m_rhs_local(0) ),
  m_pointer_to_m_cm1( &m_cell_cm1(0,0) ),
  m_pointer_to_m_c( &m_cell_c(0,0) ),
  m_pointer_to_m_cp1( &m_cell_cp1(0,0) ),
  m_wg_tol( input_reader.get_within_group_solve_tolerance() )
{    
  m_row_destination = new PetscInt[m_np];
  m_col_destination = new PetscInt[m_np];
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
    }  
  }
  else
  {
    LMFGA_STRUCTURE lmfga_type = input_reader.get_lmfga_structure();
    if( lmfga_type == GROUP_COLLAPSE)
    {
      throw Dark_Arts_Exception(TIME_MARCHER , "MF LMFGA with Group Collapse Not Coded");
      
    }
    else if( lmfga_type == NO_COLLAPSE)
    {
      throw Dark_Arts_Exception(TIME_MARCHER , "MF LMFGA without Group Collapse Not Coded");      
   
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
    
  MatCreate(PETSC_COMM_WORLD, &m_mip_global);
  MatSetSizes(m_mip_global,PETSC_DECIDE,PETSC_DECIDE,m_n_mip_sys_size,m_n_mip_sys_size);
  MatSetType(m_mip_global, MATSEQAIJ);
  /// preallocate maximum required space
  MatSeqAIJSetPreallocation(m_mip_global, 3*m_np , NULL);  
  
  /// This might use direct inversion or may do nothing at all
  // KSPCreate(PETSC_COMM_WORLD,&m_krylov_solver);
  // KSPSetType(m_krylov_solver,KSPPREONLY);
  // KSPGetPC(m_krylov_solver, &m_pre_conditioner);
  // PCSetType(m_pre_conditioner , PCLU);

  /** to set options for Petsc we use the
  PetscOptionsGetxxx(prepend string, option name string , v); commands
  
  */
  KSPCreate(PETSC_COMM_WORLD,&m_krylov_solver);
  // KSPGetPC(m_krylov_solver, &m_pre_conditioner);
  // PCSetType(m_pre_conditioner , PCGAMG);
  // KSPSetFromOptions(m_krylov_solver);
  // KSPSetTolerances(m_krylov_solver, 1.0E-8, 1.0E-16, 1.0, m_n_mip_blocks*m_np);
  
  
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
}

void Diffusion_Operator::kill_petsc_objects()
{
  delete[] m_row_destination;
  delete[] m_col_destination;
  VecDestroy( &m_mip_solution);
  VecDestroy( &m_mip_rhs);
  MatDestroy( &m_mip_global);
  KSPDestroy(&m_krylov_solver);
  return;
}

Diffusion_Operator::~Diffusion_Operator()
{

}

void Diffusion_Operator::dump_matrix()
{
  PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);
  MatView(m_mip_global, PETSC_VIEWER_STDOUT_WORLD);
  return;
}

void Diffusion_Operator::build_rhs(const Intensity_Moment_Data& phi_new, const Intensity_Moment_Data& phi_old)
{
  int cell = 0; 
  int grp = 0;
  for(int i=0; i < m_n_mip_blocks ; i++)
  {
    m_diffusion_ordering->get_cell_and_group(i,cell,grp);
    
    m_diffusion_matrix_creator->set_cell_group_information(cell,grp,m_cell_data.get_cell_width(cell) );
    m_diffusion_matrix_creator->calculate_pseudo_r_sig_a_and_pseudo_r_sig_s(m_r_sig_a,m_r_sig_s);
    
    phi_new.get_cell_angle_integrated_intensity(cell,grp,0,m_phi_new_vec);
    phi_old.get_cell_angle_integrated_intensity(cell,grp,0,m_phi_old_vec);
    
    Eigen::VectorXd temp = Eigen::VectorXd::Zero(m_np);
    temp = m_phi_new_vec - m_phi_old_vec;
    m_rhs_local = m_r_sig_s*temp;
    
    /// save into PETSc vector
    /// first fill row destination index
    for(int el= 0 ; el < m_np ; el++)
    {
      m_row_destination[el] = el + m_np*i;
    }
    m_petsc_err = VecSetValues( m_mip_rhs , m_np , m_row_destination , &m_rhs_local(0), INSERT_VALUES);
  }
  
  m_petsc_err = VecAssemblyBegin(m_mip_rhs);
  m_petsc_err = VecAssemblyEnd(m_mip_rhs);
  
  return;
} 

void Diffusion_Operator::build_matrix()
{
  int cell = 0; 
  int grp = 0;
  for(int i=0; i < m_n_mip_blocks ; i++)
  {
    /// determine cell and group for material property evaluation
    m_diffusion_ordering->get_cell_and_group(i,cell,grp);
    
    /// determine dx of current cell since V_Diffusion_Matrix_Creator::calculate_pseudo_r_sig_a_and_pseudo_r_sig_s requires dx_cell
    if(cell==0)
    {
      m_dx_c=m_cell_data.get_cell_width(0) ;
      m_dx_cp1 = m_cell_data.get_cell_width(1) ; 
    }
    else if(cell == (m_n_cell-1))
    {
      m_dx_cm1 = m_dx_c;
      m_dx_c = m_dx_cp1;
    }
    else
    {
      m_dx_cm1 = m_dx_c;
      m_dx_c = m_dx_cp1;
      m_dx_cp1 = m_cell_data.get_cell_width(cell+1);      
    }
    
    /// send cell number, group number, and dx of current cell to m_diffusion_matrix_creator
    m_diffusion_matrix_creator->set_cell_group_information( cell, grp, m_dx_c);
    
    /// calculate \f$ \mathbf{R}_{\widetilde{\Sigma}_a} \f$ and \mathbf{R}_{\widetilde{\Sigma}_s} \f$
    m_diffusion_matrix_creator->calculate_pseudo_r_sig_a_and_pseudo_r_sig_s(m_r_sig_a,m_r_sig_s);

    /// calculate the various diffusion coefficients we need at cell edges
    m_diffusion_matrix_creator->calculate_d_dependent_quantities(m_d_r_cm1, m_d_l_c , m_d_r_c , m_d_l_cp1, m_s_matrix); 

    int el_i , el_j;
    if(cell==0)
    {
      m_kappa_cm12 = m_kappa_calculator.calculate_boundary_kappa(m_dx_c, m_d_l_c);
      m_kappa_cp12 = m_kappa_calculator.calculate_interior_edge_kappa(m_dx_c, m_dx_cp1 , m_d_r_c , m_d_l_cp1) ;
      
      m_local_assembler.calculate_left_boundary_matrices(m_kappa_cm12, m_kappa_cp12 , 
        m_dx_c, m_dx_cp1, m_d_l_c , m_d_r_c , m_d_l_cp1, m_r_sig_a, m_s_matrix, m_cell_c, m_cell_cp1);  

      /// deterimne row desitanation
      for( el_i = 0 ; el_i < m_np ; el_i++)
        m_row_destination[el_i] = i*m_np + el_i;        
      
      /// cell c first
      for(el_j=0;el_j < m_np ; el_j++)
        m_col_destination[el_j] = i*m_np + el_j;
      
      // respect_sparsity_loading(m_cell_c, m_row_destination ,  m_col_destination);
      dense_dumping(m_cell_c, m_row_destination ,  m_col_destination);
      
      /// cell cp1
      for(el_j=0;el_j < m_np ; el_j++)
        m_col_destination[el_j] = (i+1)*m_np + el_j;
      
      // respect_sparsity_loading(m_cell_cp1, m_row_destination ,  m_col_destination);
      dense_dumping(m_cell_cp1, m_row_destination ,  m_col_destination);
    }
    else if(cell == (m_n_cell-1))
    {
      m_kappa_cm12 = m_kappa_cp12;
      m_kappa_cp12 = m_kappa_calculator.calculate_boundary_kappa(m_dx_c, m_d_r_c) ;
      
      m_local_assembler.calculate_right_boundary_matrices( m_kappa_cm12 , m_kappa_cp12 ,m_dx_cm1 , m_dx_c , 
        m_d_r_cm1,  m_d_l_c , m_d_r_c, m_r_sig_a, m_s_matrix, m_cell_cm1 , m_cell_c); 
        
      /// determine the row we are loading into
      for( el_i = 0 ; el_i < m_np ; el_i++)
        m_row_destination[el_i] = i*m_np + el_i;    

      /// cell c-1 first
      for(el_j=0;el_j < m_np ; el_j++)
        m_col_destination[el_j] = (i-1)*m_np + el_j;
      
      // respect_sparsity_loading(m_cell_cm1, m_row_destination ,  m_col_destination);
      dense_dumping(m_cell_cm1, m_row_destination ,  m_col_destination);
 
      /// cell c
      for(el_j=0;el_j < m_np ; el_j++)
        m_col_destination[el_j] = i*m_np + el_j;     
        
      // respect_sparsity_loading(m_cell_c, m_row_destination ,  m_col_destination);
      dense_dumping(m_cell_c, m_row_destination ,  m_col_destination);
    }
    else
    {
      m_kappa_cm12 = m_kappa_cp12; 
      m_kappa_cp12 = m_kappa_calculator.calculate_interior_edge_kappa(m_dx_c, m_dx_cp1 , m_d_r_c , m_d_l_cp1) ;
      
      m_local_assembler.calculate_interior_matrices( m_kappa_cm12 , m_kappa_cp12 ,m_dx_cm1 , m_dx_c , m_dx_cp1 , 
        m_d_r_cm1,  m_d_l_c , m_d_r_c, m_d_l_cp1 , m_r_sig_a, m_s_matrix, m_cell_cm1 , m_cell_c, m_cell_cp1); 
        
      /// Determine row of the global matrix
      for( el_i = 0 ; el_i < m_np ; el_i++)
        m_row_destination[el_i] = i*m_np + el_i; 
        
        /// cell c-1 first
      for(el_j=0;el_j < m_np ; el_j++)
        m_col_destination[el_j] = (i-1)*m_np + el_j;
        
      // respect_sparsity_loading(m_cell_cm1, m_row_destination ,  m_col_destination);
      dense_dumping(m_cell_cm1, m_row_destination ,  m_col_destination);
      /// cell c
      for(el_j=0;el_j < m_np ; el_j++)
        m_col_destination[el_j] = i*m_np + el_j;
        
      // respect_sparsity_loading(m_cell_c, m_row_destination ,  m_col_destination);
      dense_dumping(m_cell_c, m_row_destination ,  m_col_destination);
 
      /// cell cp1
      for(el_j=0;el_j < m_np ; el_j++)
        m_col_destination[el_j] = (i+1)*m_np + el_j;

      // respect_sparsity_loading(m_cell_cp1, m_row_destination ,  m_col_destination);
      dense_dumping(m_cell_cp1, m_row_destination ,  m_col_destination);
    }
  } 
 
  MatAssemblyBegin(m_mip_global,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(m_mip_global,MAT_FINAL_ASSEMBLY);
    
  return;
}

void Diffusion_Operator::set_time_data(  const double dt, const double time_stage, const double rk_a_ii)
{
  m_diffusion_matrix_creator->set_time_data( dt, time_stage, rk_a_ii);
  /// new time stage, set-up matrix
  build_matrix();
  /// associate the new matrixmak with the krylov solver
  /// no preconditioning
  KSPSetOperators(m_krylov_solver,m_mip_global,m_mip_global,SAME_NONZERO_PATTERN);  
  KSPSetType(m_krylov_solver,KSPCG);
  KSPSetFromOptions(m_krylov_solver);
  KSPSetUp(m_krylov_solver);
  
  return;
}

void Diffusion_Operator::solve_system()
{
  KSPSolve(m_krylov_solver,m_mip_rhs,m_mip_solution);
  
  return;
}

void Diffusion_Operator::update(Intensity_Moment_Data& phi_new , const Intensity_Moment_Data& phi_old)
{
  /// build the driving source (rhs)
  build_rhs(phi_new,phi_old);
  // double l2_rhs =0.;
  // VecNorm(m_mip_rhs , NORM_2 , &l2_rhs);
  // // m_mip_rhs
  // KSPSetTolerances(m_krylov_solver, m_wg_tol , 1.0E-14*l2_rhs , 1.0, 10*m_n_mip_blocks*m_np);
  
  /// invert the diffusion operator
  solve_system();
  /// map the solution back into phi_new
  map_solution_into_intensity_moment_data(phi_new);
  
  return;
}

void Diffusion_Operator::map_solution_into_intensity_moment_data(Intensity_Moment_Data& phi_new)
{
  int cell = 0; 
  int grp = 0;
  
  double* update_ptr;
  
  VecGetArray(m_mip_solution,&update_ptr);
  int n_local = 0;
  
  for(int i=0; i < m_n_mip_blocks ; i++)
  {
    m_diffusion_ordering->get_cell_and_group(i,cell,grp);
    phi_new.add_from_array_pointer( &update_ptr[n_local] , cell, grp);
    
    n_local += m_np;      
  }
  
  VecRestoreArray(m_mip_solution,&update_ptr);
  return;
}

void Diffusion_Operator::make_and_dump_rhs(const Intensity_Moment_Data& phi_new , const Intensity_Moment_Data& phi_old)
{
  build_rhs(phi_new,phi_old);

  PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);
  VecView(m_mip_rhs, PETSC_VIEWER_STDOUT_WORLD);
  return;
}

void Diffusion_Operator::after_rhs_solve_system_and_dump_solution()
{
  solve_system();

  return;
}

void Diffusion_Operator::respect_sparsity_loading(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& local_mat, 
  PetscInt row_dest[], PetscInt col_dest[] )
{
  double mat_sum = 0.;
  for(int i = 0 ; i< m_np ; i++)
  {
    for(int j=0 ; j < m_np ; j++)
    {
      mat_sum += fabs( local_mat(i,j) );
    }
  }
  mat_sum /= double(m_np*m_np);
  
  for(int i = 0 ; i< m_np ; i++)
  {
    for(int j=0 ; j < m_np ; j++)
    {
      if( fabs(local_mat(i,j)) > (1.0E-14* mat_sum) )
        MatSetValues(m_mip_global, 1 , &row_dest[i] , 1, &col_dest[j] , &local_mat(i,j) , INSERT_VALUES ); 
    }
  }
  return;
}

void Diffusion_Operator::dense_dumping(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& local_mat, 
  PetscInt row_dest[], PetscInt col_dest[] )
{
  MatSetValues(m_mip_global, m_np, m_row_destination, m_np, m_col_destination, &local_mat(0,0) , INSERT_VALUES);  
  
  return;
}