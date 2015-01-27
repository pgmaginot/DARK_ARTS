#include "Transport_Sweep.h"

Transport_Sweep::Transport_Sweep(const Fem_Quadrature& fem_quadrature, 
  const Cell_Data& cell_data, 
  Materials& materials, 
  const Angular_Quadrature& angular_quadrature, 
  const int n_stages,
  const Temperature_Data& t_old, 
  const Intensity_Data& i_old,
  const K_Temperature& kt, 
  K_Intensity& ki,
  const Temperature_Data& t_star,
  const Input_Reader& input_reader)
  :
  m_n_cells( cell_data.get_total_number_of_cells() ),
  m_n_groups( angular_quadrature.get_number_of_groups()  ),
  m_n_dir( angular_quadrature.get_number_of_dir() ),
  m_n_l_mom( angular_quadrature.get_number_of_leg_moments() ),
  m_np( fem_quadrature.get_number_of_interpolation_points() ), 
  m_ang_quad(angular_quadrature),
  m_matrix_scratch( Eigen::MatrixXd(m_np,m_np) ),
  m_vector_scratch( Eigen::VectorXd::Zero(m_np) ),
  m_local_phi(m_n_l_mom, Eigen::VectorXd::Zero(m_np) ),
  m_rhs_vec( Eigen::VectorXd::Zero(m_np) ),
  m_lhs_mat( Eigen::MatrixXd(m_np,m_np)),
  m_local_soln( Eigen::VectorXd::Zero(m_np)  ),
  m_time(-1.),
  m_psi_in(m_n_groups,m_n_dir),
  m_left_reflecting( angular_quadrature.has_left_reflection() ),
  m_k_i_sweep(false),
  m_krylov_sweep(false),
  m_sweep_type_set(false)
{
  if(m_n_groups > 1)
  {
    m_sweep_matrix_creator = std::shared_ptr<V_Sweep_Matrix_Creator> 
    (new Sweep_Matrix_Creator_MF(fem_quadrature, materials, n_stages, angular_quadrature.get_sum_w() , m_n_l_mom, m_n_groups,
      t_old, i_old, kt, ki,t_star) );
  }
  else
  {
    m_sweep_matrix_creator = std::shared_ptr<V_Sweep_Matrix_Creator> (new Sweep_Matrix_Creator_Grey(fem_quadrature, materials, n_stages, angular_quadrature.get_sum_w() ,
      m_n_l_mom, t_old, i_old, kt, ki,t_star) );
  }
  /// has a ptr to the object that actually owns s_i
  m_fixed_source = std::shared_ptr<V_Sweep_Fixed_Source> (new Sweep_Fixed_Source_Linearization(fem_quadrature,m_sweep_matrix_creator) );
  m_no_source = std::shared_ptr<V_Sweep_Fixed_Source> (new Sweep_Fixed_Source_None(fem_quadrature) );
  
  /// this guy needs a reference to k_i
  m_k_i_saver = std::shared_ptr<V_Solution_Saver> (new Solution_Saver_K_I(fem_quadrature,m_sweep_matrix_creator,angular_quadrature,ki, materials.get_c() ) );
  m_angle_integrated_saver = std::shared_ptr<V_Solution_Saver> (new Solution_Saver_Flux_Moments(fem_quadrature,angular_quadrature) );
  
  /// initialize boundary condition objects
  switch( input_reader.get_radiation_bc_type_left() )
  {
    case VACUUM_BC:
    {
      m_sweep_bc_left = std::shared_ptr<V_Transport_BC> ( new Transport_BC_Vacuum() );
      break;
    }
    case REFLECTIVE_BC:
    {
      m_sweep_bc_left = std::shared_ptr<V_Transport_BC> ( new Transport_BC_Reflective() );
      break;
    }
    case INCIDENT_BC:
    {
      if(input_reader.get_left_bc_value_type() == INCIDENT_CURRENT)
      {
        if(m_n_groups > 1)
        {
          m_sweep_bc_left = std::shared_ptr<V_Transport_BC> ( new Transport_BC_MF_Current(
              angular_quadrature,
              input_reader.get_left_bc_angle_dependence() ,
              input_reader.get_left_bc_start_time() , 
              input_reader.get_left_bc_end_time() ,
              input_reader.get_left_bc_constant() ,
              input_reader.get_left_bc_energy_dependence()
                ) );
        }
        else
        {
          m_sweep_bc_left = std::shared_ptr<V_Transport_BC> ( new Transport_BC_Grey_Current(
              angular_quadrature,
              input_reader.get_left_bc_angle_dependence() ,
              input_reader.get_left_bc_start_time() , 
              input_reader.get_left_bc_end_time() ,
              input_reader.get_left_bc_constant() 
                ) );
        }
      }
      else if(input_reader.get_left_bc_value_type() == INCIDENT_TEMPERATURE)
      {
        if(m_n_groups>1)
        {      
          m_sweep_bc_left = std::shared_ptr<V_Transport_BC> ( new Transport_BC_MF_Planckian(
            materials,
            angular_quadrature,
            input_reader.get_left_bc_angle_dependence() ,
            input_reader.get_left_bc_start_time() , 
            input_reader.get_left_bc_end_time() ,
            input_reader.get_left_bc_constant() ,
            input_reader.get_left_bc_energy_dependence()
              ) );
        }
        else
        {
          m_sweep_bc_left = std::shared_ptr<V_Transport_BC> ( new Transport_BC_Grey_Planckian(
            materials,
            angular_quadrature,
            input_reader.get_left_bc_angle_dependence() ,
            input_reader.get_left_bc_start_time() , 
            input_reader.get_left_bc_end_time()  ,
            input_reader.get_left_bc_constant()         
              ) );
        }
      }
      else
        throw Dark_Arts_Exception( SUPPORT_OBJECT ,  "Unknown left transport BC type in Transport_Sweep constructor" );
        
      break;
    }
    case MMS_BC:
    {
      m_sweep_bc_left = std::shared_ptr<V_Transport_BC> (new Transport_BC_MMS( angular_quadrature, input_reader , cell_data.get_cell_left_edge(0) ) );
      break;
    }
    case INVALID_RADIATION_BC_TYPE:
    {
      throw Dark_Arts_Exception( SUPPORT_OBJECT ,  "Made it into Transport_Sweep constructor with invalid left bc");
      break;
    }
  }
  
  switch( input_reader.get_radiation_bc_type_right() )
  {
    case VACUUM_BC:
    {
      m_sweep_bc_right = std::shared_ptr<V_Transport_BC> (new Transport_BC_Vacuum() );
      break;
    }
    case REFLECTIVE_BC:
    {
      throw Dark_Arts_Exception( SUPPORT_OBJECT ,  "Reflective boundary condition on right face of slab caught in Transport_Sweep constructor");
      break;
    }
    case INCIDENT_BC:
    {
      if(input_reader.get_right_bc_value_type() == INCIDENT_CURRENT)
      {
        if(m_n_groups > 1)
        {
          m_sweep_bc_right = std::shared_ptr<V_Transport_BC> ( new Transport_BC_MF_Current(
              angular_quadrature,
              input_reader.get_right_bc_angle_dependence() ,
              input_reader.get_right_bc_start_time() , 
              input_reader.get_right_bc_end_time() ,
              input_reader.get_right_bc_constant() ,
              input_reader.get_right_bc_energy_dependence()
                ) );
        }
        else
        {
          m_sweep_bc_right = std::shared_ptr<V_Transport_BC> ( new Transport_BC_Grey_Current(
              angular_quadrature,
              input_reader.get_right_bc_angle_dependence() ,
              input_reader.get_right_bc_start_time() , 
              input_reader.get_right_bc_end_time() ,
              input_reader.get_right_bc_constant() 
                ) );
        }
      }
      else if(input_reader.get_right_bc_value_type() == INCIDENT_TEMPERATURE)
      {
        if(m_n_groups>1)
        {      
          m_sweep_bc_right = std::shared_ptr<V_Transport_BC> (new Transport_BC_MF_Planckian(
            materials,
            angular_quadrature,
            input_reader.get_right_bc_angle_dependence() ,
            input_reader.get_right_bc_start_time() , 
            input_reader.get_right_bc_end_time(),
            input_reader.get_right_bc_constant() ,
            input_reader.get_right_bc_energy_dependence()
              ) );
        }
        else
        {
          m_sweep_bc_right = std::shared_ptr<V_Transport_BC> (new Transport_BC_Grey_Planckian(
            materials,
            angular_quadrature,
            input_reader.get_right_bc_angle_dependence() ,
            input_reader.get_right_bc_start_time() , 
            input_reader.get_right_bc_end_time()  ,      
            input_reader.get_right_bc_constant()  ) );
        }
      }
      else
      {
        throw Dark_Arts_Exception( SUPPORT_OBJECT , "Invalid Incident_Value type for right boundary in Transport_Sweep constructor");
      }
      break;
    }
    case MMS_BC:
    {
      m_sweep_bc_right = std::shared_ptr<V_Transport_BC> (new Transport_BC_MMS( angular_quadrature, input_reader , 
        cell_data.get_cell_left_edge(m_n_cells-1) + cell_data.get_cell_width(m_n_cells-1) ) );
      break;
    }
    case INVALID_RADIATION_BC_TYPE:
    {
      throw Dark_Arts_Exception( SUPPORT_OBJECT ,"Made it into Transport_Sweep constructor with invalid right bc");
      break;
    }
  }
  
  /// to set the time, 
  m_sweep_saver = m_angle_integrated_saver;
}

void Transport_Sweep::set_ard_phi_ptr(Intensity_Moment_Data* ard_phi_ptr)
{
  if(!ard_phi_ptr)
    throw Dark_Arts_Exception( SUPPORT_OBJECT , "Attempting to set ard_phi_tr to a NULL ptr");
  
  m_sweep_matrix_creator->set_ard_phi_ptr(ard_phi_ptr);
  return;
}

void Transport_Sweep::set_sweep_type(const bool is_krylov, const bool is_get_k_i)
{
  if(is_krylov)
  {
    m_sweep_source = m_no_source;
    m_krylov_sweep = true;
  }
  else
  {
    m_sweep_source = m_fixed_source;
    m_krylov_sweep = false;
  }
  
  if(is_get_k_i)
  {
    if(is_krylov)
    {  
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Cannot calculate k_i during a Krylov mode sweep");
    }
    
    m_k_i_sweep = true;
    m_sweep_saver = m_k_i_saver;
  }
  else
  {
    m_k_i_sweep = false;
    m_sweep_saver = m_angle_integrated_saver; 
  }
  
  m_sweep_type_set = true;
  
  return;
}

void Transport_Sweep::sweep_mesh(const Intensity_Moment_Data& phi_old, Intensity_Moment_Data& phi_new)
{
  /** perform a single transport sweep across the mesh
    * Do not save the full intensity vector
    * Full intensity vector is only needed when calculating \f$ k_I \f$, and is only needed locally to that solve
    * this means there will be a separate transport sweep to calculate \f$ k_I \f$
    
    Loop over:
      cells
        groups
          directions
  */
  
  if(!m_sweep_type_set)
    throw Dark_Arts_Exception(TRANSPORT, "Must set transport sweep type prior to performing a sweep");
    
  if(!m_k_i_sweep)
    phi_new.clear_angle_integrated_intensity();
    
  /// starting and ending number of cells
  int cell_end, cell_start;
  /// increment or decrement depending on direction set
  int incr = 0;
  
  int dir, d_offset;
  double mu;
  
  /// dir_set = 0, loop over negative mu
  /// dir_set = 1, loop over positive mu    
  for(int dir_set = 0; dir_set < 2 ; dir_set++)
  {
    if(dir_set == 0)
    {
      cell_start = m_n_cells-1;
      cell_end = 0;
      incr = -1;
      d_offset = 0;      
      get_neg_mu_boundary_conditions(m_krylov_sweep);
    }
    else
    {
      cell_start = 0;
      cell_end = m_n_cells - 1;
      incr = +1;
      d_offset = m_n_dir/2;
      get_pos_mu_boundary_conditions(m_krylov_sweep);
    }
        
    for( int cell = cell_start ; cell != (cell_end+incr) ; cell += incr)
    {
      m_sweep_matrix_creator->update_cell_dependencies(cell); 
        
      for(int grp = 0; grp < m_n_groups ; grp++)
      {
        /// avoid calculating cross sections as much as possible, call m_sweep_matrix_creator->update_group_dependencies
        m_sweep_matrix_creator->update_group_dependencies(grp);
          
        /// populate V_Sweep_Matrix_Creator m_k_i_r_sig_t and m_k_i_r_sig_s_zero
        if(m_k_i_sweep)
          m_sweep_matrix_creator->calculate_k_i_quantities();  
        
        /// get all the moments of this group's angular flux moments
        phi_old.get_all_moments(m_local_phi,cell,grp);
        
        for(int d = 0; d<(m_n_dir/2) ; d++)
        {
          dir = d_offset + d;
          mu = m_ang_quad.get_mu(dir);          
          m_sweep_matrix_creator->update_direction_dependencies(dir);
        
          /// build LHS matrix
          /// overwrite LHS with \f$ \mathbf[L} \f$
          m_sweep_matrix_creator->construct_l_matrix(mu,m_lhs_mat);
          
          /// add in \f$ \mathbf{R}_{\sigma_t} \f$
          m_sweep_matrix_creator->get_r_sig_t(m_matrix_scratch);   
          
          m_lhs_mat += m_matrix_scratch; 
          
          /// build RHS vector
          /// overwrite with incident flux contribution
          m_sweep_matrix_creator->construct_f_vector(mu,m_rhs_vec);
          
          m_rhs_vec *= m_psi_in(grp,dir);     
          
          /// add in source moments
          m_sweep_source->get_source(m_vector_scratch);          
          m_rhs_vec += m_vector_scratch;
          
          /// add in scattering moments
          for(int l=0;l<m_n_l_mom;l++)
          {          
            m_sweep_matrix_creator->get_r_sig_s(m_matrix_scratch,l);
            m_rhs_vec += m_ang_quad.get_leg_poly(dir,l)*m_matrix_scratch*m_local_phi[l];
          }
          /// get the local solution
          m_local_soln = m_lhs_mat.partialPivLu().solve(m_rhs_vec);
                    
          /// save the moments of the local solutions (or calculate k_I), and update outflow
          m_sweep_saver->save_local_solution(phi_new,m_local_soln,m_psi_in,cell,grp,dir);
          
        } /// group loop        
      } /// direction loop
    } /// cell loop
  } /// direction set loop (positive vs negative mu) 
    
  m_sweep_type_set = false;
  return;
}

void Transport_Sweep::get_neg_mu_boundary_conditions(const bool is_krylov)
{
  double val = -1.;
  for(int d=0 ; d< m_n_dir/2 ; d++)
  {
    for(int g=0; g< m_n_groups ; g++)
    {
      if(is_krylov)
      {
        val = 0.;
      }
      else
      {
        val = m_sweep_bc_right->get_boundary_condition(m_ang_quad.get_mu(d), g, m_time );
      }
      m_psi_in(g,d) = val;
    }
  }
}

void Transport_Sweep::get_pos_mu_boundary_conditions(const bool is_krylov)
{
  double val = -1.;
  
  for(int d=0 ; d < m_n_dir/2 ; d++)
  {
    for(int g=0; g< m_n_groups ; g++)
    {
      if(is_krylov)
      {
        val = 0.;
      }
      else
      {
        if(m_left_reflecting)
        {
          /// enforce reflection.  Assumes quadrature set is arranged from most negative mu to most positive mu.
          val = m_psi_in(g , m_n_dir/2 - d -1);
        }
        else
        {
          val = m_sweep_bc_left->get_boundary_condition(m_ang_quad.get_mu(d+m_n_dir/2), g , m_time);
        }
      }
      m_psi_in(g,d+m_n_dir/2) = val;
    }
  }
}

void Transport_Sweep::set_time_data( const double dt, const double time_stage, const std::vector<double>& rk_a_of_stage_i, const int stage )
{
  m_time = time_stage;
  m_sweep_matrix_creator->set_time_data(dt,time_stage,rk_a_of_stage_i,stage);
  m_sweep_saver->set_stage(stage);
  return;
}