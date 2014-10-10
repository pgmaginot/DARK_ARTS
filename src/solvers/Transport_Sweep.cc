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
  const Temperature_Data& t_star)
  :
  m_n_cells{ cell_data.get_total_number_of_cells() },
  m_n_groups{ angular_quadrature.get_number_of_groups()  },
  m_n_dir{ angular_quadrature.get_number_of_dir() },
  m_n_l_mom{ angular_quadrature.get_number_of_leg_moments() },
  m_np{ fem_quadrature.get_number_of_interpolation_points() }, 
  m_ang_quad(angular_quadrature),
  m_matrix_scratch{ Eigen::MatrixXd(m_np,m_np) },
  m_vector_scratch{ Eigen::VectorXd::Zero(m_np) },
  m_local_phi(m_n_l_mom, Eigen::VectorXd::Zero(m_np) ),
  m_rhs_vec{ Eigen::VectorXd::Zero(m_np) },
  m_lhs_mat{ Eigen::MatrixXd(m_np,m_np) },
  m_local_soln{ Eigen::VectorXd::Zero(m_np)  },
  m_time{-1.},
  m_psi_in(m_n_groups,m_n_dir)
{
  if(m_n_groups > 1)
  {
    m_sweep_matrix_creator = std::shared_ptr<V_Sweep_Matrix_Creator> 
    (new Sweep_Matrix_Creator_MF(fem_quadrature, materials, n_stages, angular_quadrature.get_sum_w() , m_n_l_mom,
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
  m_k_i_saver = std::shared_ptr<V_Solution_Saver> (new Solution_Saver_K_I(fem_quadrature,m_sweep_matrix_creator,angular_quadrature,ki) );
  m_angle_integrated_saver = std::shared_ptr<V_Solution_Saver> (new Solution_Saver_Flux_Moments(fem_quadrature,angular_quadrature) );
}

void Transport_Sweep::set_ard_phi_ptr(Intensity_Moment_Data* ard_phi_ptr)
{
  m_sweep_matrix_creator->set_ard_phi_ptr(ard_phi_ptr);
  return;
}

void Transport_Sweep::sweep_mesh(const Intensity_Moment_Data& phi_old, Intensity_Moment_Data& phi_new, const bool is_krylov, const bool is_get_k_i)
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
  if(is_krylov)
  {
    m_sweep_source = m_no_source;
  }
  else
  {
    m_sweep_source = m_fixed_source;
  }
  
  if(is_get_k_i)
  {
    if(is_krylov)
    {
      std::cerr << "Cannot calculate k_i with a Krylov mode sweep\n";
      exit(EXIT_FAILURE);
    }
    m_sweep_saver = m_k_i_saver;
  }
  else
  {
    m_sweep_saver = m_angle_integrated_saver; 
    phi_new.clear_angle_integrated_intensity();
  }
  
  
  
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
    }
    else
    {
      cell_start = 0;
      cell_end = m_n_cells - 1;
      incr = 1;
      d_offset = m_n_dir/2;
    }
    
    /// load boundary conditions / inflow intensities
    get_boundary_conditions(is_krylov);
    
    for( int cell = cell_start ; cell != (cell_end+incr) ; cell += incr)
    {
      m_sweep_matrix_creator->update_cell_dependencies(cell);
      for(int grp = 0; grp < m_n_groups ; grp++)
      {
        /// avoid calculating cross sections as much as possible, call m_sweep_matrix_creator->update_group_dependencies
        m_sweep_matrix_creator->update_group_dependencies(grp);
        /// populate V_Sweep_Matrix_Creator m_k_i_r_sig_t and m_k_i_r_sig_s_zero
        if(is_get_k_i)
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
          m_rhs_vec *= m_psi_in(dir,grp);          
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
  
  return;
}

void Transport_Sweep::get_boundary_conditions(const bool is_krylov)
{
  double val = -1.;
  for(int d=0 ; d< m_n_dir ; d++)
  {
    for(int g=0; g< m_n_groups ; g++)
    {
      if(is_krylov)
      {
        val = m_ang_quad.calculate_boundary_conditions(d, g, m_time);
      }
      else
      {
        val = 0.;
      }
      m_psi_in(g,d) = val;
    }
  }
}

void Transport_Sweep::set_time_data( const double dt, const double time_stage, const std::vector<double>& rk_a_of_stage_i, const int stage )
{
  m_time = time_stage;
  m_sweep_matrix_creator->set_time_data(dt,time_stage,rk_a_of_stage_i,stage);
  return;
}