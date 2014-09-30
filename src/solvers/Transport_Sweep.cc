#include "Transport_Sweep.h"

Transport_Sweep::Transport_Sweep(const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* materials, 
  const Angular_Quadrature& angular_quadrature, const int n_stages)
  :
  m_n_cells{ cell_data->get_total_number_of_cells() },
  m_n_groups{ angular_quadrature.get_number_of_groups()  },
  m_n_dir{ angular_quadrature.get_number_of_dir() },
  m_n_l_mom{ angular_quadrature.get_number_of_leg_moments() },
  m_np{ fem_quadrature.get_number_of_interpolation_points() }, 
  m_sum_w{ angular_quadrature.get_sum_w() },
  m_ang_quad{&angular_quadrature},
  m_matrix_scratch{ Eigen::MatrixXd(m_np,m_np) },
  m_rhs_vec{ Eigen::VectorXd::Zero(m_np) },
  m_lhs_mat{ Eigen::MatrixXd(m_np,m_np) },
  m_local_soln{ Eigen::VectorXd::Zero(m_np)  }
{
  if(m_n_groups > 1)
  {
    m_sweep_matrix_creator = std::shared_ptr<V_Sweep_Matrix_Creator> (new Sweep_Matrix_Creator_MF(fem_quadrature, materials, n_stages, angular_quadrature.get_sum_w() ) );
  }
  else
  {
    m_sweep_matrix_creator = std::shared_ptr<V_Sweep_Matrix_Creator> (new Sweep_Matrix_Creator_Grey(fem_quadrature, materials, n_stages, angular_quadrature.get_sum_w()) );
  }
 
}

void Transport_Sweep::set_ard_phi_ptr(Intensity_Moment_Data* ard_phi_ptr)
{
  m_sweep_matrix_creator->set_ard_phi_ptr(ard_phi_ptr);
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
  
  /// starting and ending number of cells
  int cell_end, cell_start;
  /// increment or deceremnt depending on direction set
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
    for( int cell = cell_start ; cell != (cell_end+incr) ; cell += incr)
    {
      for(int d = 0; d<(m_n_dir/2) ; d++)
      {
        dir = d_offset + d;
        // mu = 
        
        for(int grp = 0; grp < m_n_groups ; grp++)
        {
        
        } /// group loop        
      } /// direction loop
    } /// cell loop
  } /// direction set loop (positive vs negative mu)
  
  
  
  
  return;
}

