#include "Adaptive_Check_T_Change_Volumetric.h"

Adaptive_Check_T_Change_Volumetric::Adaptive_Check_T_Change_Volumetric(
    const Temperature_Data& t_old , const K_Temperature& k_t , const Time_Data& time_data,
      const Cell_Data& cell_data , const Fem_Quadrature& fem_quadrature, const Input_Reader& input_reader)
  :
  V_Adaptive_Check() , 
  m_t_old(t_old),
  m_k_t(k_t),
  m_cell_data(cell_data),
  m_fem_quadrature(fem_quadrature) ,
  m_n_stages(time_data.get_number_of_stages() ),
  m_n_cells(m_cell_data.get_total_number_of_cells() ),
  m_n_el(fem_quadrature.get_number_of_interpolation_points() ),
  m_n_qp( fem_quadrature.get_number_of_source_points() ),
  m_n_cell_per_grouping( input_reader.get_cells_per_adaptive_grouping() ),
  m_n_groupings( (m_n_cells/m_n_cell_per_grouping) ) ,  
  m_k_t_vec(Eigen::VectorXd::Zero(m_n_el) ) , 
  m_sum_k_t(Eigen::VectorXd::Zero(m_n_el) ) , 
  m_t_old_vec(Eigen::VectorXd::Zero(m_n_el) ),
  m_t_new_vec(Eigen::VectorXd::Zero(m_n_el) ),
  m_t_change_goal( 1.2*input_reader.get_t_change_adaptive_goal() ),
  m_offset( input_reader.get_t_change_floor() )
{
  if( (m_n_cells%m_n_cell_per_grouping) != 0)
    throw Dark_Arts_Exception(INPUT, "Number of cells must be cleanly divided by number of cells per grouping for T_Change_Volumetric adaptivity");
  
  m_dfem_at_qp.resize(m_n_qp*m_n_el,0.);
  m_qp_w_at_dfem.resize(m_n_qp,0.);
  
  fem_quadrature.get_dfem_at_source_points( m_dfem_at_qp );
  fem_quadrature.get_source_weights( m_qp_w_at_dfem );
  
  m_sdirk_b.resize(m_n_stages,0.);

  for(int s=0 ; s< m_n_stages ; s++)
    m_sdirk_b[s] = time_data.get_b(s);
}


bool Adaptive_Check_T_Change_Volumetric::adaptive_check(const double dt, double& adapt_criteria) 
{ 
  adapt_criteria = 0.;
  
  for(int s=0 ; s< m_n_stages ; s++)
  {
    m_sdirk_b[s] *= dt;
  }
  
  for(int grp = 0 ; grp < m_n_groupings ; grp++)
  {
    double delta_grp = 0.;
    double denom_grp = 0.;
    for(int c = 0 ; c < m_n_cell_per_grouping ; c++)
    {
      
      double delta_cell = 0.;
      double denom_cell = 0.;
      int cell = grp*m_n_cell_per_grouping + c;
      
      m_sum_k_t = Eigen::VectorXd::Zero(m_n_el);
      for(int s = 0 ; s < m_n_stages ; s++)
      {
        m_k_t.get_kt( cell, s , m_k_t_vec);
        m_sum_k_t += m_sdirk_b[s]*m_k_t_vec;
      }
      m_t_old.get_cell_temperature(cell,m_t_old_vec);
      m_t_new_vec = m_t_old_vec + m_sum_k_t;
      /// have t_old, and the difference between t_old and t_next
      /// now we need to calculate the L2 norm contribution of both over this grouping
      m_fem_quadrature.evaluate_variable_at_quadrature_pts(m_t_old_vec , m_dfem_at_qp , m_t_old_evals);
      m_fem_quadrature.evaluate_variable_at_quadrature_pts(m_t_new_vec , m_dfem_at_qp , m_t_new_evals);
      
      double dx = m_cell_data.get_cell_width(cell);
      for(int q = 0 ; q < m_n_qp ; q++)
      {
        delta_cell += m_qp_w_at_dfem[q]*( m_t_new_evals[q] - m_t_old_evals[q])*( m_t_new_evals[q] - m_t_old_evals[q] )  ;
        denom_cell += m_qp_w_at_dfem[q] * ( m_t_new_evals[q] + m_t_old_evals[q] + m_offset)*( m_t_new_evals[q] + m_t_old_evals[q] + m_offset);
      }
      delta_cell *= dx/2.;
      denom_cell *= dx/2.;
      
      delta_grp += delta_cell;
      denom_grp += denom_cell;
    }    
    
    double err = sqrt(delta_grp)/sqrt(denom_grp);
    adapt_criteria = std::max(adapt_criteria , err ) ;   
  }
    
  for(int s=0 ; s< m_n_stages ; s++)
  {
    m_sdirk_b[s] /= dt;
  }
  bool too_big = (adapt_criteria > m_t_change_goal);
  return too_big;
}