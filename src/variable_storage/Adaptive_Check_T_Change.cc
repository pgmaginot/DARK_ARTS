#include "Adaptive_Check_T_Change.h"

Adaptive_Check_T_Change::Adaptive_Check_T_Change(const Temperature_Data& t_old , const K_Temperature& k_t ,  const Time_Data& time_data,
  const int n_stages, const int n_cells, const int n_el_cell, const double goal)
  :
  V_Adaptive_Check() , 
  m_t_old(t_old),
  m_k_t(k_t),
  m_n_stages(n_stages),
  m_n_cells(n_cells),
  m_n_el(n_el_cell),
  m_k_t_vec(Eigen::VectorXd::Zero(m_n_el) ) , 
  m_t_old_vec(Eigen::VectorXd::Zero(m_n_el) ),
  m_t_change_goal( 1.2*goal )
{
  // std::cout << "Adaptive_Check_T_Change created" << std::endl;
  m_sdirk_b.resize(m_n_stages,0.);

  for(int s=0 ; s< m_n_stages ; s++)
    m_sdirk_b[s] = time_data.get_b(s);
}


bool Adaptive_Check_T_Change::adaptive_check(const double dt, double& adapt_criteria) 
{ 
  adapt_criteria = 0.;
  Eigen::VectorXd sum_k_t = Eigen::VectorXd::Zero(m_n_el);
  
  // std::cout << "dt in check: " << dt << std::endl;
  for(int s=0 ; s< m_n_stages ; s++)
    m_sdirk_b[s] *= dt;
  
  for(int c = 0 ; c < m_n_cells ; c++)
  {
    sum_k_t = Eigen::VectorXd::Zero(m_n_el);
    for(int s = 0 ; s < m_n_stages ; s++)
    {
      m_k_t.get_kt( c, s , m_k_t_vec);
      sum_k_t += m_sdirk_b[s]*m_k_t_vec;
    }
    m_t_old.get_cell_temperature(c,m_t_old_vec);
    /// have t_old, and the difference between t_old and t_next
    for(int el = 0 ; el < m_n_el ; el++)
    {
      double err = fabs( sum_k_t(el) )/( std::max( fabs(m_t_old_vec(el)) +  fabs(m_t_old_vec(el) + sum_k_t(el) )   , 2.0E-3) );
      // std::cout << "err: " << err << std::endl;
      adapt_criteria = std::max(adapt_criteria , err ) ;    
    }
    
  }
  
  adapt_criteria *= 2.;
  
  for(int s=0 ; s< m_n_stages ; s++)
  {
    m_sdirk_b[s] /= dt;
  }
    
  
  // std::cout << "Adapt criteria after being calculated: " << std::scientific << std::setprecision(5) << adapt_criteria << std::endl;
  bool too_big = (adapt_criteria > m_t_change_goal);
  return too_big;
}