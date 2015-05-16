#include "Adaptive_Check_T_Change.h"

Adaptive_Check_T_Change::Adaptive_Check_T_Change(const Temperature_Data& t_old , const K_Temperature& k_t ,  
  const int n_stages, const int n_cells, const int n_el_cell)
  :
  V_Adaptive_Check() , 
  m_t_old(t_old),
  m_k_t(k_t),
  m_n_stages(n_stages),
  m_n_cells(n_cells),
  m_n_el(n_el_cell)
{
  m_sdirk_b.resize(m_n_stages,0.);

}


bool Adaptive_Check_T_Change::adaptive_check(const int stage, const double dt) 
{
  return false;
}