#ifndef Adaptive_Check_T_Change_h
#define Adaptive_Check_T_Change_h


/** @file   Adaptive_Check_T_Change.h
  *   @author pmaginot
  *   @brief To use when there is no adaptive time check for any variable type
 */
#include "V_Adaptive_Check.h"
#include "Temperature_Data.h"
#include "K_Temperature.h"


class Adaptive_Check_T_Change : public V_Adaptive_Check
{
public:
  Adaptive_Check_T_Change(
    const Temperature_Data& t_old , const K_Temperature& k_t ,
    const int n_stages, const int n_cells, const int n_el_cell);
    
  virtual ~Adaptive_Check_T_Change(){}

  bool adaptive_check(const int stage, const double dt) override;
  
private:  
  const Temperature_Data& m_t_old;
  const K_Temperature& m_k_t;
  const int m_n_stages;
  const int m_n_cells;
  const int m_n_el;
  std::vector<double> m_sdirk_b;
};

#endif