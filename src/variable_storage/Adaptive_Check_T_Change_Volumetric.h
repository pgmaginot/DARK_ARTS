#ifndef Adaptive_Check_T_Change_Volumetric_h
#define Adaptive_Check_T_Change_Volumetric_h


/** @file   Adaptive_Check_T_Change.h
  *   @author pmaginot
  *   @brief To use when there is no adaptive time check for any variable type
 */
#include "V_Adaptive_Check.h"
#include "Temperature_Data.h"
#include "K_Temperature.h"
#include "Time_Data.h"
#include "Fem_Quadrature.h"
#include "Cell_Data.h"
#include <algorithm>
#include <iomanip>


class Adaptive_Check_T_Change_Volumetric : public V_Adaptive_Check
{
public:
  Adaptive_Check_T_Change_Volumetric(
    const Temperature_Data& t_old , const K_Temperature& k_t , const Time_Data& time_data,
      const Cell_Data& cell_data , const Fem_Quadrature& fem_quadrature, const Input_Reader& input_reader);
    
  virtual ~Adaptive_Check_T_Change_Volumetric(){}

  bool adaptive_check(const double dt, double& adapt_criteria) override;
  
private:  
  const Temperature_Data& m_t_old;
  const K_Temperature& m_k_t;
  const Cell_Data& m_cell_data;
  const Fem_Quadrature& m_fem_quadrature;
  
  const int m_n_stages;
  const int m_n_cells;
  const int m_n_el;
  const int m_n_qp;  
  const int m_n_cell_per_grouping;
  const int m_n_groupings;
  Eigen::VectorXd m_k_t_vec;
  Eigen::VectorXd m_sum_k_t;
  Eigen::VectorXd m_t_old_vec;
  Eigen::VectorXd m_t_new_vec;
  const double m_t_change_goal;
  const double m_offset;
  std::vector<double> m_dfem_at_qp;
  std::vector<double> m_qp_w_at_dfem;
  std::vector<double> m_t_old_evals;
  std::vector<double> m_t_new_evals;
  std::vector<double> m_sdirk_b;
};

#endif