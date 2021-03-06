#ifndef Time_Data_h
#define Time_Data_h

#include "Input_Reader.h"
#include "Inputs_Allowed.h"

#include "DT_Calculator_Ramp.h"
#include "DT_Calculator_Exponential.h"
#include "DT_Calculator_Vector.h"
#include "DT_Calculator_Temperature_Change.h"
#include <memory>

class Time_Data
{
public:
  /* ***************************************************
  *
  *   Public Functions
  *
  *************************************************** */
  Time_Data(const Input_Reader&  input_reader);
  virtual ~Time_Data(){}
  
  int get_number_of_stages(void) const;
  
  double get_a(const int stage, const int index) const;
  double get_b(const int stage) const;
  double get_c(const int stage) const;
  
  double get_dt(const int step, const double time_now, const double dt_old, const double adapt_criteria);
  
  double get_t_start(void) const;
  double get_t_end(void) const;
  double get_dt_min(void) const;
  double get_dt_max(void) const;
    
protected:  
  int m_number_stages;
  
  const TIME_SOLVER m_time_solver = INVALID_TIME_SOLVER;
  
  const double m_dt_min;
  const double m_dt_max;
  const double m_t_end;
  const double m_t_start;
  /**
   SDIRK time integration basics
   \f{eqnarray}{
        y^{n+1} &=& y^{n} + \Delta t \sum_{s=1}^{N_{stages}}{ b_s k_s} \\
        k_s &=& f(t^n + c_s \Delta t , y_s) \\
        y_s &=& y^n + \Delta t \sum_{j=1}^{s-1}{a_{s,j} k_j } + \Delta t a_{s,s} k_s 
   \f}
  */
  int m_current_dump_goal;
  const int m_n_extra_dumps;
  std::vector<double> m_times_to_dump;
  
  
  /// vector of length n_stages
  std::vector<double> m_b;
  /// vector of length n_stages
  std::vector<double> m_c;
  /// lower triangular matrix (represented as a vector, of size(n_stages*(n_stages+1)/2)
  std::vector<double> m_a;   
  
  void fill_sdirk_vectors(void);
  
/* ***************************************************
  *
  *   Protected Functions
  *
  *************************************************** */
  
  std::shared_ptr<V_DT_Calculator> m_calculate_dt;  
};

#endif