/**
  A custom data structure that will hold the information relevant to the Temperature_Updates
*/
#ifndef Err_Temperature_h
#define Err_Temperature_h

#include "Eigen/Dense"
#include <vector>

class Err_Temperature{
public:
  Err_Temperature(const int n_el);
  virtual ~Err_Temperature(){}
  void set_error(const int c, const int el, const double err, Eigen::VectorXd& delta);
  void clear(void);
  int get_cell_with_worst_err(void) const;
  int get_element_with_worst_err(void) const;
  double get_worst_err(void) const;
  void get_big_delta(Eigen::VectorXd& vec) const;
  
  
  void set_small_number(const double val);
  double get_small_number(void) const;
  
private:
  const int m_np;
  Eigen::VectorXd big_delta;
  double error;
  int cell, element;  
  
  double small_number;
};

#endif
