/**
  A custom data structure that will hold the result of comparing two Intensity_Moment_Data objects
*/
#ifndef Err_Phi_h
#define Err_Phi_h

class Err_Phi{
public:
  Err_Phi(void);
  ~Err_Phi(){}
  void set_error(int c, int g, int l, double err);
  void clear(void);
  int get_cell_with_worst_err(void) const;
  int get_group_with_worst_err(void) const;
  int get_legendre_moment_with_worst_err(void) const;
  double get_worst_err(void) const;
private:
  double error;
  int cell, group, leg_moment;  
};

#endif
