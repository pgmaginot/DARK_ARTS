#ifndef Psi_In_h
#define Psi_In_h

#include <vector>

class Psi_In
{
 public:
 
 Psi_In(const int n_groups, const int n_dir);
 ~Psi_In(){}
 
 double operator()(int group, int dir) const;
 
 void set(const int group, const int dir, const double psi_in);
 
 private:
  
  const int m_n_groups;
  const int m_n_dir;
  
  std::vector<double> m_outflows;
  
  int location(const int group, const int dir) const;
  
};

#endif
