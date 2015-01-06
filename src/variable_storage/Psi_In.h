#ifndef Psi_In_h
#define Psi_In_h

#include <vector>
#include "Dark_Arts_Exception.h"

class Psi_In
{
 public:
 
 Psi_In(const int n_groups, const int n_dir);
 virtual ~Psi_In(){}
 
 double& operator()(int group, int dir);
 
 // void set(const int group, const int dir, const double psi_in);
 
 private:
  
  const int m_n_groups;
  const int m_n_dir;
  
  std::vector<double> m_outflows;
  
  int location(const int group, const int dir) const;
  
};

#endif
