#ifndef MIP_Kappa_Calculator_h
#define MIP_Kappa_Calculator_h

#include <memory>
#include "Input_Reader.h" 

class MIP_Kappa_Calculator
{
public:
  /// Only able to initialize if given an Input_Reader object
  /// constructor defined in Fem_Quadrature.cc
  MIP_Kappa_Calculator(const int p_ord , const double z_mip);
    
  virtual ~MIP_Kappa_Calculator(){}
      
  double calculate_interior_edge_kappa(const double dx_l, const double dx_r , const double  d_l , const double d_r) const;
  
  double calculate_boundary_kappa(const double dx , const double d) const;
    
protected:
  /** ****************************************************************
    * Variables that are initialzed in the constructor initialization list
    ****************************************************************
   */   
  const double m_ord;
  const double m_z_mip;    
};

#endif