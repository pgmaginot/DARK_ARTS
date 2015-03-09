#ifndef V_Diffusion_Ordering_h
#define V_Diffusion_Ordering_h

#include "Cell_Data.h"
#include "Angular_Quadrature.h"
#include "Dark_Arts_Exception.h"

/** @file   V_Diffusion_Ordering.h
  *   @author pmaginot
  *   @brief Provide a base class that for a given matrix block, determines the cell and group number needed
 */
 

class V_Diffusion_Ordering
{
public:
  /// Only able to initialize if given an Input_Reader object
  /// constructor defined in Fem_Quadrature.cc
  V_Diffusion_Ordering(const Cell_Data& cell_data, const Angular_Quadrature& angular_quadrature, const int n_mip_loops);
    
  virtual ~V_Diffusion_Ordering(){}
    
  virtual void get_cell_and_group(const int block_i , const int mip_loop_number,  int& cell , int& group) = 0;
  
  void check_bounds(const int cell_num , const int group_num);
  
protected:
  /** ****************************************************************
    * Variables that are initialzed in the constructor initialization list
    ****************************************************************
   */   
  const int m_n_groups;
  const int m_n_cells;
  const int m_n_mip_loops;
};

#endif