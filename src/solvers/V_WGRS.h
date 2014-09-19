#ifndef V_WGRS_h
#define V_WGRS_h

#include "Input_Reader.h"
#include "Fem_Quadrature.h"
#include "Cell_Data.h"
#include "Materials.h"
#include "Angular_Quadrature.h"


/** @file   V_WGRS.h
  *   @author pmaginot
  *   @brief provide function interfaces for within group radiation solvers
 */

class V_WGRS
{
public:
  V_WGRS(const Input_Reader& input_reader,
    const Fem_Quadrature& fem_quadrature, Cell_Data* cell_data, Materials* materials, 
    const Angular_Quadrature& angular_quadrature, const int n_stages);
    
  virtual ~V_WGRS(){}

  virtual void solve(const Temperature_Data& t_star, Intensity_Moment_Data& phi) = 0;
protected:

};

#endif