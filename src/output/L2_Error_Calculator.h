#ifndef L2_Error_Calculator_h
#define L2_Error_Calculator_h

#include "Input_Reader.h"
#include "Temperature_Data.h"
#include "Intensity_Data.h"
#include "Intensity_Moment_Data.h"
#include "Dark_Arts_Exception.h"

/** @file   L2_Error_Calculator.h
  *   @author pmaginot
  *   @brief Declare Input_Reader class
  *   @class Input_Reader
 */
enum Data_Type{ PHI , INTENSITY, TEMPERATURE}; 
 
class L2_Error_Calculator
{
public:
  L2_Error_Calculator(const Angular_Quadrature& angular_quadrature,
    const Fem_Quadrature& fem_quadrature, 
    const Cell_Data& cell_data, 
    const Input_Reader& input_reader);
    
  virtual ~L2_Error_Calculator(){}
  
  double calculate_l2_error(const double time_eval, const Intensity_Moment_Data& phi) const;
  double calculate_l2_error(const double time_eval, const Intensity_Data& intensity) const;
  double calculate_l2_error(const double time_eval, const Temperature_Data& temperature) const;
  
  double calculate_cell_avg_error(const double time_eval, const Intensity_Moment_Data& phi) const;
  double calculate_cell_avg_error(const double time_eval, const Intensity_Data& intensity) const;
  double calculate_cell_avg_error(const double time_eval, const Temperature_Data& temperature) const;
  
protected:
  /// for temperature and phi
  double calculate_local_l2_error(const Eigen::VectorXd& numeric_at_dfem, const double xL , 
    const double dx , const double time, const Data_Type data ) const;  
  double calculate_local_cell_avg_error(const Eigen::VectorXd& numeric_at_dfem , const double xL , 
    const double dx , const double time,const Data_Type data ) const;
  
  /// for intensity
  double calculate_local_l2_error(const Eigen::VectorXd& numeric_at_dfem, const double xL , 
    const double dx , const double time, const int dir, const Data_Type data ) const;  
  double calculate_local_cell_avg_error(const Eigen::VectorXd& numeric_at_dfem , const double xL , 
    const double dx , const double time, const int dir, const Data_Type data ) const;
  
  const int m_n_dfem;
  const int m_n_cells;
  const int m_n_dir;
  
  const int m_n_quad_pts;
  
  const Cell_Data& m_cell_data;
  const Angular_Quadrature& m_angular_quadrature;
  const Fem_Quadrature& m_fem_quadrature;
  MMS_Temperature m_analytic_temperature;
  MMS_Intensity m_analytic_intensity;
  
  std::vector<double> m_integration_quadrature_pts;
  std::vector<double> m_integration_quadrature_wts;
  
  std::vector<double> m_dfem_at_integration_pts;
};


#endif