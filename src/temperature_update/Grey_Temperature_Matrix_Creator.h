#ifndef Grey_Temperature_Matrix_Creator_h
#define Grey_Temperature_Matrix_Creator_h

/** @file   Grey_Temperature_Matrix_Creator.h
  *   @author pmaginot
  *   @brief Declare the Temperautre_Update class that will update a Temperature_Object given an Intensity_Object
 */
 #include "V_Temperature_Matrix_Creator.h"

class Grey_Temperature_Matrix_Creator : public V_Temperature_Matrix_Creator
{
public:
  Grey_Temperature_Matrix_Creator(const Fem_Quadrature& fem_quadrature, const Cell_Data& cell_data, Materials& material, 
    const Angular_Quadrature& angular_quadrature, const int n_stages, const Temperature_Data& t_old, const Intensity_Moment_Data& phi);
    
  virtual ~Grey_Temperature_Matrix_Creator(){}
    
  void calculate_update_quantities(const int cell, const Eigen::VectorXd& t_star, const K_Temperature& k_t,
    Eigen::MatrixXd& coefficient , Eigen::VectorXd& rhs) override;
    
  void calculate_k_t(const int cell, const Eigen::VectorXd& t_star, Eigen::VectorXd& k_t) override;
    
protected:
  Eigen::VectorXd scratch1;
  Eigen::VectorXd scratch2;
  Eigen::MatrixXd scratch_mat;
};

#endif