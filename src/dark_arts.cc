#include "Input_Reader.h"
#include "Fem_Quadrature.h"
#include "Cell_Data.h"
#include "Time_Stepper.h"
#include "Angular_Quadrature.h"
#include "Intensity_Data.h"
#include "Materials.h"

#include "Eigen/Dense"

int main(int argc, char** argv)
{
  std::cout << "argc = " << argc << '\n'; 
  for(int i = 0; i < argc; i++) 
    std::cout << "argv[" << i << "] = " << argv[i] << '\n'; 
  
  Input_Reader input_reader;
  bool input_parsed = false;
    
  input_parsed = input_reader.read_xml(argv[1]);
    
  if(!input_parsed)
  {
    std::cerr << "Error reading input" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  /// Initialize a Quadrule object to be able to get all of the quadrature we need
  Quadrule_New quad_fun;
  
  /// Initialize FEM data
  /// get all interpolation points, quadrature formuals, matrix formation routines, etc.
  Fem_Quadrature fem_quadrature( input_reader , quad_fun);
    
  /// Initalize cell data (dx, xL, xR, x_ip, material_number)
  Cell_Data cell_data( input_reader );
  
  /// Initialize time-stepping scheme (SDIRK method)
  Time_Stepper time_stepper( input_reader );
  
  /// Initialize angular quadrature data.  Will include number of: directions, groups, and legendre moments.
  /// will also include evaluations of Legendre polynomials
  Angular_Quadrature angular_quadrature( input_reader , quad_fun );
  
  /// Initialize intensity and angle integrated intensity of previous time step
  Intensity_Data intensity_old( cell_data, angular_quadrature, fem_quadrature);
  
  /// Create a Materials object that contains all opacity, heat capacity, and source objects
  Materials materials( input_reader, fem_quadrature , &cell_data);
  
  
  std::shared_ptr<V_Matrix_Construction>
  Eigen::MatrixXd result(3,3);  
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> diag1(3);
  Eigen::Matrix3d base;
  
  base(0,0) = 11.;
  base(0,1) = 12.;
  base(0,2) = 13.;
  base(1,0) = 21.;
  base(1,1) = 22.;
  base(1,2) = 23.;
  base(2,0) = 31.;
  base(2,1) = 32.;
  base(2,2) = 33.;
  
  diag1.diagonal()[0] = 0.1;
  diag1.diagonal()[1] = 0.2;
  diag1.diagonal()[2] = 0.3;
  
  std::cout << "base matrix= "<< std::endl << base << std::endl;  
  
  std::cout << "diagonal matrix="<< std::endl << diag1.diagonal() << std::endl;
  
  
  result = base * diag1;
  std::cout << "base*diag1" << std::endl << result << std::endl;
  
  result = diag1*base;
  std::cout << "diag1*base" << std::endl << result << std::endl;
  
  // Eigen::Matrix3d dense_diagonal;
  // dense_diagonal(0,0) = 0.1;
  // dense_diagonal(1,1) = 0.2;
  // dense_diagonal(2,2) = 0.3;
  
  // result = base * dense_diagonal;
  // std::cout << "base*dense_diagonal" << std::endl << result << std::endl;
  
  // result = dense_diagonal*base;
  // std::cout << "dense_diagonal*base" << std::endl << result << std::endl;
  
  // result = base + dense_diagonal;
  // std::cout << "base+dense_diagonal" << std::endl << result << std::endl;
  
  // result = dense_diagonal + base;
  // std::cout << "dense_diagonal+base" << std::endl << result << std::endl;
  
  // std::cout << "Diagonal of dense_diagonal" << std::endl << dense_diagonal.diagonal() << std::endl;
  
  // Eigen::VectorXd eig_type_vec(3);
  // eig_type_vec(0) = 0.01;
  // eig_type_vec(1) = 0.02;
  // eig_type_vec(2) = 0.02;
  
  
  
  // std::cout << "Return of function" << std::endl << mat
  
 
  

  
  
  
  // std::vector<double> dbl_vec(3,0.);
  // dbl_vec[0] = 0.01;
  // dbl_vec[1] = 0.02;
  // dbl_vec[2] = 0.03;
  
  
  
  return 0;
}