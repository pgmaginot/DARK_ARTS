#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <Eigen/Dense>

#include "Dark_Arts_Exception.h"

int main(int argc, char** argv)
{
  int val = 0;
  
  try{
    const int n_p = 4;
    std::vector< std::vector<double> > values_to_save(n_p, std::vector<double>(n_p,0.));
    
    Eigen::MatrixXd col_major_mat = Eigen::MatrixXd::Zero(n_p,n_p);
    for(int i=0; i <n_p ; i++)
    {
      for(int j=0; j<n_p ;j++)
      {
        values_to_save[i][j] = rand();
        col_major_mat(i,j) = values_to_save[i][j];
      }
    }
    
    std::cout << "This is col_major_mat:\n" << col_major_mat;
    std::cout << "\ncol_major_mat is: " << col_major_mat.IsRowMajor << " with respect to being row major\n";
    
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> row_major_mat(n_p,n_p);
    std::cout << "This is the unitialized row_major_mat :\n" << row_major_mat;
    std::cout << "\n row_major_mat is: " << row_major_mat.IsRowMajor << " with respect to being row major\n";
    
    row_major_mat = Eigen::MatrixXd::Zero(n_p,n_p);
    std::cout << "Now I've zeroedrow_major_mat with MatrixXd.  \n row_major_mat is: " << row_major_mat.IsRowMajor << " with respect to being row major\n";  
      
    row_major_mat = col_major_mat;
    std::cout << "Now I've set row_major_mat equal to a col_major_mat matrix.  \n row_major_mat is: " << row_major_mat.IsRowMajor << " with respect to being row major\n";  
    
    std::cout << "This is row_major_mat now:\n" << row_major_mat << std::endl;
    
    double *col_maj_ptr = &col_major_mat(0,0);
    double *row_maj_ptr = &row_major_mat(0,0);
    
    for(int i = 0 ; i < n_p*n_p ; i++)
    {
      std::cout << "row_maj_ptr["<<i<<"]: " << row_maj_ptr[i] << std::endl;
    }
    std::cout << std::endl;
    for(int i = 0 ; i < n_p*n_p ; i++)
    {
      std::cout << "col_maj_ptr["<<i<<"]: " << col_maj_ptr[i] << std::endl;
    }
  }
  catch(const Dark_Arts_Exception& da)
  {
    val = -1;
    da.testing_message();
  }
  
  return val;
}
