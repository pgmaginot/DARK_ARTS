#ifndef Cv_Polynomial_Temperature_h
#define Cv_Polynomial_Temperature_h

#include "VCv.h"
#include "Input_Reader.h"
#include <vector>
#include <stdlib.h>
#include <math.h>

class Cv_Polynomial_Temperature: public VCv
{

public:
  Cv_Polynomial_Temperature(const Input_Reader& input_reader, 
  const int mat_num);
  virtual ~Cv_Polynomial_Temperature(){}

  double get_cv(const double position, const double temperature) override;

private:
  const int m_high_poly_degree;
  std::vector<double> m_poly_coeff;
};

#endif