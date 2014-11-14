#ifndef Absorption_Opacity_Polynomial_Space_h
#define Absorption_Opacity_Polynomial_Space_h

#include "VAbsorption_Opacity.h"
#include "Input_Reader.h"
#include <vector>
#include <stdlib.h>
#include <math.h>

class Absorption_Opacity_Polynomial_Space: public VAbsorption_Opacity
{

public:
  Absorption_Opacity_Polynomial_Space(const Input_Reader& input_reader, 
  const int mat_num);
  virtual ~Absorption_Opacity_Polynomial_Space(){}

  double get_absorption_opacity(const int group, 
    const double temperature, const double position) override;

private:
  const int m_high_poly_degree;
  std::vector<double> m_poly_coeff;
};

#endif