#ifndef Scattering_Opacity_Polynomial_Space_h
#define Scattering_Opacity_Polynomial_Space_h

#include "VScattering_Opacity.h"
#include "Input_Reader.h"
#include <vector>
#include <stdlib.h>
#include <math.h>

class Scattering_Opacity_Polynomial_Space: public VScattering_Opacity
{

public:
  Scattering_Opacity_Polynomial_Space(const Input_Reader& input_reader, 
  const int mat_num);
  virtual ~Scattering_Opacity_Polynomial_Space(){}

  double get_scattering_opacity(const int l_mom, const int group, const double temperature, const double position) override;

private:
  const int m_high_poly_degree;
  std::vector<double> m_poly_coeff;
};

#endif