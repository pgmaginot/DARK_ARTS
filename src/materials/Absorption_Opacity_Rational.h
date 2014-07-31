#ifndef Absorption_Opacity_Rational_h
#define Absorption_Opacity_Rational_h

#include "VAbsorption_Opacity.h"
#include "Input_Reader.h"
#include <vector>
#include <stdlib.h>
#include <math.h>

class Absorption_Opacity_Rational: public VAbsorption_Opacity
{
public:
  Absorption_Opacity_Rational(const Input_Reader& input_reader, 
  const int mat_num);
  virtual ~Absorption_Opacity_Rational();

  double get_absorption_opacity(const int group, const double temperature, const double position) override;
private:
  /// \f$ \sigma_a = \frac{m_{const} }{m_{offset} + T^{m_p} }  \f$
  double m_const = -1.0;
  double m_offset = -1.0;
  int m_p = -1;
};

#endif