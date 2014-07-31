#ifndef Scattering_Opacity_Rational_h
#define Scattering_Opacity_Rational_h

#include "VScattering_Opacity.h"
#include "Input_Reader.h"
#include <vector>
#include <stdlib.h>
#include <math.h>

class Scattering_Opacity_Rational: public VScattering_Opacity
{

public:
  Scattering_Opacity_Rational(const Input_Reader& input_reader, 
  const int mat_num);
  virtual ~Scattering_Opacity_Rational();

  double get_scattering_opacity(const int l_mom, const int group, const double temperature, const double position) override;

private:
  /// \f$ \sigma_s = \frac{m_{const} }{m_{offset} + T^{m_p} }  \f$
  double m_const = -1.0;
  double m_offset = -1.0;
  int m_p = -1;
};

#endif