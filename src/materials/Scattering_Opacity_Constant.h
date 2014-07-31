#ifndef Scattering_Opacity_Constant_h
#define Scattering_Opacity_Constant_h

#include "VScattering_Opacity.h"
#include "Input_Reader.h"
#include <vector>
#include <stdlib.h>
#include <math.h>

class Scattering_Opacity_Constant: public VScattering_Opacity
{

public:
  Scattering_Opacity_Constant(const Input_Reader& input_reader, 
  const int mat_num);
  virtual ~Scattering_Opacity_Constant();

  double get_scattering_opacity(const int l_mom, const int group, const double temperature, const double position) override;

private:
  /// \f$ \sigma_s = m_{const} \f$
  double m_const = -1.0;
};

#endif