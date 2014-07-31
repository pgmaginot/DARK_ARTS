#ifndef Absorption_Opacity_Constant_h
#define Absorption_Opacity_Constant_h

#include "VAbsorption_Opacity.h"
#include "Input_Reader.h"
#include <vector>
#include <stdlib.h>
#include <math.h>

class Absorption_Opacity_Constant: public VAbsorption_Opacity
{
public:
  Absorption_Opacity_Constant(const Input_Reader& input_reader, 
  const int mat_num);
  virtual ~Absorption_Opacity_Constant();

  double get_absorption_opacity(const int group, 
    const double temperature, const double position) override;
private:
  /// \f$ \sigma_a = m_{const} \f$
  double m_const = -1.0;
};

#endif