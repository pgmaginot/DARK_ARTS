#ifndef Absorption_Opacity_Adams_Novak_Model_h
#define Absorption_Opacity_Adams_Novak_Model_h

#include "VAbsorption_Opacity.h"
#include "Input_Reader.h"
#include <vector>
#include <stdlib.h>
#include <math.h>

class Absorption_Opacity_Adams_Novak_Model: public VAbsorption_Opacity
{
public:
  Absorption_Opacity_Adams_Novak_Model(const Input_Reader& input_reader);
  virtual ~Absorption_Opacity_Adams_Novak_Model();

  double get_absorption_opacity(const int group, const double temperature, const double position) override;
private:

};

#endif