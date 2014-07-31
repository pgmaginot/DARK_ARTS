#ifndef Scattering_Opacity_Table_h
#define Scattering_Opacity_Table_h

#include "VScattering_Opacity.h"
#include "Input_Reader.h"
#include <vector>
#include <stdlib.h>
#include <math.h>

class Scattering_Opacity_Table: public VScattering_Opacity
{
public:
  Scattering_Opacity_Table(const Input_Reader& input_reader, 
  const int mat_num);
  virtual ~Scattering_Opacity_Table();

  double get_scattering_opacity(const int l_mom, const int group, 
    const double temperature, const double position) override;
private:  
  std::string m_filename;
};

#endif