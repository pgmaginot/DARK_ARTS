#ifndef Absorption_Opacity_Table_h
#define Absorption_Opacity_Table_h

#include "VAbsorption_Opacity.h"
#include "Input_Reader.h"
#include <vector>
#include <stdlib.h>
#include <math.h>

class Absorption_Opacity_Table: public VAbsorption_Opacity
{
public:
  Absorption_Opacity_Table(const Input_Reader& input_reader, 
  const int mat_num);
  virtual ~Absorption_Opacity_Table();

  double get_absorption_opacity(const int group, 
    const double temperature, const double position) override;
private:  
  std::string m_filename;
};

#endif