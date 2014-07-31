#ifndef Materials_h
#define Materials_h

#include <vector>
#include <memory>
#include <stdlib.h>

#include "Input_Reader.h"
#include "VAbsorption_Opacity.h"
#include "VScattering_Opacity.h"
#include "VCv.h"
#include "VSource_T.h"
#include "VSource_I.h"
#include "Absorption_Opacity_Rational.h"
#include "Scattering_Opacity_Rational.h"
#include "Absorption_Opacity_Constant.h"
#include "Scattering_Opacity_Constant.h"
#include "Absorption_Opacity_Table.h"
#include "Scattering_Opacity_Table.h"
#include "Source_I_Constant.h"
#include "Source_T_Constant.h"
#include "Cv_Constant.h"

class Materials
{

public:
  Materials(const Input_Reader& input_reader);
  ~Materials(){}

  void load_materials(const Input_Reader& input_reader);
  /**
    * Want to store opacity, cv, source objects
    * Return (virtual) base class pointers to give outsiders access to material properties
  */  
private:
  /** Vectors of pointers to the various opacity, cv, and source objects we will create in
   * the Materials class constructor (based off what we read in the input file)
  */
  std::vector<std::shared_ptr<VAbsorption_Opacity>> m_abs_opacities;
  std::vector<std::shared_ptr<VScattering_Opacity>> m_scat_opacities;
  std::vector<std::shared_ptr<VCv>> m_cv;
  std::vector<std::shared_ptr<VSource_T>> m_source_t;
  std::vector<std::shared_ptr<VSource_I>> m_source_i;
  
  int m_num_materials = -1;
};

#endif