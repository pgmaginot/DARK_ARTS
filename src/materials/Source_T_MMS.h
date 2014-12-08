#ifndef Source_T_MMS_h
#define Source_T_MMS_h

#include "VSource_T.h"
#include "Input_Reader.h"
#include "MMS_Intensity.h"
#include "MMS_Temperature.h"
#include "VCv.h"
#include "VAbsorption_Opacity.h"
#include "Angular_Quadrature.h"
#include "Planck.h"

class Source_T_MMS : public VSource_T
{

public:
  Source_T_MMS(const Input_Reader& input_reader, 
    const Angular_Quadrature& angular_quadrature,
    std::vector<std::shared_ptr<VAbsorption_Opacity> >& abs_op, 
    std::vector<std::shared_ptr<VCv> >& cv, 
    const int mat_num, const Planck& planck) ;
  virtual ~Source_T_MMS(){}

  double get_temperature_source(const double position, const double time) override;
private:
  const Planck& m_planck;
  const double m_sn_w;
  MMS_Intensity m_intensity;
  MMS_Temperature m_temperature;
    
  std::shared_ptr<VAbsorption_Opacity> m_abs_op;
  std::shared_ptr<VCv> m_cv;
};

#endif