#ifndef Cv_Constant_h
#define Cv_Constant_h

#include "VCv.h"
#include "Input_Reader.h"
#include <vector>
#include <stdlib.h>
#include <math.h>

class Cv_Constant: public VCv
{
public:
  Cv_Constant(const Input_Reader& input_reader, 
  const int mat_num);
  virtual ~Cv_Constant();

  double get_cv(const double position, const double temperature) override;
private:
  /// \f$ C_v = m_{const} \f$
  double m_const = -1.0;
};

#endif
