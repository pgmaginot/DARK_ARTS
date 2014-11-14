#ifndef Cv_Rational_h
#define Cv_Rational_h

#include "VCv.h"
#include "Input_Reader.h"
#include <vector>
#include <stdlib.h>
#include <math.h>

class Cv_Rational: public VCv
{
public:
  Cv_Rational(const Input_Reader& input_reader, 
  const int mat_num);
  virtual ~Cv_Rational(){}

  double get_cv(const double position, const double temperature) override;
private:
  /// \f$ C_v = m_{const} \f$
  const double m_const;
  const int m_power;
  const double m_offset;
};

#endif
