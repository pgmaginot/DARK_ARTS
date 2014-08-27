// Copyright (c) 2000-2008, Texas Engineering Experiment Station (TEES), a
// component of the Texas A&M University System.
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are not permitted without specific prior written permission
// from TEES.

// If written permission is obtained for redistribution or further use, the
// following conditions must be met:

// 1) Redistributions of source code must retain the above copyright notice,
// this list of conditions and the disclaimer below.

// 2) Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions, and the disclaimer below in the documentation and/or
// other materials provided with the distribution.

// 3) Neither the name of TEES, the name of the Texas A&M University System, nor
// the names of its contributors may be used to endorse or promote products
// derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS AS IS
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.


//==============================================================================
//
// Planck.h
//
// 10/22/08 MPA
//
// Header file for Planck class.
//
// 06/17/2009 WDH
//
// Added integrate_B_grey() and integrate_dBdT_grey() methods to correctly
// return acT^4 and 4acT^3 respectfully if the problem is grey.  In the future,
// we need to modify integrate_B() and integrate_dBdT() to return the correct
// values for a grey problem.
//
//==============================================================================


#include <vector>
#include <cassert>
#include <exception>
#include <sstream>

#include <cmath>
#include <iostream>

// #ifndef math_h
// #define math_h
// #include "math.h"
// #endif
// #ifndef stdio_h
// #define stdio_h
// #include "stdio.h"
// #endif

//#include "CommonException.h"



using std::vector;
using std::stringstream;

#ifndef Planck_h
#define Planck_h

class Planck
{

private:
  // accuracy parameter
  double accuracy;

  // physical constants, intialized in constructor
  double h;   // Planck constant
  double k;   // Boltzmann constant
  double c;   // speed of light
  double a;   // radiation constant a
  double pi;  // Pi


  // points and weights are for integration by gaussian quadrature
  vector<long double> points;
  vector<long double> weights;

  // sets guassian quadrature for integration
  void gauss_quad();

public:

  Planck(double accuracy_parameter = 1e-15);//  {gauss_quad(accuracy)};
  ~Planck() {}; // empty destructor

  // to find the planck function B at a specific temperature and energy
  double get_B(double T, double E);

  // returns the partial of B with respect to temperature
  double get_dBdT(double T, double E);

  // to find Bg; pass temp, e_min, e_max
  double integrate_B(double T, double E_min, double E_max);

  //Get Temp for corresponding Energy
  double getTempForEnergy(double E);

  // grey-case Bg
  double integrate_B_grey(double T);

  // to find (dB/dT)g
  double integrate_dBdT(double T, double E_min, double E_max);

  // grey-case (dB/dT)g
  double integrate_dBdT_grey(double T);

};
inline double Planck::get_B(double T, double E)
{
  if( T == 0 )
    return 0;
  if( T < 0 || E < 0 )
  {
    stringstream str;
    str << "Planckian evaluation requires temperature and energy to be "
        << "positive.";
    //throw CommonException(str, CET_INTERNAL_ERROR);
  }

  //                     2 E^3
  //    B(E, T) = ----------------------
  //              h^3 c^2 (e^(E/kT) - 1)

  double tempVal= 2. * pow(E, 3) * pow(h, -3) * pow(c, -2) / (exp(E/(k*T)) - 1.);
  //double tempVal1=2 * pow(E, 3) * pow(h, -3) * pow(c, -2) / (exp(E/(k*T)) - 1);
  return tempVal;

} // Planck::get_B()

inline double Planck::get_dBdT(double T, double E)
{
  if( T == 0 )
    return 0;
  if( T < 0 || E < 0 )
  {
    stringstream str;
    str << "Planckian evaluation requires temperature and energy to be "
        << "positive.";
    //throw CommonException(str, CET_INTERNAL_ERROR);
  }

  //    dB(E, T)       2      E^4      e^(E/kT)
  //    -------- = ---------  ---  ----------------
  //       dT      h^3 c^2 k  T^2  (e^(E/kT) - 1)^2

  return 2 * pow(h,-3) * pow(c,-2) * pow(k,-1)
    * pow(E, 4) * pow(T,-2) * exp(E/(k*T)) * pow(exp(E/(k*T)) - 1, -2);

} // Planck::get_dBdT()

#endif