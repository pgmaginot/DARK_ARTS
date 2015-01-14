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

/*!
  \file global_constants.h
  \brief Constants and simple manipulation routines used
  throughout the application.
*/

#ifndef _global_constants_h
#define _global_constants_h

const unsigned int MAX_OUTPUT_BUFFER_SIZE     = 8388604; // 8MB
const unsigned int DEFAULT_OUTPUT_BUFFER_SIZE = 2097151; // 2MB
const unsigned int IO_CORE_LIMIT              = 16;

/*!
  \brief Some universal constants (primarily required for radiation transport)
*/
const double PLANCK_CONSTANT      = 4.13566733e-15;         // eV s
const double BOLTZMANN_CONSTANT   = 8.617343e-5;            // eV K^-1
const double SPEED_OF_LIGHT       = 299792458.0*100;        // cm s^-1
const double PI                   = 3.14159265358979323846;
const double FOUR_PI              = 4.0*PI;
const double RADIATION_CONSTANT_A = 4.722181875e-03;        // eV/cm^3-K^4
const double MEV2J                = 1.602176487e-13;        // J/MeV
const double EV2J                 = 1.602176487e-19;        // J/eV
const double NATURAL_LOG_2        = 0.693147181;            // natural log of 2
const double AVOGADROS_NUMBER     = 6.02214129e23;          // atoms/mol
const double AVOGADROS_NUMBER_SM  = 6.02214129e-01;         // atoms-bn/cm-mol
const double BN_PER_CM2           = 1.0e-24;                // bn/cm^2

// RADIATION_CONSTANT_A = (8.0*pow(PI,5)*pow(BOLTZMANN_CONSTANT,4))/
//  (15.0*pow(PLANCK_CONSTANT,3)*pow(SPEED_OF_LIGHT,3));   // eV/cm^3-K^4

#endif
