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

#ifndef opacities_h
#define opacities_h

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cassert>
#include <bitset>
#include <exception>
#include <stdexcept>

#include "global_constants.h"
#include "Planck.h"
#include "Edits.h"
#include "CommonException.h"
#include "basic_macros.h"
#include "Kind_Info.h"
#include "Logger.h"
#include <stapl/runtime.hpp>

// =============================================================================
//    enums of file format and opacity type
// =============================================================================

#define OPACITY_MAX 1.0e20

enum OPAC_STATUS { OS_VALID, OS_LOW_T, OS_HIGH_T, OS_LOW_RHO, OS_HIGH_RHO };


enum opac_type
{
  OT_grey_ross              , // grey rosseland averaged opacities
  OT_grey_planck            , // grey planck averaged
  OT_mg_ross                , // multigroup rosseland
  OT_mg_planck              , // multigroup planck
  OT_model_ross             , // rosseland averaged 'model' or analytic opacity
  OT_model_planck           , // planck averaged 'model' or analytic opacity
  OT_abc                    , // sigma = A + B * T ^ C

  OT_COUNT                    // a way to size a bitset
  // *** DO NOT ADD TYPES AFTER OT_COUNT !!!
};

// =============================================================================
//    Opacities class
// =============================================================================

class Opacities
{
protected:

  // name for the object
  std::string id;

  // units, defined as an enum in Common/Kind_Info.h
  UNIT_TYPE opac_units;
  // true if microscopic, false if macroscopic
  bool microscopic;

  // a way for the object to know what info it's maintaining
  std::bitset<OT_COUNT> opac_store;
  opac_type OT_default;

  //    Opacity data storage members
  // ==================================
  real8 mod_const;
  std::vector<std::vector<real8> > grey_planck;
  std::vector<std::vector<real8> > grey_ross;
  std::vector<std::vector<std::vector<real8> > > planck;
  std::vector<std::vector<std::vector<real8> > > ross;
  std::vector<real8> abc;

  //    Independent variable data members
  // =======================================
  std::vector<real8> temps;
  std::vector<real8> densities;
  std::vector<real8> group_structure;
  size_t num_temps;
  size_t num_dens;
  size_t num_groups;

  void resize(std::vector<std::vector<real8> >& vec, size_t a, size_t b)
  {
    vec.clear();  vec.resize(a);
    for( size_t i=0; i!=a; ++i )
      vec[i].resize(b, 0.0);
  }
  void resize(std::vector<std::vector<std::vector<real8> > >& v, size_t a,
    size_t b,size_t c)
  {
    v.clear();  v.resize(a);
    for( size_t i=0; i!=a; ++i )
    {
      v[i].resize(b);
      for( size_t j=0; j!=b; ++j )  v[i][j].resize(c, 0.0);
    }
  }

  //    Hashing members and functions
  // ===================================

  std::vector<short unsigned int> hash_temp;
  std::vector<short unsigned int> hash_dens;
  real8 temp_space;
  real8 dens_space;

  void set_hash_tables();

  real8 f_temp_hash(real8 T) {
    return temp_space * log(T/temps[0]);
  }
  real8 f_dens_hash(real8 rho) {
    return dens_space * log(rho/densities[0]);
  }

  int get_t_index(real8 temperature, OPAC_STATUS& t_err);
  int get_d_index(real8 density, OPAC_STATUS& d_err);

  //    Internal member functions
  // ===============================

  // File reading
  void read_file(std::string filename);
  void compress(std::vector<std::vector<std::vector<real8> > >& fine_r,
    std::vector<std::vector<std::vector<real8> > >& fine_p,
    std::vector<real8>& fine_group_bounds);

  // Retrieval auxillaries
  void interpolate(std::vector<real8>& vec,
    std::vector<std::vector<std::vector<real8> > >* opac,
    real8 temp, real8 dens, size_t t_index, size_t d_index, size_t g_min,
    size_t g_max, bool T_interp, bool rho_interp);

  real8 get_grey_opacity(real8 temperature, real8 m_density,
    OPAC_STATUS& error, OPAC_STATUS& d_err, opac_type type = OT_COUNT);

   real8 get_model_opacity(real8 temperature, real8 E_min = 0,
    real8 E_max = 1.0e48, opac_type type = OT_COUNT);

  void get_model_opacity(std::vector<real8>& vec, real8 temperature, size_t g_min,
    size_t g_max, opac_type type = OT_COUNT);

  real8 get_abc_opacity(real8 temperature) {
    return abc[0] + abc[1] * std::pow( temperature, abc[2] );
  }

public:

  //    Constructors and destructors
  // ===========================================================================

  // dummy constructor for BaseComponent when running neutronics
  Opacities(): opac_units(UT_cm2), microscopic(true), OT_default(OT_model_planck),
    mod_const(0.0), num_temps(0), num_dens(0), num_groups(0), temp_space(0.0),
    dens_space(0.0)
  {}

  // real constructor
  Opacities(std::string _id) : id(_id), opac_units(UT_cm2), microscopic(true),
    OT_default(OT_model_planck), mod_const(0.0), num_temps(0), num_dens(0),
    num_groups(0), temp_space(0.0), dens_space(0.0)
  {};

  Opacities(const Opacities& other, const std::string _id /*= "NULL"*/);

  Opacities& operator=(const Opacities& other);

  void define_type(stapl::typer&);

  ~Opacities() {};

  //    Initialization methods
  // ============================

  void set_default(opac_type _default) {
    OT_default = _default;
    opac_store[OT_default] = true;
  }
  void set_model(real8 _mod_const) {
    set_default(OT_model_planck);
    mod_const = _mod_const;
  }
  void set_abc(real8 a, real8 b, real8 c) {
    set_default(OT_abc);
    abc.clear(); abc.resize(3);
    abc[0] = a; abc[1] = b; abc[2] = c;
  }
  void initialize(std::string filename, opac_type _default);
  void set_group_structure(std::vector<real8> g_s);
  void store(opac_type t) {
    opac_store[t] = true;
  }

  //    Get and Set functions
  // ===========================================================================

  std::string get_id()                          { return id; }
  size_t get_num_groups()                       { return num_groups; }
  void get_group_bounds(std::vector<real8>& vec)     { vec = group_structure; }
  UNIT_TYPE get_units()                         { return opac_units; }
  bool is_microscopic()                         { return microscopic; }

  std::vector<real8> get_group_bounds(size_t group_number);
  std::vector<int> get_group_nums(real8 E_min = 0, real8 E_max = 1.0e48);

  // output funcitons
  void full_write(std::ostream& write);
  void hash_write(std::ostream& hout);

  //    Opacity retrieval
  // ===========================================================================

  void get_opacities(std::vector<real8>& vec, real8 temperature, real8 m_density,
    OPAC_STATUS& t_err, OPAC_STATUS& d_err, opac_type type = OT_COUNT)
  {
    get_opacities(vec, temperature, m_density, t_err, d_err, 0, num_groups,
      type);
  }
  void get_opacities(std::vector<real8>& vec, real8 temperature, real8 m_density,
    OPAC_STATUS& t_err, OPAC_STATUS& d_err, size_t g_min, size_t g_max,
    opac_type type = OT_COUNT);

};

// =============================================================================
//     multigroup retrieval
// =============================================================================

inline void Opacities::get_opacities(std::vector<real8> &vec, real8 temperature,
  real8 m_density, OPAC_STATUS& t_err, OPAC_STATUS& d_err, size_t g_min,
  size_t g_max, opac_type type)
{
  // check for valid input arguments
  // ===============================
  if( type == OT_COUNT )
    type = OT_default;
  if( !opac_store[type] )
  {
    std::stringstream str;
    str << "A request was made for non-existant opacity data.";
    throw CommonException(str, CET_INTERNAL_ERROR);
  }
  if( temperature < 0.0 )
  {
    std::stringstream str;
    str << "The supplied temperature for opacity evaluation must be "
        << "non-negative.";
    throw CommonException(str, CET_INTERNAL_ERROR);
  }
  if( m_density < 0.0 )
  {
    std::stringstream str;
    str << "The supplied density for opacity evaluation must be "
        << "non-negative.";
    throw CommonException(str, CET_INTERNAL_ERROR);
  }
  if( g_max > num_groups )
    g_max = num_groups;

  // call appropriate function
  // =========================
  std::vector<real8> grey_sig;
  std::vector<std::vector<std::vector<real8> > >* opac = NULL;
  switch( type )
  {
    case OT_model_ross:
    case OT_model_planck:
      get_model_opacity(vec, temperature, g_min, g_max, type);
      return; break;
    case OT_grey_ross:
    case OT_grey_planck:
      grey_sig.push_back(get_grey_opacity(temperature, m_density, t_err, d_err,
        type));
      vec = grey_sig;  return;
      break;
    case OT_mg_ross:
      opac = &ross;
      break;
    case OT_mg_planck:
      opac = &planck;
      break;
    case OT_abc:
      grey_sig.push_back( get_abc_opacity(temperature) );
      vec = grey_sig;  return;
    default:
      std::stringstream str;
      str << "Invalid opacity type.";
      throw CommonException(str, CET_INTERNAL_ERROR);
      break;
  }

  // carry on with multigroup retrieval
  // ==================================
  bool interp_over_T = true, interp_over_rho = true;
  int t_index, d_index;

  t_index = get_t_index(temperature, t_err);
  if( t_index < 0 )
  {
    interp_over_T = false;
    t_index = -t_index - 1;
  }

  d_index = get_d_index(m_density, d_err);
  if( d_index < 0 )
  {
    interp_over_rho = false;
    d_index = -d_index - 1;
  }

  // interpolate and return
  // ======================
  interpolate(vec, opac, temperature, m_density, t_index,
    d_index, g_min, g_max, interp_over_T, interp_over_rho);

}

// interpolation routine for multigroup
// ====================================

inline void Opacities::interpolate(std::vector<real8> &vec,
  std::vector<std::vector<std::vector<real8> > >* opac, real8 temp, real8 dens,
  size_t t_index, size_t d_index, size_t g_min, size_t g_max, bool T_interp,
  bool rho_interp)
{
  vec.clear();
  // this is a simple bilinear interpolation over log(rho) and log(T)
  // we will probably implement a better scheme at some point
  if( T_interp && rho_interp)
  {
    for (size_t g=g_min; g<g_max; g++)
    {
      real8 d2_d1 = log(densities[d_index+1] / densities[d_index]);
      real8 d_d1 = log(dens / densities[d_index]);
      real8 f_d = d_d1 / d2_d1;

      // interpolate over density at lower temperature
      real8 a = (1 - f_d) * (*opac)[t_index][d_index][g]
                     + f_d * (*opac)[t_index][d_index+1][g];

      // interpolate over density at higher temperature
      real8 b = (1 - f_d) * (*opac)[t_index+1][d_index][g]
                     + f_d * (*opac)[t_index+1][d_index+1][g];

      // interpolate over temperature
      real8 t2_t1 = log(temps[t_index+1] / temps[t_index]);
      real8 t_t1 = log(temp / temps[t_index]);
      real8 f_t = t_t1 / t2_t1;

      vec.push_back( (1 - f_t) * a + f_t * b );
    }
  }
  else if( T_interp && ! rho_interp )
  {
    for (size_t g=g_min; g<g_max; g++)
    {
      real8 a = (*opac)[t_index][d_index][g];
      real8 b = (*opac)[t_index+1][d_index][g];

      // interpolate over temperature
      real8 t2_t1 = log(temps[t_index+1] / temps[t_index]);
      real8 t_t1 = log(temp / temps[t_index]);
      real8 f_t = t_t1 / t2_t1;

      vec.push_back( (1 - f_t) * a + f_t * b );
    }
  }
  else if( ! T_interp && rho_interp )
  {
    for (size_t g=g_min; g<g_max; g++)
    {
      real8 d2_d1 = log(densities[d_index+1] / densities[d_index]);
      real8 d_d1 = log(dens / densities[d_index]);
      real8 f_d = d_d1 / d2_d1;

      real8 a = (1 - f_d) * (*opac)[t_index][d_index][g]
                     + f_d * (*opac)[t_index][d_index+1][g];

      vec.push_back( a );
    }
  }
  else if( ! T_interp && ! rho_interp )
  {
    for (size_t g=g_min; g<g_max; g++)
    {
      vec.push_back( (*opac)[t_index][d_index][g] );
    }
  }

}


// =============================================================================
//    grey retrieval
// =============================================================================

inline real8 Opacities::get_grey_opacity(real8 temperature, real8 m_density,
  OPAC_STATUS& t_err, OPAC_STATUS& d_err, opac_type type)
{
  bool interp_over_T = true, interp_over_rho = true;
  int t_index, d_index;

  t_index = get_t_index(temperature, t_err);
  if( t_index < 0 )
  {
    interp_over_T = false;
    t_index = -t_index - 1;
  }

  d_index = get_d_index(m_density, d_err);
  if( d_index < 0 )
  {
    interp_over_rho = false;
    d_index = -d_index - 1;
  }
  if( interp_over_T )
    assert( temps[t_index] <= temperature && temperature <= temps[t_index+1] );
  if( interp_over_rho)
    assert( densities[d_index]<=m_density&&m_density<=densities[d_index+1] );

  // interpolate
  real8 value;
  std::vector<std::vector<real8> >* opac;
  switch( type )
  {
    case OT_grey_planck:    opac = &grey_planck;    break;
    case OT_grey_ross:      opac = &grey_ross;      break;
    default:
      std::stringstream str;  WHERE(str);
      str << "Invalid opacity type requested.";
      throw CommonException( str, CET_INTERNAL_ERROR );
      break;
  }

  if( interp_over_T )
  {
    if( interp_over_rho )
    {
      real8 d2_d1 = log(densities[d_index+1] / densities[d_index]);
      real8 d_d1 = log(m_density / densities[d_index]);
      real8 f_d = d_d1 / d2_d1;
      real8 a = (1 - f_d) * (*opac)[t_index][d_index]
                     + f_d * (*opac)[t_index][d_index+1];
      real8 b = (1 - f_d) * (*opac)[t_index+1][d_index]
                     + f_d * (*opac)[t_index+1][d_index+1];
      real8 t2_t1 = log(temps[t_index+1] / temps[t_index]);
      real8 t_t1 = log(temperature / temps[t_index]);
      real8 f_t = t_t1 / t2_t1;
      value = (1 - f_t) * a + f_t * b;
    }
    else
    {
      real8 a = (*opac)[t_index][d_index];
      real8 b = (*opac)[t_index+1][d_index];
      real8 t2_t1 = log(temps[t_index+1] / temps[t_index]);
      real8 t_t1 = log(temperature / temps[t_index]);
      real8 f_t = t_t1 / t2_t1;
      value = (1 - f_t) * a + f_t * b;
    }
  }
  else // don't interpolate over T
  {
    if( interp_over_rho )
    {
      real8 d2_d1 = log(densities[d_index+1] / densities[d_index]);
      real8 d_d1 = log(m_density / densities[d_index]);
      real8 f_d = d_d1 / d2_d1;
      value = (1 - f_d) * (*opac)[t_index][d_index]
                  + f_d * (*opac)[t_index][d_index+1];
    }
    else
    {
      value = (*opac)[t_index][d_index];
    }
  }
  return value;
}

// =============================================================================
//    model retrieval
// =============================================================================

/*!
   The model opacity  given in "Asymptotic Analysis of a Computational Method
   for Time- and Frequency- Dependent Radiative Transfer", Adams and Nowak,
   Journal of Computational Physics 146, 366-403 (1998) is:

    \f$ \sigma(E,T) = \sigma_0 frac(1-exp(-E/kT),E^3) \f$.

    To find multigroup (Planck weighted) model opacities, we need

    \f$ \int_E1^E2 {\sigma(E,T)*B(E,T)}dE / \int_E1^E2 {B(E,T)}dE \f$

    which yields

    \f$ 2*\sigma_0/(h^3 c^2)*(exp(-E1/kT)-exp(-E2/kT) /
    \int_E1^E2 {B(E,T)}dE \f$.
*/

inline real8 Opacities::get_model_opacity(real8 temperature, real8 E_min,
  real8 E_max, opac_type type)
{
  // Note that our units for B are NOT per steradian

  Logger msg;
  real8 mod;

  if( type == OT_COUNT )
    type = OT_default;

  assert( temperature >= 0.0 );
  assert( E_min >= 0.0 );
  assert( E_max >= E_min );

  if( type == OT_model_ross )
    msg.log(WARN) << "Rosseland-weighted model opacity has not "
                  << "been implemented; Planck values are being used." << std::endl;

  real8 K  = BOLTZMANN_CONSTANT; // eV K^-1
  real8 H  = PLANCK_CONSTANT;    // eV s
  real8 C  = SPEED_OF_LIGHT;     // cm s^-1
  real8 Pi = PI;

  Planck planck1;
  real8 Bg = planck1.integrate_B(temperature, E_min, E_max);

  if (Bg == 0)
  {
    mod = OPACITY_MAX;
    // I need to do this integral properly; but 10^6 is better than:
    // real8 E = .5*(E_min + E_max);
    // mod = mod_const*(1-exp(-E/(k*temperature))*std::pow(E,-3));
  }
  else
  {
    real8 opacity = 8 * Pi * mod_const * std::pow(H,-3.0) * std::pow(C,-2.0) *
     (exp(-E_min/(K*temperature)) - exp(-E_max/(K*temperature)));
    mod = opacity / Bg ;
  }
  return mod;
}

inline void Opacities::get_model_opacity(std::vector<real8> &vec,
  real8 temperature, size_t g_min, size_t g_max, opac_type type)
{
  // Note that our units for B are NOT per steradian

  Logger msg;
  vec.clear();

  if( type == OT_model_ross )
    msg.log(WARN) << "Rosseland-weighted model opacity has not "
                  << "been implemented; Planck values are being used." << std::endl;

  const real8 H  = PLANCK_CONSTANT;      // eV s
  const real8 K  = BOLTZMANN_CONSTANT;   // eV K^-1
  const real8 C  = SPEED_OF_LIGHT;       // cm s^-1
  const real8 Pi = PI;
  const real8 A  = RADIATION_CONSTANT_A; // eV cm^-3 K^-4
  const real8 F  = 1; // some sort of units constant

  if( num_groups == 1 )
  {
    real8 grey_mod = 8 * Pi * mod_const * K * F / (std::pow(H,3.0) *
      std::pow(C,3.0) * std::pow(temperature,3.0) * A);
    vec.push_back(grey_mod);
  }
  else
  {
    Planck planck1;

    for( size_t i=g_min; i<g_max; i++)
    {
      real8 Bg = planck1.integrate_B(temperature, group_structure[i+1],
        group_structure[i]);

      if (Bg == 0)
      {
        real8 E = .5*(group_structure[i] + group_structure[i+1]);
        vec.push_back(mod_const*(1-exp(-E/(K*temperature)))*std::pow(E,-3.0));
      }
      else
      {
        real8 opacity = 8 * Pi * mod_const * std::pow(H,-3.0) *
          std::pow(C,-2.0) * K * temperature *
          (exp(-group_structure[i+1]/(K*temperature)) -
          - exp(-group_structure[i]/(K*temperature)));
        vec.push_back(opacity/Bg);
      }
    }
  }
}

// =============================================================================
//    index retrieval and interpolation decisions
// =============================================================================

inline int Opacities::get_t_index(real8 temperature, OPAC_STATUS& t_err)
{
  // if interpolation is unnecessary (or impossible), we return a negative value
  // as indication to the calling function.  since we need this information even
  // (or especially) for a 0 index, we add -1 to any negative index.  mpa.
  int t_index;
  t_err = OS_VALID;

  if( num_temps == 1 )  // interpolation impossible
  {
    // indicates 0 index and no interpolation
    t_index = -1;
  }
  else if( temperature < temps[0] )  // we don't extrapolate
  {
    t_err = OS_LOW_T;
    t_index = -1;
  }
  else if( temperature > temps.back() )
  {
    t_err = OS_HIGH_T;
    // indicates index of temps.size() - 1 and no interpolation
    t_index = -num_temps;
  }
  else
  {
    size_t hash_index = size_t(floor(f_temp_hash(temperature)));
    if( hash_index >= hash_temp.size() ) {
      std::stringstream str;  WHERE(str);
      str << "Error with temperature hash in Opacities class:" << std::endl
        << "Requested temperature was " << temperature << "K.  Maximum "
        << "temperature is " << temps.back() << "K.";
      throw CommonException( str, CET_INTERNAL_ERROR );
    }
    t_index = hash_temp[hash_index];
    // this should be an if; i just have to demonstrate that that's robust
    while (temps[t_index] > temperature)
      t_index--;
    if( temperature == temps[t_index] )
      t_index = -t_index - 1;
  }

  return t_index;
}

inline int Opacities::get_d_index(real8 density, OPAC_STATUS& d_err)
{
  // if interpolation is unnecessary (or impossible), we return a negative value
  // as indication to the calling function.  since we need this information even
  // (or especially) for a 0 index, we add -1 to any negative index.  mpa.
  int d_index;
  d_err = OS_VALID;

  if( num_dens == 1 )  // interpolation impossible
  {
    // indicates 0 index and no interpolation
    d_index = -1;
  }
  else if( density < densities[0] )  // we don't extrapolate
  {
    d_err = OS_LOW_RHO;
    d_index = -1;
  }
  else if( densities.back() < density )
  {
    d_err = OS_HIGH_RHO;
    // indicates index of densities.size() - 1 and no interpolation
    d_index = -num_dens;
  }
  else
  {
    d_index = hash_dens[size_t(floor(f_dens_hash(density)))];
    // this should be an if; i just have to demonstrate that that's robust
    while( density < densities[d_index] )
      d_index--;
    if( density == densities[d_index] )
      d_index = -d_index - 1;
  }

  return d_index;
}

// =============================================================================
//    inline methods
// =============================================================================

inline std::vector<real8> Opacities::get_group_bounds(size_t group_number)
{
  Logger msg;
  std::vector<real8> bounds (2, 0.0);
  if( num_groups == 1 )
  {
    msg.log(WARN) << "Group bounds requested in non-multigroup "
                  << "opacities object " << id << "." << std::endl;
    bounds[0] = 0;
    bounds[1] = 1.0e48;
  }
  else
  {
    assert( group_number < num_groups );
    bounds[0] = group_structure[group_number];
    bounds[1] = group_structure[group_number+1];
  }
  return bounds;
}

inline std::vector<int> Opacities::get_group_nums(real8 E_min, real8 E_max)
{
  std::vector<int> group_nums;
  if( num_groups > 1 )
  {
    for (size_t i=0; i<num_groups; i++)
    {
      if ( E_min < group_structure[i] && group_structure[i+1] < E_max)
        group_nums.push_back(i);
    }
  }
  else
  {
    group_nums.push_back(0);
  }
  return group_nums;
}

#endif
