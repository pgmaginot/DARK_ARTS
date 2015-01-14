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
// #include "Edits.h"
// #include "basic_macros.h"

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

  OT_COUNT                    // a way to size a bitset
  // *** DO NOT ADD TYPES AFTER OT_COUNT !!!
};

enum UNIT_TYPE
{
  // Length
  UT_cm,          // centimeter
  UT_m,           // meter

  // Mass
  UT_g,           // gram
  UT_kg,          // kilogram

  // Density
  UT_Npbcm,       // number / barn-centimeter
  UT_gpcm3,       // g / cm^3         ( gram per centimeter-cubed )
  UT_kgpm3,       // kg / m^3         ( kilogram per meter-cubed )
  UT_Npcm3,       // number / cm^3    ( number of parts per centimeter-cubed )
  UT_Npm3,        // number / m^3     ( number of parts per meter-cubed )

  // Temperature
  UT_K,           // Kelvin
  UT_eV,          // electron-Volt

  // Specific Heat
  UT_eVpgK,       // eV / (g-K)

  // Scalar flux
  UT_npcm2s,      // number / (cm^2-s)

  // Arbsorption/Emission Rate Density
  UT_eVpcm3s,     // eV / (cm^3-s)

  // Thermal radiation opacity
  UT_cm2pg,       // cm^2 / g         ( centimeters-squared per gram )
  UT_m2pkg,       // m^2 / kg         ( meters-squared per kilogram )
  UT_icm,         // cm^-1            ( inverse centimeter [macroscopic] )
  UT_im,          // m^-1             ( inverse meter [macroscopic] )

  // Neutron/Gamma cross sections
  UT_barns,       // barn = 1.0e-24 cm^2 = 1 angstrom^2
  UT_cm2          // cm^2
};

// =============================================================================
//    Opacities class
// =============================================================================

class Opacities
{
protected:
  Planck planck_func;
  
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
  /// for Adams-Novak model opacity
  double mod_const;
  /// single energy group data
  std::vector<std::vector<double> > grey_planck;
  std::vector<std::vector<double> > grey_ross;
  /// multigroup data
  std::vector<std::vector<std::vector<double> > > planck;
  std::vector<std::vector<std::vector<double> > > ross;
  std::vector<double> abc;

  //    Independent variable data members
  // =======================================
  std::vector<double> temps;
  std::vector<double> densities;
  std::vector<double> group_structure;
  size_t num_temps;
  size_t num_dens;
  size_t num_groups;

  /// inline vector of vector resizing functions
  void resize(std::vector<std::vector<double> >& vec, size_t a, size_t b)
  {
    vec.clear();  vec.resize(a);
    for( size_t i=0; i!=a; ++i )
      vec[i].resize(b, 0.0);
  }
  void resize(std::vector<std::vector<std::vector<double> > >& v, size_t a,
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
  double temp_space;
  double dens_space;

  void set_hash_tables();

  double f_temp_hash(double T) {
    return temp_space * log(T/temps[0]);
  }
  double f_dens_hash(double rho) {
    return dens_space * log(rho/densities[0]);
  }

  int get_t_index(double temperature, OPAC_STATUS& t_err);
  int get_d_index(double density, OPAC_STATUS& d_err);

  //    Internal member functions
  // ===============================

  // File reading
  void read_file(std::string filename);
  void compress(std::vector<std::vector<std::vector<double> > >& fine_r,
    std::vector<std::vector<std::vector<double> > >& fine_p,
    std::vector<double>& fine_group_bounds);

  // Retrieval auxillaries
  void interpolate(std::vector<double>& vec,
    std::vector<std::vector<std::vector<double> > >* opac,
    double temp, double dens, size_t t_index, size_t d_index, size_t g_min,
    size_t g_max, bool T_interp, bool rho_interp);

  double get_grey_opacity(double temperature, double m_density,
    OPAC_STATUS& error, OPAC_STATUS& d_err, opac_type type = OT_COUNT);

public:

  //    Constructors and destructors
  // ===========================================================================

  // dummy constructor for BaseComponent when running neutronics
  Opacities(const Input_Reader& input_reader): 
    planck_func(1.0E-12, input_reader, 2.0) , 
    opac_units(UT_cm2), microscopic(true), OT_default(OT_mg_ross),
    mod_const(0.0), num_temps(0), num_dens(0), num_groups(0), temp_space(0.0),
    dens_space(0.0)
  {}

  // real constructor
  Opacities(const Input_Reader& input_reader, std::string _id) : 
    planck_func(1.0E-12, input_reader, 2.0) ,
    id(_id), opac_units(UT_cm2), microscopic(true),
    OT_default(OT_mg_ross), mod_const(0.0), num_temps(0), num_dens(0),
    num_groups(0), temp_space(0.0), dens_space(0.0)
  {};

  ~Opacities() {};

  //    Initialization methods
  // ============================

  void set_default(opac_type _default) {
    OT_default = _default;
    opac_store[OT_default] = true;
  }

  void initialize(std::string filename, opac_type _default);
  void set_group_structure(std::vector<double> g_s);
  void store(opac_type t) {
    opac_store[t] = true;
  }

  //    Get and Set functions
  // ===========================================================================

  std::string get_id()                          { return id; }
  size_t get_num_groups()                       { return num_groups; }
  void get_group_bounds(std::vector<double>& vec)     { vec = group_structure; }
  UNIT_TYPE get_units()                         { return opac_units; }
  bool is_microscopic()                         { return microscopic; }

  std::vector<double> get_group_bounds(size_t group_number);
  std::vector<int> get_group_nums(double E_min = 0, double E_max = 1.0e48);

  // output funcitons
  void full_write(std::ostream& write);
  void hash_write(std::ostream& hout);

  //    Opacity retrieval
  // ===========================================================================

  void get_opacities(std::vector<double>& vec, double temperature, double m_density,
    OPAC_STATUS& t_err, OPAC_STATUS& d_err, opac_type type = OT_COUNT)
  {
    get_opacities(vec, temperature, m_density, t_err, d_err, 0, num_groups,
      type);
  }
  void get_opacities(std::vector<double>& vec, double temperature, double m_density,
    OPAC_STATUS& t_err, OPAC_STATUS& d_err, size_t g_min, size_t g_max,
    opac_type type = OT_COUNT);

};

// =============================================================================
//     multigroup retrieval
// =============================================================================
/// return a vector of evaluated opacities in vec evalauted at temperature and m_density
inline void Opacities::get_opacities(std::vector<double> &vec, double temperature,
  double m_density, OPAC_STATUS& t_err, OPAC_STATUS& d_err, size_t g_min,
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
    throw Dark_Arts_Exception(SUPPORT_OBJECT, str);
  }
  if( temperature < 0.0 )
  {
    std::stringstream str;
    str << "The supplied temperature for opacity evaluation must be "
        << "non-negative.";
    throw Dark_Arts_Exception(SUPPORT_OBJECT, str);
  }
  if( m_density < 0.0 )
  {
    std::stringstream str;
    str << "The supplied density for opacity evaluation must be "
        << "non-negative.";
    throw Dark_Arts_Exception(SUPPORT_OBJECT, str);
  }
  if( g_max > num_groups )
    g_max = num_groups;

  // call appropriate function
  // =========================
  std::vector<double> grey_sig;
  std::vector<std::vector<std::vector<double> > >* opac = NULL;
  switch( type )
  {
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
    default:
      std::stringstream str;
      str << "Invalid opacity type.";
      throw Dark_Arts_Exception(SUPPORT_OBJECT, str);
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
/// vec holds the multigroup opacities evalauted at temperature temp
inline void Opacities::interpolate(std::vector<double> &vec,
  std::vector<std::vector<std::vector<double> > >* opac, double temp, double dens,
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
      double d2_d1 = log(densities[d_index+1] / densities[d_index]);
      double d_d1 = log(dens / densities[d_index]);
      double f_d = d_d1 / d2_d1;

      // interpolate over density at lower temperature
      double a = (1 - f_d) * (*opac)[t_index][d_index][g]
                     + f_d * (*opac)[t_index][d_index+1][g];

      // interpolate over density at higher temperature
      double b = (1 - f_d) * (*opac)[t_index+1][d_index][g]
                     + f_d * (*opac)[t_index+1][d_index+1][g];

      // interpolate over temperature
      double t2_t1 = log(temps[t_index+1] / temps[t_index]);
      double t_t1 = log(temp / temps[t_index]);
      double f_t = t_t1 / t2_t1;

      vec.push_back( (1 - f_t) * a + f_t * b );
    }
  }
  else if( T_interp && ! rho_interp )
  {
    for (size_t g=g_min; g<g_max; g++)
    {
      double a = (*opac)[t_index][d_index][g];
      double b = (*opac)[t_index+1][d_index][g];

      // interpolate over temperature
      double t2_t1 = log(temps[t_index+1] / temps[t_index]);
      double t_t1 = log(temp / temps[t_index]);
      double f_t = t_t1 / t2_t1;

      vec.push_back( (1 - f_t) * a + f_t * b );
    }
  }
  else if( ! T_interp && rho_interp )
  {
    for (size_t g=g_min; g<g_max; g++)
    {
      double d2_d1 = log(densities[d_index+1] / densities[d_index]);
      double d_d1 = log(dens / densities[d_index]);
      double f_d = d_d1 / d2_d1;

      double a = (1 - f_d) * (*opac)[t_index][d_index][g]
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

inline double Opacities::get_grey_opacity(double temperature, double m_density,
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
  double value;
  std::vector<std::vector<double> >* opac;
  switch( type )
  {
    case OT_grey_planck:    opac = &grey_planck;    break;
    case OT_grey_ross:      opac = &grey_ross;      break;
    default:
      std::stringstream str; 
      str << "Invalid opacity type requested.";
      throw Dark_Arts_Exception(SUPPORT_OBJECT , str);
      break;
  }

  if( interp_over_T )
  {
    if( interp_over_rho )
    {
      double d2_d1 = log(densities[d_index+1] / densities[d_index]);
      double d_d1 = log(m_density / densities[d_index]);
      double f_d = d_d1 / d2_d1;
      double a = (1 - f_d) * (*opac)[t_index][d_index]
                     + f_d * (*opac)[t_index][d_index+1];
      double b = (1 - f_d) * (*opac)[t_index+1][d_index]
                     + f_d * (*opac)[t_index+1][d_index+1];
      double t2_t1 = log(temps[t_index+1] / temps[t_index]);
      double t_t1 = log(temperature / temps[t_index]);
      double f_t = t_t1 / t2_t1;
      value = (1 - f_t) * a + f_t * b;
    }
    else
    {
      double a = (*opac)[t_index][d_index];
      double b = (*opac)[t_index+1][d_index];
      double t2_t1 = log(temps[t_index+1] / temps[t_index]);
      double t_t1 = log(temperature / temps[t_index]);
      double f_t = t_t1 / t2_t1;
      value = (1 - f_t) * a + f_t * b;
    }
  }
  else // don't interpolate over T
  {
    if( interp_over_rho )
    {
      double d2_d1 = log(densities[d_index+1] / densities[d_index]);
      double d_d1 = log(m_density / densities[d_index]);
      double f_d = d_d1 / d2_d1;
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
//    index retrieval and interpolation decisions
// =============================================================================

inline int Opacities::get_t_index(double temperature, OPAC_STATUS& t_err)
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
      std::stringstream str;  
      str << "Error with temperature hash in Opacities class:" << std::endl
        << "Requested temperature was " << temperature << "K.  Maximum "
        << "temperature is " << temps.back() << "K.";
       throw Dark_Arts_Exception(SUPPORT_OBJECT , str);
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

inline int Opacities::get_d_index(double density, OPAC_STATUS& d_err)
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

inline std::vector<double> Opacities::get_group_bounds(size_t group_number)
{
  std::vector<double> bounds (2, 0.0);
  if( num_groups == 1 )
  {
    std::stringstream err;
    err << "Group bounds requested in non-multigroup "
                  << "opacities object " << id << "." << std::endl;
    throw Dark_Arts_Exception(SUPPORT_OBJECT, err);
    
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

inline std::vector<int> Opacities::get_group_nums(double E_min, double E_max)
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
