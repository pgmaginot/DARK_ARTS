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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include "Opacities.h"

using std::string;
using std::ifstream;
using std::endl;
using std::setprecision;
using std::setw;
using std::ios;
using std::streambuf;
using std::pow;
using std::stringstream;
using std::ostream;
using std::vector;

//    Initialize from data file
// ===============================
void Opacities::initialize(string filename, opac_type _def)
{
  set_default(_def);
  read_file(filename);
  if( num_temps > 1 || num_dens > 1 )
    set_hash_tables();
}

void Opacities::set_group_structure(vector<real8> g_s)
{
  // Invalid group structure
  // =======================
  if( g_s.size() < 2 )
  {
    stringstream str; WHERE(str);
    str << "Only one group boundary passed to Opacities object.";
    throw CommonException( str, CET_INTERNAL_ERROR );
  }

  // Grey problem
  // ============
  else if( g_s.size() == 2 )
  {
    // For non-vector types, just set object members
    if( OT_default == OT_model_ross || OT_default == OT_model_planck
        || OT_default == OT_abc )
    {
      group_structure = g_s;
      num_groups = group_structure.size() - 1;
      return;
    }
    // Multigroup data, compress to grey
    if( group_structure.size() > 2 )
    {
      vector<vector<vector<real8> > > fine_r = ross;
      vector<vector<vector<real8> > > fine_p = planck;
      vector<real8> fine_groups = group_structure;
      group_structure = g_s;
      num_groups = group_structure.size() - 1;
      compress(fine_r, fine_p, fine_groups);
      // Populate grey structures
      if( OT_default == OT_grey_ross )
      {
        for( size_t t=0; t!=num_temps; ++t ) {
          for( size_t d=0; d!=num_dens; ++d )
            grey_ross[t][d] = ross[t][d][0];
        }
      }
      if( OT_default == OT_grey_planck )
      {
        for( size_t t=0; t!=num_temps; ++t ) {
          for( size_t d=0; d!=num_dens; ++d )
            grey_planck[t][d] = planck[t][d][0];
        }
      }
    }
    // Set object members
    group_structure = g_s;
    num_groups = group_structure.size() - 1;
  }

  // Multigroup
  // ==========
  else if( g_s.size() > 2 )
  {
    // Incompatible opacity types
    if( OT_default == OT_abc )
    {
      stringstream str;
      str << "ABC opacities are not yet compatible with multigroup.";
      throw CommonException( str, CET_DATA_ERROR );
    }
    else if( OT_default == OT_grey_ross || OT_default == OT_grey_planck )
    {
      stringstream str; WHERE(str);
      str << "Grey opacities specified for a multigroup problem.";
      throw CommonException( str, CET_INPUT_ERROR );
    }

    // Model opacity types (no data to alter)
    else if( OT_default == OT_model_planck || OT_default == OT_model_ross )
    {
      group_structure = g_s;
      num_groups = group_structure.size() - 1;
    }

    // Multigroup data
    else
    {
      // Missing data
      if( ( OT_default == OT_mg_ross && ross.empty() )
        || ( OT_default == OT_mg_planck && planck.empty() ) )
      {
        stringstream str; WHERE(str);
        str << "Missing multigroup data.";
        throw CommonException( str, CET_DATA_ERROR );
      }
      // Compress
      if( group_structure.size() > g_s.size() )
      {
        vector<vector<vector<real8> > > fine_r = ross;
        vector<vector<vector<real8> > > fine_p = planck;
        vector<real8> fine_groups = group_structure;
        group_structure = g_s;
        num_groups = group_structure.size() - 1;
        compress(fine_r, fine_p, fine_groups);
      }
      else if( group_structure.size() == g_s.size() )
      {
        num_groups = group_structure.size() - 1;
        for( size_t g=0; g!=num_groups; ++g )
        {
          if( group_structure[g] != g_s[g] ) {
            stringstream str; WHERE(str);
            str << "User supplied group structure must match material "
              << "data.  Mismatch in group "<< g << ".";
            throw CommonException( str, CET_DATA_ERROR );
          }
        }
      }
      else
      {
        stringstream str;  WHERE(str);
        str << "User supplied group boundaries must be a subset of "
          << "the group boundaries in the material data.";
        throw CommonException( str, CET_DATA_ERROR );
      }
    } // multigroup data

  } // multigroup problem

  // Discard unused data
  if( ! opac_store[OT_grey_ross] )    grey_ross.clear();
  if( ! opac_store[OT_grey_planck] )  grey_planck.clear();
  if( ! opac_store[OT_mg_ross] )      ross.clear();
  if( ! opac_store[OT_mg_planck] )    planck.clear();

}

// === Copy Constructor ===

Opacities::Opacities(const Opacities& other, const string _id)
{
  if (this != &other)
  {
    id = _id;
    opac_store = other.opac_store;
    OT_default = other.OT_default;
    num_groups = other.num_groups;
    group_structure = other.group_structure;
    num_temps = other.num_temps;
    num_dens = other.num_dens;
    mod_const = other.mod_const;
    abc = other.abc;
    microscopic = other.microscopic;
    opac_units = other.opac_units;
    if( opac_store[OT_grey_ross] )
      grey_ross = other.grey_ross;
    if( opac_store[OT_grey_planck] )
      grey_planck = other.grey_planck;
    if( opac_store[OT_mg_ross] )
      ross = other.ross;
    if( opac_store[OT_mg_planck] )
      planck = other.planck;
    if( num_temps > 1 )
    {
      temps = other.temps;
      hash_temp = other.hash_temp;
      temp_space = other.temp_space;
    }
    if( num_dens > 1 )
    {
      densities = other.densities;
      hash_dens = other.hash_dens;
      dens_space = other.dens_space;
    }
  }
}

// === Assignment Operator ===

Opacities& Opacities::operator=(const Opacities& other)
{
  if (this != &other)
  {
    if( id.empty() )
      id = other.id;
    opac_store = other.opac_store;
    OT_default = other.OT_default;
    num_groups = other.num_groups;
    group_structure = other.group_structure;
    num_temps = other.num_temps;
    num_dens = other.num_dens;
    mod_const = other.mod_const;
    abc = other.abc;
    microscopic = other.microscopic;
    opac_units = other.opac_units;
    if( opac_store[OT_grey_ross] )
      grey_ross = other.grey_ross;
    if( opac_store[OT_grey_planck] )
      grey_planck = other.grey_planck;
    if( opac_store[OT_mg_ross] )
      ross = other.ross;
    if( opac_store[OT_mg_planck] )
      planck = other.planck;
    if( num_temps > 1 )
    {
      temps = other.temps;
      hash_temp = other.hash_temp;
      temp_space = other.temp_space;
    }
    if( num_dens > 1 )
    {
      densities = other.densities;
      hash_dens = other.hash_dens;
      dens_space = other.dens_space;
    }
  }
  return *this;
}


// =============================================================================
//                           MEMBER FUNCTIONS
// =============================================================================

void Opacities::define_type(stapl::typer& t)
{
  t.member(id);
  t.member(opac_units);
  t.member(microscopic);
  t.member(opac_store);
  t.member(OT_default);
  t.member(mod_const);
  t.member(grey_planck);
  t.member(grey_ross);
  t.member(planck);
  t.member(ross);
  t.member(abc);
  t.member(temps);
  t.member(densities);
  t.member(group_structure);
  t.member(num_temps);
  t.member(num_dens);
  t.member(num_groups);
  t.member(hash_temp);
  t.member(hash_dens);
  t.member(temp_space);
  t.member(dens_space);
}

void Opacities::set_hash_tables()
{
  if( num_temps > 1 )
  {
    // initialize temperature hashing function
    //   setting the spacing to minimize sharing buckets
    size_t min = 0;
    for( size_t i=0; i<num_temps-1; i++ )
    {
      if ( temps[i+1] / temps[i] < temps[min+1] / temps[min] )
        min = i;
    }
    temp_space = 1 / log( temps[min+1] / temps[min]);

    // set temperature hashing table
    hash_temp.resize(0);
    real8 j=0;
    for( size_t i=0; i<num_temps; i++ )
    {
      real8 f_value = f_temp_hash(temps[i]);
      while (j <= f_value)
      {
        hash_temp.push_back(i);
        j++;
      }
    }
    // to ensure that hash_temp covers the full space of f_temp_hash
    hash_temp.push_back(num_temps-1);
  }   // if( num_temps > 1 )

  if( num_dens > 1 )
  {
    // initialize density hashing function
    size_t min = 0;
    for( size_t i=0; i<num_dens-1; i++ )
    {
      if ( densities[i+1] / densities[i] < densities[min+1] / densities[min] )
        min = i;
    }
    dens_space = 1 / log( densities[min+1] / densities[min]);

    // set density hashing table
    hash_dens.resize(0);
    real8 j=0;
    for( size_t i=0; i<num_dens; i++ )
    {
      real8 f_value = f_dens_hash(densities[i]);
      while( j <= f_value )
      {
        hash_dens.push_back(i);
        j++;
      }
    }
    // to ensure that we encapsulate full space of f_dens_hash in hash_dens
    hash_dens.push_back(num_dens-1);
  }   // if( num_dens > 1 )

} // Opacities::set_hash_tables()


// =============================================================================
//    output functions
// =============================================================================

void Opacities::hash_write(ostream& hout)
{
  if( num_temps == 1 )
  {
    hout << "There is no temperature hash because there is only one "
          << "temperature\n";
  }
  else
  {
    hout<< "f_temp_hash(T) = ln(T/temps[0]) * 1 / (ln( [minimum ratio] ) )"
          << "\n               = ln(T/" << temps[0] << ") * " << temp_space
          << "\n\nHash table\n==========\n";
    for( size_t i=0; i!=hash_temp.size(); ++i )
    {
      hout<< "hash_temp[" << setw(3) << i << "] = " << setw(3) << hash_temp[i];
      hout << "     temps[" <<setw(3) << hash_temp[i] << "] = "
            << temps[hash_temp[i]] << " K" << endl;
    }
    hout << endl << endl;
  }

  if( num_dens == 1 )
  {
    hout << "There is no density hash because there is only one density\n";
  }
  else
  {
    hout << "f_dens_hash(rho) = ln(rho/densities[0]) * 1 / "
          << "(ln( [minimum ratio] ) )\n                 = ln(rho/"
          << densities[0] << ") * " << dens_space << endl;
    hout << "\nHash table\n==========\n";
    for( size_t i=0; i!=hash_dens.size(); ++i )
    {
      hout<< "hash_dens[" << setw(3) << i << "] = " << setw(3) << hash_dens[i];
      hout << "     densities[" <<setw(3) << hash_dens[i] << "] = "
            << densities[hash_dens[i]] << " g/cc" << endl;
    }
    hout << endl << endl;
  }

}

void Opacities::full_write(ostream& write)
{
  write << "Opacities object " << id << ":" << endl;
  write << "==============================" << endl;

  // Model opacities
  if( OT_default == OT_model_ross || OT_default == OT_model_planck )
  {
    write << "Model opacity constant sigma_0 = " << mod_const
      << " ( eV^3 cm^2 ) / g" << endl;
    return;
  }

  // ABC opacities
  else if( OT_default == OT_abc )
  {
    write << "A = " << abc[0] << " cm^2 / g" << endl;
    write << "B = " << abc[1] << " cm^2 / ( g K^[" << abc[2] << "] )" << endl;
    write << "C = " << abc[2] << " [dimensionless]" << endl;
    return;
  }

  // Planck and Rosseland opacities
  size_t t_index, d_index, g_index;
  write << "Temperatures -\n";
  for (t_index=0; t_index<temps.size(); t_index++)
    write << setw(12) << temps[t_index] << " K = "
          << setw(12) << temps[t_index]*8.617343e-5 << " eV\n";
  write << "\nDensities (g/cc) -\n";
  for (d_index=0; d_index<densities.size(); d_index++)
    write << densities[d_index] << endl;
  if( opac_store[OT_mg_ross] || opac_store[OT_mg_planck] )
  {
    write << "\nCoarse group bounds (eV) -\n";
    for (g_index=0; g_index<=num_groups; g_index++)
      write << group_structure[g_index] << endl;
  }


  if( opac_store[OT_grey_ross] )
  {
    write << "\nRosseland Grey Opacities (cm^2/g) -\n"
          << "===================================\n";
    for (t_index=0; t_index<grey_ross.size(); t_index++)
    {
      write << "T = " << temps[t_index] << "K = " << 8.617343e-5*temps[t_index]
            << "eV\n";
      for (d_index=0; d_index<grey_ross[t_index].size(); d_index++)
      {
        write << "rho = " << setw(8) << densities[d_index] << " g/cc - "
              << setw(8) << grey_ross[t_index][d_index] << endl;
      }
      write << endl;
    }
    write << endl << endl;
  }

  if( opac_store[OT_grey_planck] )
  {
    write << "\nPlanck Grey Opacities (cm^2/g) -\n"
          << "===================================\n";
    for (t_index=0; t_index<grey_planck.size(); t_index++)
    {
      write << "T = " << temps[t_index] << "K = " << 8.617343e-5*temps[t_index]
            << "eV\n";
      for (d_index=0; d_index<grey_planck[t_index].size(); d_index++)
      {
        write << "rho = " << setw(8) << densities[d_index] << " g/cc - "
              << setw(8) << grey_planck[t_index][d_index] << endl;
      }
      write << endl;
    }
    write << endl << endl;
  }

  if( opac_store[OT_mg_ross] )
  {
    write << "\nRosseland Multigroup Opacities (cm^2 / g) -\n"
         << "===================================\n";
    for (t_index=0; t_index<ross.size(); t_index++)
    {
      for (d_index=0; d_index<ross[t_index].size(); d_index++)
      {
        write << "\nT = "<<temps[t_index] << "K = "
              << 8.617343e-5*temps[t_index]
              << "eV, rho = " << densities[d_index] << " g/cc\n";
        for (g_index=0; g_index<num_groups; g_index++)
          write << ross[t_index][d_index][g_index] << endl;
      }
    }
    write << endl << endl;
  }

  if( opac_store[OT_mg_planck] )
  {
    write << "\nPlanck Multigroup Opacities (cm^2 / g) -\n"
         << "===================================\n";
    for (t_index=0; t_index<planck.size(); t_index++)
    {
      for (d_index=0; d_index<planck[t_index].size(); d_index++)
      {
        write << "\nT = "<<temps[t_index] << "K = "
              << 8.617343e-5*temps[t_index]
              << "eV, rho = " << densities[d_index] << " g/cc\n";
        for (g_index=0; g_index<num_groups; g_index++)
          write << planck[t_index][d_index][g_index] << endl;
      }
    }
    write << endl << endl;
  }

} // Opacities::full_write()


// ============================================================================
//    File reading
// ============================================================================

void Opacities::read_file(string filename)
{
  ifstream fin( filename.c_str() );
  if( !fin.is_open() ) {
    stringstream str; WHERE(str);
    str << "Unable to open data file " << filename << ".";
    throw CommonException( str, CET_INPUT_ERROR );
  }

  // Factor for conversion to the units being used by PDT
  //   microscopic - cm^2/g
  //   macroscopic - cm^-1
  real8 conversion_factor;

  char line[250];
  string word;

  fin.getline(line, 250);
  fin.getline(line, 250);

  // "This file is a " multigroup / single-group
  fin >> word >> word >> word >> word;
  string multi;
  fin >> multi;

  // opacity / neutron / gamma / coupled neutron-gamma
  string type;
  fin >> type;
  fin.getline(line, 250);

  if( type != "opacity" ) {
    stringstream str; WHERE(str);
    str << "Opacities object trying to read non-opacity file "
      << filename << ".";
    throw CommonException( str, CET_DATA_ERROR );
  }

  fin >> num_temps >> word >> num_dens >> word >> word >> num_groups;
  fin.getline(line, 250);

  temps.resize( num_temps );
  densities.resize( num_dens );
  group_structure.resize( num_groups + 1 );

  size_t n_processes;
  fin >> n_processes;
  fin.getline(line, 250);
  fin.getline(line, 250);
  fin.getline(line, 250);

  // "All cross sections are in units of xxx"
  string micro, units;
  fin >> micro;
  fin >> word >> word >> word >> word >> word >> word;
  fin >> units;
  /* ***FIXME***
   * ***UNITS***
   * Right now I am implementing two contradictory responses to the unit issue.
   * The conversion_factor approach is quick, easy, and reasonably sensible.
   * The unit tracking is what I prefer, and I'm leaving it because it's a good
   * framework.  However, while the conversion_factor is being used, the label
   * of units is wrong.  We simply need to pick one and stick with it.
   */
  if( micro == "Microscopic" )
  {
    if     ( units == "cm^2/g." ) {
      opac_units = UT_cm2pg;
      conversion_factor = 1.0;
    }
    else
    {
      stringstream str;
      str << "Unexpected units \"" << units << "\" in PDT material data "
        << "file " << filename << ".";
      throw CommonException( str, CET_INTERNAL_ERROR );
    }
    microscopic = true;
  }
  else if( micro == "Macroscopic" )
  {
    if     ( units == "cm^-1." ) {
      opac_units = UT_icm;
      conversion_factor = 1.0;
    }
    else if( units == "m^-1." ) {
      // opac_units = UT_im;
      opac_units = UT_icm;
      conversion_factor = 0.01;
    }
    else
    {
      stringstream str;
      str << "Invalid units \"" << units << "\" in PDT material data "
        << "file " << filename << ".";
      throw CommonException( str, CET_INTERNAL_ERROR );
    }
    microscopic = false;
  }
  else
  {
    stringstream str;
    str << "PDT material data file " << filename << " is not labeled "
      << "as microscopic or macroscopic.\nIt is probably an old file; try "
      << "generating a new one and running again.";
    throw CommonException( str, CET_DATA_ERROR );
  }

  fin.getline(line, 250);
  fin.getline(line, 250);
  fin.getline(line, 250);

  for( size_t t=0; t!=num_temps; ++t )
    fin >> temps[t];
  fin.getline(line, 250);
  fin.getline(line, 250);
  fin.getline(line, 250);

  for( size_t d=0; d!=num_dens; ++d )
    fin >> densities[d];
  fin.getline(line, 250);
  fin.getline(line, 250);
  fin.getline(line, 250);

  for( size_t g=0; g!=num_groups+1; ++g )
    fin >> group_structure[g];
  fin.getline(line, 250);
  fin.getline(line, 250);
  fin.getline(line, 250);
  fin.getline(line, 250);

  resize(grey_ross, num_temps, num_dens);
  resize(grey_planck, num_temps, num_dens);
  resize(ross, num_temps, num_dens, num_groups);
  resize(planck, num_temps, num_dens, num_groups);

  for( size_t t=0; t!=num_temps; ++t )
  {
    for( size_t d=0; d!=num_dens; ++d )
    {
      for( size_t p=0; p!=n_processes; ++p )
      {
        int mt_num;
        fin >> word >> mt_num;
        fin.getline(line, 250);
        if     ( mt_num == 3001 ) {
          real8 num;
          fin >> num;
          grey_ross[t][d] = num * conversion_factor;
        }
        else if( mt_num == 3002 ) {
          real8 num;
          fin >> num;
          grey_planck[t][d] = num * conversion_factor;
        }
        else if( mt_num == 3011 ) {
          for( size_t g=0; g!=num_groups; ++g ){
            real8 num;
            fin >> num;
            ross[t][d][g] = num * conversion_factor;
          }
        }
        else if( mt_num == 3012 ) {
          for( size_t g=0; g!=num_groups; ++g ) {
            real8 num;
            fin >> num;
            planck[t][d][g] = num * conversion_factor;
          }
        }
        else {
          stringstream str; WHERE(str);
          str << "Unrecognized process identifier " << mt_num
            << " in file " << filename << ".";
          throw CommonException( str, CET_DATA_ERROR );
        }
      } // processes
      fin.getline(line, 250);
      fin.getline(line, 250);
      fin.getline(line, 250);
      fin.getline(line, 250);
    } // densities
  } // temps

  if( !fin.eof() ) { // fin.fail() ) {
    stringstream str;
    str << "Problem reading opacity data file " << filename << ".";
    throw CommonException( str, CET_DATA_ERROR );
  }
  fin.close();

}

// ============================================================================
//    Fine-group compression
// ============================================================================

void Opacities::compress(vector<vector<vector<real8> > >& fine_r,
  vector<vector<vector<real8> > >& fine_p, vector<real8>& fine_group_bounds)
{
  // ======= NOTE ===========
  // the user-specified group bounds must coincide with fine group bounds
  // ========================

  bool VALID_GROUP_STRUCTURE = true;
  for( size_t i=0; i!=num_groups+1; ++i )
  {
    bool match = false;
    for( size_t j=0; j!=fine_group_bounds.size()+1; ++j )
    {
      if( group_structure[i] == fine_group_bounds[j] )
      {
        match = true;
        break;
      }
    }
    if( !match )
      VALID_GROUP_STRUCTURE = false;
  }
  if( ! VALID_GROUP_STRUCTURE )
    throw std::invalid_argument(
      "Opacities::compress - Group boundaries must coincide with data's.");

  size_t t_index, d_index, f_index, g_index, min, max;
  Planck planckobject(1e-12);

  // all averaging will take place at specific temp-density couplets
  for (t_index=0; t_index<temps.size(); ++t_index)
  {
    for (d_index=0; d_index<densities.size(); ++d_index)
    {
      ross[t_index][d_index].resize( num_groups, 0.0 );
      planck[t_index][d_index].resize( num_groups, 0.0 );

      min = 0;
      while (fine_group_bounds[min] > group_structure[0])
        min++;
      max = min;

      // loop over all coarse groups
      for (g_index=0; g_index<num_groups; ++g_index)
      {

        // find all fine groups that fall into current group
        //   fine_group_bounds are lower bounds, and at the end
        //   of this while, max should index the first fine group
        //   of the _next_ group.
        while (fine_group_bounds[max] > group_structure[g_index+1])
          max++;

        // numerators and denominators for rosseland and planck
        real8 r_num = 0;
        real8 r_den = 0;
        real8 p_num = 0;
        real8 p_den = 0;

        for(f_index=min; f_index!=max; ++f_index)
        {
          // dBdT for the fine group
          real8 dBfdT = planckobject.integrate_dBdT(temps[t_index],
                    fine_group_bounds[f_index+1], fine_group_bounds[f_index]);
          r_num += dBfdT / fine_r[t_index][d_index][f_index];
          r_den += dBfdT;

          // B for the fine group
          real8 Bf = planckobject.integrate_B(temps[t_index],
                  fine_group_bounds[f_index+1], fine_group_bounds[f_index]);
          p_num += Bf * fine_p[t_index][d_index][f_index];
          p_den += Bf;

        } // fine groups

        // Rosseland compression
        if (r_num != 0)
          ross[t_index][d_index][g_index] = r_den / r_num;
        else
        {
          for (f_index=min; f_index!=max; ++f_index)
            r_num += 1 / fine_r[t_index][d_index][f_index];
          // assuming non-zero opacities
          ross[t_index][d_index][g_index] = (min - max) / r_num;
        }

        // Planck compression
        if (p_den != 0)
          planck[t_index][d_index][g_index] = p_num / p_den;
        else
        {
          p_num = 0;
          for (f_index=min; f_index!=max; ++f_index)
            p_num += fine_p[t_index][d_index][f_index];
          planck[t_index][d_index][g_index] = p_num / (min - max );
        }

        min = max;

      } // coarse groups
    } // densities
  } // temperatures

} // Opacities::compress()

// end of file
