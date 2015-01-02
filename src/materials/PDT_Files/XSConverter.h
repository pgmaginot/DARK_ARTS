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


// =============================================================================
//
// 01/18/2013 hayes
//
// Added support for converting ISOTXS text files to PDT format
//
  /*!
    \file   Converter.h
    \author Michael Adams, Aaron Roney, Bruno Turcksin
    \date   June 2010

    \brief
     This file describes a class (and its auxiliaries) used to convert data
   files to the standard PDT format.  The program that actually does the
   conversion is in Converter.cc.
  */
// =============================================================================

#ifndef XSCONVERTER_H
#define XSCONVERTER_H

#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <string.h>
#include <stdio.h>
#include <sstream>

#include "Opacities.h"
#include "CrossSections.h"
#include "CommonException.h"
#include "MT.h"

/*!
  \enum mf_format
  \brief An enumeration of valid material file formats.
*/
enum mf_format {
  MF_PDT, MF_TOPS, MF_SESAME, MF_CSV, MF_CX_TEST, MF_AMPX_NEUTRON,
  MF_AMPX_COUPLED, MF_MATXS, MF_CRASH, MF_CEPXS, MF_CEPXS_BFP, MF_DRAGON,
  MF_ISOTXS
};

/*!
  \enum mat_data_type
  \brief An enumeration of valid material data types.
*/
enum mat_data_type {
  MDT_THERMAL, MDT_NEUTRON, MDT_GAMMA, MDT_COUPLED_NG, MDT_ELECTRON
};

/*!
  \enum mat_data_units
  \brief An enumeration of valid material data units.
*/
enum mat_data_units {
  MDU_cmspg,      // cm^2 / g               , opacity units
  MDU_m2pkg,      // m^2 / kg               , opacity units
  MDU_im,         // m^-1                   , macroscopic opacity units
  MDU_barns,      // barn = 1e-24 cm^2      , most n0 cross sections
  MDU_cms,        // cm^2                   , for n0 test problems
  MDU_icm         // cm^-1                  , for electron transport
};

/*!
  \struct Datafile_Info
  \brief A struct containing descriptors of a material data file.
*/
struct Datafile_Info
{
  mf_format format;
  mat_data_type type;
  mat_data_units units;
  std::string sourcefile;
  bool multigroup;
  bool microscopic;
};

// Struct for split/coalesce commands
// ==================================
struct command
{
  std::string commandType;
  std::string groupType;
  int operator1;
  int operator2;
};


// =============================================================================
// =============================================================================
/*!
  \class Converter
  \brief This class inherits Opacities in order to store and rewrite opacity
    data.
*/
// =============================================================================
// =============================================================================

class Converter : public Opacities, public CrossSections
{
public:
  // Container of data sets by process identifiers ("MT" numbers)
  std::vector<MT_base*> mt_vec;

  // Numbers of opacity, neutron, gamma,electron, coupled, and transfer
  // processes
  size_t np_opacity;
  size_t np_neutron;
  size_t np_gamma;
  size_t np_electron;
  size_t np_positron;
  size_t np_coupled;
  size_t np_transfer;

  // Pointers to manage multiple inheritance issues
  std::vector<double>* temps;
  std::vector<double>* densities;
  std::vector<double>* group_structure;
  size_t* num_temps;
  size_t* num_dens;
  size_t* num_groups;

  // Populated for coupled neutron gamma data for CEPXS or for CEPXS-BFP
  std::vector<double> n_egs;
  std::vector<double> g_egs;
  std::vector<double> e_egs;
  std::vector<double> p_egs;

  // for CrossSections types
  size_t zero_density;
  std::vector<double> empty_densities;

  // Auxiliaries for Opacities objects
  void reverse(std::vector<double>& vec);
  void reverse_opacities();

  // Called by group_manipulation() to split or coalesce groups
  void group_split(const int& group, const double& energy,
    const bool& isNeutron, const bool& isCoupled);
  void group_combine(const int& group1, const int& group2,
    const bool& isNeutron, const bool& isCoupled);


public:
  Datafile_Info info;

  Converter() : Opacities(), CrossSections() {
    this->mod_const = 0;
    this->opac_store[OT_grey_ross] = this->opac_store[OT_grey_planck] = true;
    this->populate_mt_map();
    this->cx_scalar.resize(MT_N_SCALAR_COUNT);
    this->cx_single.resize(MT_1D_COUNT);
    this->cx_transfer.resize(MT_COUNT - MT_1D_COUNT);
    X_OFF = MT_1D_COUNT + 1;
    np_opacity=np_neutron=np_gamma=np_coupled=np_transfer = 0;
  }

  Converter(const Converter& other)
  : Opacities(other), CrossSections(other),
      mt_vec( other.mt_vec ),
      np_opacity( other.np_opacity ),
      np_neutron( other.np_neutron ),
      np_gamma( other.np_gamma ),
      np_electron( other.np_electron ),
      np_positron( other.np_positron ),
      np_coupled( other.np_coupled ),
      np_transfer( other.np_transfer ),
      n_egs( other.n_egs ),
      g_egs( other.g_egs ),
      e_egs( other.e_egs ),
      p_egs( other.p_egs ),
      zero_density( other.zero_density ),
      empty_densities( other.empty_densities ),
      info( other.info )
  {
    if( other.temps == &(other.CrossSections::temps ) )  set_ptrs_cx();
    if( other.temps == &(other.Opacities::temps ) )      set_ptrs_opac();
    *temps = *(other.temps);
    *densities = *(other.densities);
    *group_structure = *(other.group_structure);
    *num_temps = *(other.num_temps);
    *num_dens = *(other.num_dens);
    *num_groups = *(other.num_groups);
  }

  ~Converter() {
    for( std::vector<MT_base*>::iterator mt=mt_vec.begin(); mt!=mt_vec.end(); ++mt)
    { if( *mt!=NULL )  delete *mt; }
  }

  // Allows user to split and coalesce energy groups
  void group_manipulation(const std::string& cmd = "");

  //    PDT format write and read
  // ===============================
  void write_pdt(std::string filename);
  void read_pdt(std::string filename);
  void read_pdt_opac(std::string filename);
  void read_pdt_cx(std::string filename);

  //    Opacities functionality
  // =============================
  void set_ptrs_opac();

  /// CRASH format
  void read_crash(std::string filename, double opac_max);

  /// TOPS format
  void read_tops(std::string filename, double opac_max);
  void if_warnings(std::string file, std::vector<std::vector<std::vector<double> > >& ross,
    std::vector<std::vector<std::vector<double> > >& planck);
  void fix_energy(double e_max) {
    if( group_structure->back() >= e_max ) {
      std::stringstream str;
      str << "Highest group upper bound must be greater than lower bound. "
          << "Check file and verify your data.";
      throw CommonException( str, CET_INPUT_ERROR );
    }
    group_structure->push_back( e_max );
    reverse(*group_structure);
  }

  /// SESAME format
  void read_sesame(std::string filename, double opac_max);

  /// For one particular file, kind of pointless
  void read_csv(std::string filename);

  //    CrossSections functionality
  // =================================
  void set_ptrs_cx();

  /// AMPX neutron
  void read_ampx_neutron(std::string filename, float clyde);

  /// AMPX coupled
  void read_ampx_coupled(std::string filename, float clyde);

  // AMPX auxiliaries
  void skip_record(std::fstream& ampx);
  void read_record_nine(std::fstream& ampx, int num_processes, int num_groups);
  void read_record_three(std::fstream& ampx);
  void read_scat_matrix(std::fstream& ampx);

  /// MATXS format (NJOY data)
  void read_matxs(std::string filename, std::string pdtfile, int opt_combine);
  std::map<std::string, size_t> get_hollerith_1d_map();
  std::map<std::string, size_t> get_hollerith_trans_map();

  /// ISOTXS format
  void read_isotxs(std::string filename, std::string pdtfile);

  /// DRAGON format
  void read_dragon(std::string filename, std::string pdtfile);
  void read_edi(std::ifstream& fin, std::string pdtfile);

  /// CEPXS format
  void read_cepxs(std::string filename);

  /// CEPXS-BFP format
  void read_cepxs_bfp(std::string filename);
  /// Used by read_cepxs_bfp
  void read_energy_boundary_cepxs_bfp(std::ifstream &fin);
  /// Used by read_cepxs_bfp
  void read_xs_cepxs_bfp(std::ifstream &fin, int mom, int length_number,
    int n_elemts);
  /// Transform char in doubles
  std::vector<double> read_line_xs(const char* data, int length_number);
  /// Compute particle type, energy group and type of cross section.
  std::vector<int> index(int i, int j, int n_elemts);
  /// Compute the index for the scattering cross sections
  std::vector<int> scattering_index(int row, int column, int particle);
  /// Change the format of the cross section to the MT_base format
  void store(bool BFP);
  /// Total number of groups in CEPXS or CEPSX-BFP
  int n_groups;
  /// Number of entries in CEPXS-BFP
  int n_entries;
  /// Position of the self-scattering cross section in CEPXS-BFP
  int self_scatter_position;
  /// Anisotropy order in CEPXS or CEPXS-BFP
  int L_max;
  /// Number of groups for each particle type (gamma, electron, positron)
  std::vector<int> n_groups_particles;
  /// Number of types of particles in CEPXS or CEPXS-BFP
  int n_particles;
  /// Energy cross sections for CEPXS or CEPSX-BFP
  std::vector<std::vector<double> > energy_xs;
  /// Momentum transfer for CEPXS-BFP
  std::vector<std::vector<double> > alpha;
  /// Stopping power for CEPXS-BFP
  std::vector<std::vector<double> > stopping_power;
  /// Absorption cross section for CEPXS-BFP
  std::vector<std::vector<double> > absorption_xs;
  /// Total cross section for CEPXS or CEPXS-BFP
  std::vector<std::vector<double> > total_xs;
  /// Scattering cross section for CEPXS-BFP
  std::vector<std::vector<std::vector<double> > > scattering_xs;

  /// Test problem format
  void read_test(std::string filename);
};

#endif
