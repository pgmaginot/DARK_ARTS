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
  \file BaseComponent.h
  \brief Offers a basic component representation.

  All component types will be derived from this one.  Functionality will be
  added in a hierarchical way on the derivation DAG.

  This class simply provides an interface to material data.  the most common
  calls will be get_sigt(), get_transfer_cx(), and get_opacities; more specific
  calls can be made by MT number.  CX_MT_id is an enum type declared in
  Common/CrossSections.h.
*/

#ifndef _BaseComponent_h
#define _BaseComponent_h

// System includes
#include <vector>
#include <iostream>
#include <fstream>
#include <cassert>
#include <bitset>

// PDT includes
#include "CrossSections.h"
#include "Opacities.h"
#include "Edits.h"
#include "CommonException.h"
#include "BaseDepletionEvent.h"

enum BC_TYPE
{
  BCT_VACUUM, BCT_CX, BCT_OPAC
};

struct SpecificHeat
{
  bool valid;
  real8 A, B, C;
  UNIT_TYPE cv_units;

  void define_type(stapl::typer& t);

  void set(real8 a, real8 b, real8 c)  {
    A = a;  B = b;  C = c;  valid = true;
  };

  real8 get_Cv(real8 T) {
    if( ! valid )
    {
      std::stringstream str;
      str << "Uninitialized value for specific heat.";
      throw CommonException(str, CET_DATA_ERROR);
    }
    return A + B * pow(T, C);
  };

};


/*!
  \class BaseComponent
  \brief The base class for all types of components.
*/

class BaseComponent
{
protected:

  // The identifier for the component
  std::string id;

  // The type of the component
  BC_TYPE data_type;

  // Global component ID
  int gID;

  // Atomic and element number
  int Anum, Znum;
  bool dontDepleteMe;

  int num_groups;
  int scat_order;
  UNIT_TYPE md_units;
  bool md_microscopic;
  bool md_fissionable;
  bool md_derivedMTAbs;

  // CrossSections object handles all cross sections (CrossSections.h, .cc)
  CrossSections* crosssections;

  // Opacities object handles all opacities (Opacities.h, .cc)
  Opacities* opacities;

  // Analytic specific heat representation
  SpecificHeat Cv;

  // Vector of sensitivity of the QOI to SCALAR parameters
  std::vector<real8> dQOI_dSCALAR;

  // Vector of sensitivity of the QOI to transfer parameters
  // sig(g'->g,m) = dQOI_dTRANSFER[param][g][m][g']
  std::vector<std::vector<std::vector<std::vector<real8> > > > dQOI_dTRANSFER;

  // Vector of sensitivity of the QOI to ONED parameters
  std::vector<std::vector<real8> > dQOI_dONED;

  // A list of depletion reactions that result in a production of this component
  std::vector<ReactionEvent> reactionEvents;

  // A list of decay reactions that result in a production of this component
  std::vector<DecayEvent> decayEvents;

public:

  BaseComponent(std::string _id, int scat, SpecificHeat _Cv)
    : id(_id), gID(0), Anum(0), Znum(0), dontDepleteMe(false), num_groups(0),
      scat_order(scat), md_units(UT_barns), md_microscopic(true),
      md_fissionable(false), md_derivedMTAbs(false), crosssections(NULL),
      opacities(NULL), Cv(_Cv)
  {}

  BaseComponent(const BaseComponent& other, std::string _id) : id(_id)
  {
    if( this != &other )
    {
      gID                = other.gID;
      Anum               = other.Anum;
      Znum               = other.Znum;
      dontDepleteMe      = other.dontDepleteMe;
      *crosssections     = *(other.crosssections);
      *opacities         = *(other.opacities);
      Cv                 = other.Cv;
      num_groups         = other.num_groups;
      scat_order         = other.scat_order;
      md_microscopic     = other.md_microscopic;
      md_fissionable     = other.md_fissionable;
      md_derivedMTAbs    = other.md_derivedMTAbs;
      md_units           = other.md_units;
      dQOI_dSCALAR       = other.dQOI_dSCALAR;
      dQOI_dTRANSFER     = other.dQOI_dTRANSFER;
      dQOI_dONED         = other.dQOI_dONED;
      reactionEvents     = other.reactionEvents;
      decayEvents        = other.decayEvents;
    }
  };

  BaseComponent& operator=(const BaseComponent& other)
  {
    if( this != &other )
    {
      id                 = other.id;
      gID                = other.gID;
      Anum               = other.Anum;
      Znum               = other.Znum;
      dontDepleteMe      = other.dontDepleteMe;
      crosssections      = other.crosssections;
      opacities          = other.opacities;
      Cv                 = other.Cv;
      num_groups         = other.num_groups;
      scat_order         = other.scat_order;
      md_microscopic     = other.md_microscopic;
      md_fissionable     = other.md_fissionable;
      md_derivedMTAbs    = other.md_derivedMTAbs;
      md_units           = other.md_units;
      dQOI_dSCALAR       = other.dQOI_dSCALAR;
      dQOI_dTRANSFER     = other.dQOI_dTRANSFER;
      dQOI_dONED         = other.dQOI_dONED;
      reactionEvents     = other.reactionEvents;
      decayEvents        = other.decayEvents;
    }
    return *this;
  };

  virtual ~BaseComponent() {}

  void define_type(stapl::typer& t);

  // Access functions
  // ================
  inline std::string get_component_id()            { return id; }
  inline void set_component_id(std::string _id)    { id = _id; }
  inline BC_TYPE get_type()                        { return data_type; }
  inline std::string get_cx_type()                 { return crosssections->get_type(); }
  inline int get_num_groups()                      { return num_groups; }
  inline int get_scat_order()                      { return scat_order; }
  inline UNIT_TYPE get_units()                     { return md_units; }
  inline UNIT_TYPE get_Cv_units()                  { return Cv.cv_units; }
  inline bool is_microscopic()                     { return md_microscopic; }
  inline bool has_fission()                        { return md_fissionable; }
  inline bool has_derivedMTAbs()                   { return md_derivedMTAbs; };

  real8 get_Cv(real8 T) {return Cv.get_Cv(T); };

  virtual void set_group_structure(std::vector<real8> g_s) = 0;

  virtual void dump_cross_sections(std::ostream& os) = 0;

  // Virtual functions for depletion data
  // ===========================

  // Function to write depletion sources
  virtual void writeDeplEvents(std::ostream& os) = 0;

  // Functions to set and return the decay constant
  virtual real8 get_lambda() = 0;

  // Function to return energy per fission
  virtual real8 get_EperFission() = 0;

  // Functions to set and return global component ID
  virtual void set_gID(int id) = 0;
  virtual int get_gID() = 0;

  // Functions to set and return Z and A numbers
  virtual void set_Anum(int A) = 0;
  virtual void set_Znum(int Z) = 0;
  virtual int get_Anum() = 0;
  virtual int get_Znum() = 0;

  // Functions to set and return the bool which says whether to deplete me or not
  virtual bool get_dontDepleteMe() = 0;
  virtual void set_dontDepleteMe(bool doMe) = 0;

  // Functions to set store flags for cross-sections and to
  // call data integrity checks
  virtual bool checkStoreTrue(CX_MT_id process) = 0;
  virtual void checkDepletionData() = 0;

  // Virtual retrieval functions
  // ===========================
  virtual void get_siga(std::vector<real8>& vec, real8 T, real8 m_density) = 0;
  virtual void get_sigf(std::vector<real8>& vec, real8 T, real8 m_density) = 0;
  virtual void get_opacities(std::vector<real8>& vec, real8 temperature,
    real8 density, OPAC_STATUS& t_err, OPAC_STATUS& d_err,
    opac_type type = OT_COUNT) = 0;
  virtual void get_sigt(std::vector<real8>& vec, real8 T, real8 m_density) = 0;
  virtual void get_nusigf(std::vector<real8>& vec, real8 T=300.0) = 0;
  virtual void get_nubar(std::vector<real8>& vec, real8 T=300.0) = 0;
  virtual void get_chi(std::vector<real8>& vec, real8 T=300.0) = 0;
  virtual void get_specific_cx(std::vector<real8>& vec, CX_MT_id mt, real8 T,
    real8 m_density) = 0;
  /// Get the stopping power used for electron-bfp transport
  virtual void get_stopping_power(std::vector<real8>& vec) = 0;
  /// Get the width of the energy groups
  virtual void get_energy_width(std::vector<real8>& vec) = 0;
  /// get the energy deposition cross section needed to compute the dose
  virtual void get_energy_deposition(std::vector<real8>& vec) = 0;
  virtual void get_transfer_cx(std::vector<std::vector<std::vector<real8> > >& vec,
    CX_MT_id MT, real8 T = 300) = 0;
  virtual void get_transfer_cx(std::vector<std::vector<std::vector<real8> > >& vec,
    CX_MT_id MT, size_t g_min, size_t g_max, real8 T = 300) = 0;
  virtual void get_transfer_cx(std::vector<std::vector<std::vector<real8> > >& vec,
    CX_MT_id MT, size_t n_mom, real8 T = 300) = 0;

  // Return a reference to dQOI_dSCALAR
  virtual std::vector<real8>& dQdSCALAR() = 0;

  // Return a reference to dQOI_dTRANSFER
  virtual std::vector<std::vector<std::vector<std::vector<real8> > > >& dQdTRANSFER() = 0;

  // Return a reference to dQOI_dONED
  virtual std::vector<std::vector<real8> >& dQdONED() = 0;

  // Return access to the depletion events
  virtual std::vector<ReactionEvent>& rxEvent() = 0;
  virtual std::vector<DecayEvent>& dkEvent() = 0;

  // Compute the production rate of a component
  virtual real8 computeRxProdRate(size_t ind, std::vector<real8>& phi) = 0;

  // Retrieve the XS associated with a particular depletion reaction
  virtual void retrieveRxXs(size_t ind, std::vector<real8>& sig) = 0;
};

// =============================================================================
// VacuumComponent Class
// =============================================================================

class VacuumComponent : public BaseComponent
{
public:

  VacuumComponent(const std::string _id, SpecificHeat _Cv, int _scat_order,
    std::bitset<EDITS_SIZE>* _edits): BaseComponent(_id, _scat_order, _Cv)
  {
    Cv.A = 1.0e-30;  Cv.B = 0.0;  Cv.C = 0.0;  Cv.valid = true;
    md_microscopic = false;
    this->data_type = BCT_VACUUM;
  };

  ~VacuumComponent() {};

  void set_group_structure(std::vector<real8> g_s) {
    num_groups = g_s.size() - 1;
  }

  void dump_cross_sections(std::ostream& os) {
    os << "Component " << id << " is a vacuum; all cross sections are zero."
      << std::endl;
  }

  void get_opacities(std::vector<real8>& vec, real8 temperature, real8 density,
    OPAC_STATUS& t_err, OPAC_STATUS& d_err, opac_type type = OT_COUNT) {
    vec.clear();  vec.resize( num_groups, 1.0e-30 );
  }

  void writeDeplEvents(std::ostream& os) {
    std::stringstream str;
    str << "Vacuum component should not be trying to write depletion events.";
    throw CommonException(str, CET_TYPE_ERROR);
  }

  virtual void get_sigt(std::vector<real8>& vec, real8 T, real8 m_density) {
    vec.clear();  vec.resize( num_groups, 1.0e-30 );
  }

  virtual void get_siga(std::vector<real8>& vec, real8 T, real8 m_density) {
    vec.clear();  vec.resize( num_groups, 1.0e-30 );
  }

  virtual void get_sigf(std::vector<real8>& vec, real8 T, real8 m_density) {
    vec.clear();  vec.resize( num_groups, 1.0e-30 );
  }

  virtual void get_nusigf(std::vector<real8>& vec, real8 T=300.0) {
    vec.clear();  vec.resize( num_groups, 1.0e-30 );
  }

  virtual void get_nubar(std::vector<real8>& vec, real8 T=300.0) {
     vec.clear();  vec.resize( num_groups, 1.0e-30 );
  }

  virtual void get_chi(std::vector<real8>& vec, real8 T=300.0) {
    vec.clear();  vec.resize( num_groups, 1.0e-30 );
  }

  virtual void get_specific_cx(std::vector<real8>& vec, CX_MT_id MT, real8 T=300.0,
                               real8 density=0.0) {
    vec.clear();  vec.resize( num_groups, 1.0e-30 );
  }

  virtual void get_stopping_power(std::vector<real8>& vec) {
    vec.clear();  vec.resize( num_groups, 1.0e-30 );
  }

  virtual void get_energy_width(std::vector<real8>& vec) {
    vec.clear();  vec.resize( num_groups, 1.0e-30 );
  }

  virtual void get_energy_deposition(std::vector<real8>& vec) {
    vec.clear();  vec.resize( num_groups, 1.0e-30 );
  }

  void get_transfer_cx(std::vector<std::vector<std::vector<real8> > >& vec,
    CX_MT_id MT, real8 T = 300)
  {
    vec.clear();  vec.resize( num_groups, std::vector<std::vector<real8> >
      ( scat_order+1, std::vector<real8>( num_groups, 0.0 ) ) );
  }

  void get_transfer_cx(std::vector<std::vector<std::vector<real8> > >& vec,
    CX_MT_id MT, size_t g_min, size_t g_max, real8 T = 300)
  {
    size_t ng = g_max - g_min;
    vec.clear();  vec.resize( ng, std::vector<std::vector<real8> >
      ( scat_order+1, std::vector<real8>( num_groups, 0.0 ) ) );
  }

  void get_transfer_cx(std::vector<std::vector<std::vector<real8> > >& vec,
    CX_MT_id MT, size_t n_mom, real8 T = 300)
  {
    vec.clear();  vec.resize( num_groups, std::vector<std::vector<real8> >
      ( n_mom-1, std::vector<real8>( num_groups, 0.0 ) ) );
  }

  real8 get_lambda() {
    std::stringstream str;
    str << "Vacuum Component tried to access decay constant";
    throw CommonException(str, CET_TYPE_ERROR);
  }

  real8 get_EperFission() {
    std::stringstream str;
    str << "Vacuum Component tried to access energy per fission";
    throw CommonException(str, CET_TYPE_ERROR);
  }

  void set_gID(int id) { gID=id;}
  int get_gID() {return gID;}

  void set_Anum(int A) { Anum=0;}
  void set_Znum(int Z) { Znum=0;}
  int get_Anum() {return Anum;}
  int get_Znum() {return Znum;}

  void set_dontDepleteMe(bool doMe) {dontDepleteMe = doMe;}
  bool get_dontDepleteMe() {return dontDepleteMe;}
  bool checkStoreTrue(CX_MT_id process) { return false;}
  void checkDepletionData() {}

  std::vector<real8>& dQdSCALAR() {
    std::stringstream str;
    str << "Vacuum component attempted to access dQOI_dSCALAR";
    throw CommonException(str, CET_TYPE_ERROR);
  }

  std::vector<std::vector<std::vector<std::vector<real8> > > >& dQdTRANSFER() {
    std::stringstream str;
    str << "Vacuum component attempted to access dQOI_dTRANSFER";
    throw CommonException(str, CET_TYPE_ERROR);
  }

  std::vector<std::vector<real8> >& dQdONED() {
    std::stringstream str;
    str << "Vacuum component attempted to access dQOI_dONED";
    throw CommonException(str, CET_TYPE_ERROR);
  }

  std::vector<ReactionEvent>& rxEvent() {
    std::stringstream str;
    str << "Vacuum component attempted to access list "
        << "of reaction depletion events";
    throw CommonException(str, CET_TYPE_ERROR);
    return reactionEvents;
  }

  std::vector<DecayEvent>& dkEvent() {
    std::stringstream str;
    str << "Vacuum component attempted to access list "
        << "of decay depletion events";
    throw CommonException(str, CET_TYPE_ERROR);
    return decayEvents;
  }

  real8 computeRxProdRate(size_t ind, std::vector<real8>& phi) {
    std::stringstream str;
    str << "Vacuum component attempted to call the depletion related "
        << "routine computeRxProdRate(...),";
    throw CommonException(str, CET_TYPE_ERROR);
    return 0.0;
  }

  void retrieveRxXs(size_t ind, std::vector<real8>& sig) {
    std::stringstream str;
    str << "Vacuum component attempted to call the depletion related "
        << "routine retrieveRxXs(...),";
    throw CommonException(str, CET_TYPE_ERROR);
  }
};

// =============================================================================
// NeutronicsComponent Class
// =============================================================================

class CXComponent : public BaseComponent
{
public:
  using BaseComponent::crosssections;

  CXComponent(const std::string _id, SpecificHeat _Cv, int _scat_order,
    std::bitset<EDITS_SIZE>* _edits): BaseComponent(_id, _scat_order, _Cv)
  {
    crosssections = new CrossSections(_id, _scat_order, _edits);
    this->data_type = BCT_CX;
  };

  void initialize_neutron(std::string filename, bool isDpl, bool isKeig) {
    crosssections->initialize_neutron(filename, isDpl, isKeig);
    this->set_values();
  }

  void initialize_coupled(std::string filename) {
    crosssections->initialize_coupled(filename);
    this->set_values();
  }

  void initialize_gamma(std::string filename) {
    crosssections->initialize_gamma(filename);
    this->set_values();
  }

  void initialize_electron(std::string filename, bool bfp,
    std::string const &reorder)
  {
    crosssections->initialize_electron(filename, bfp, reorder);
    this->set_values();
  }

  void initialize_k_eigen(std::string filename, bool isDpl, bool isKeig) {
    crosssections->initialize_k_eigen(filename, isDpl, isKeig);
    this->set_values();
  }

  void set_values() {
    md_units = crosssections->get_units();
    md_microscopic = crosssections->is_microscopic();
    md_fissionable = crosssections->has_fission();
    md_derivedMTAbs = crosssections->has_derivedMTAbs();
  }

  void set_group_structure(std::vector<real8> g_s) {
    crosssections->set_group_structure(g_s);
    this->num_groups = g_s.size() - 1;
  }

  ~CXComponent() { if (crosssections) delete crosssections; };

  // Inherited virtual
  // =================
  void get_sigt(std::vector<real8>& vec, real8 T = 300, real8 m_density = -1) {
    crosssections->get_cross_sections(vec, MT_total, T);
  };

  void get_siga(std::vector<real8>& vec, real8 T = 300, real8 m_density = -1) {
    crosssections->get_cross_sections(vec, MT_absorption, T);
  };

  void get_sigf(std::vector<real8>& vec, real8 T = 300, real8 m_density = -1) {
    crosssections->get_cross_sections(vec, MT_fission, T);
  };

  void get_nusigf(std::vector<real8>& vec, real8 T = 300) {
    crosssections->get_cross_sections(vec, MT_nu_sig_f, T);
  };

  void get_nubar(std::vector<real8>& vec, real8 T = 300) {
    crosssections->get_cross_sections(vec, MT_nubar, T);
  };

  void get_chi(std::vector<real8>& vec, real8 T = 300) {
    crosssections->get_cross_sections(vec, MT_chi, T);
  };

  /// Inherited virtual function
  void get_stopping_power(std::vector<real8>& vec) {
    crosssections->get_cross_sections(vec, MT_gamma_stopping_power,1);
  };

  /// Inherited vitual function
  void get_energy_width(std::vector<real8>& vec) {
    crosssections->get_energy_width(vec);
  };

  /// Inherited virtual function
  void get_energy_deposition(std::vector<real8>& vec) {
    crosssections->get_cross_sections(vec, MT_gamma_energy_deposition,1);
  };

  void get_opacities(std::vector<real8>& vec, real8 T, real8 rho,
    OPAC_STATUS& t_err, OPAC_STATUS& d_err, opac_type type)
  {
    std::stringstream str;
    str << "Opacities are not available for neutronics.";
    throw CommonException(str, CET_TYPE_ERROR);
  };

  void dump_cross_sections(std::ostream& os) {
    crosssections->print_all_cross_sections(os);
  };

  void writeDeplEvents(std::ostream& os)
  {
    os << "Depletion data for component "
       << this->get_component_id() << ":" << std::endl;

    if( dontDepleteMe)
      os << "\tThis nuclide WILL NOT be depleted" << std::endl
         << "\tGlobal Component ID: " << gID << std::endl
         << "\tZ-Number: " << Znum << std::endl
         << "\tA-Number: " << Anum << std::endl
         << "\tDecay Constant, lambda=" << this->get_lambda() << "s^-1."
         << std::endl;

    if (md_derivedMTAbs)
      os << "\tThe removal cross-section has been derived from the "
         << "\n\t\ttotal and scattering cross-sections." << std::endl;
    else
      os << "\tThe removal cross-section has been input explicitly "
         << "\n\t\tor was not derived from available data." << std::endl;

    if (crosssections->checkStoreTrue(MT_E_per_fission))
      os << "\tRecoverable energy per Fission: "
         << this->get_EperFission() * 1.0e-6 << " MeV." << std::endl;
    else
      os << "\tEnergy Release per Fission: n/a." << std::endl;

    size_t nrx = reactionEvents.size();
    for(size_t rxn=0; rxn<nrx; rxn++) {
      if (reactionEvents[rxn].get_MT_id() == MT_absorption)
        continue;
      if (reactionEvents[rxn].is_active())
        os << "\tChild component, "
           << reactionEvents[rxn].get_childID()
           << ", is formed via the "
           << reactionEvents[rxn].get_eventID()
           << " reaction with branching ratio "
           << reactionEvents[rxn].get_branchingRatio()
           << "." << std::endl;
    }

    size_t ndk = decayEvents.size();
    for (size_t dk=0; dk<ndk; dk++) {
      if (decayEvents[dk].is_active())
        os << "\tChild component, "
           << decayEvents[dk].get_childID()
           << ", is formed via decay with branching ratio "
           << decayEvents[dk].get_branchingRatio()
           << "." << std::endl;
    }
  }

  real8 get_lambda()
  {
    real8 lambda = 0.0;
    if (crosssections->checkStoreTrue(MT_half_life))
    {
      std::vector<real8> t12;
      crosssections->get_cross_sections(t12, MT_half_life, 1.0);
      lambda = NATURAL_LOG_2 / t12[0];
    }
    return lambda;
  }

  real8 get_EperFission()
  {
    std::vector<real8> epf;
    crosssections->get_cross_sections(epf, MT_E_per_fission, 1.0);
    return epf[0];
  }

  void set_gID(int id) { gID = id;}

  int get_gID() {return gID;}

  void set_Anum(int A) { Anum = A;}

  void set_Znum(int Z) { Znum = Z;}

  int get_Anum() {return Anum;}

  int get_Znum() {return Znum;}

  void set_dontDepleteMe(bool doMe) {dontDepleteMe = doMe;}

  bool get_dontDepleteMe() {return dontDepleteMe;}

  bool checkStoreTrue(CX_MT_id process) {
    return crosssections->checkStoreTrue(process);
  }

  void checkDepletionData()
  {
    Logger msg;
    //Set the flag to store the absoprtion cross-section to true -- its required
    crosssections->setStoreTrue(MT_absorption);
    //Check to see if the supported reactions have the correct cross-sections
    size_t numRx = reactionEvents.size();
    for (size_t rx=0; rx<numRx; rx++)
    {
      CX_MT_id thisID = reactionEvents[rx].get_MT_id();
      if (!crosssections->checkStoreTrue(thisID))
      {
        msg.log(WARN) << "PDT detected that components "
                      << reactionEvents[rx].get_parentID()
                      << " and "
                      << reactionEvents[rx].get_childID()
                      << " are a parent/child pair for depletion event "
                      << reactionEvents[rx].get_eventID()
                      << ",\n\t but PDT did not find the "
                      << reactionEvents[rx].get_eventID()
                      << " cross-sections for parent nuclide "
                      << reactionEvents[rx].get_parentID()
                      << "\n\tPDT will not process the reaction."
                      << std::endl;
        reactionEvents[rx].set_active(false);
      }
    }
    std::string tmp("Depletion data from input");
    crosssections->call_check_data(tmp);
  }

  std::vector<real8>& dQdSCALAR() {
    return dQOI_dSCALAR;
  }

  std::vector<std::vector<std::vector<std::vector<real8> > > >& dQdTRANSFER(){
    return dQOI_dTRANSFER;
  }

  std::vector<std::vector<real8> >& dQdONED(){
    return dQOI_dONED;
  }

  std::vector<ReactionEvent>& rxEvent() {return reactionEvents;}

  std::vector<DecayEvent>& dkEvent() {return decayEvents;}

  real8 computeRxProdRate(size_t rxInd, std::vector<real8>& phi)
  {
    size_t nGroup = phi.size();
    std::vector<real8> sig;
    get_specific_cx(sig,reactionEvents[rxInd].get_MT_id());
    real8 prodRate = 0.0;
    for (size_t g=0; g<nGroup; g++)
      prodRate += phi[g] * sig[g];
    prodRate *= reactionEvents[rxInd].get_branchingRatio();
    // Special case for absorption
    if (reactionEvents[rxInd].get_MT_id() == MT_absorption)
      prodRate *= -1.0;
    return prodRate;
  }

  void retrieveRxXs(size_t rxInd, std::vector<real8>& sig)
  {
    get_specific_cx(sig,reactionEvents[rxInd].get_MT_id());
  }

  // Specific to NeutronicsComponent
  // ===============================
  void get_specific_cx(std::vector<real8>& vec, CX_MT_id MT, real8 T = 300.,
    real8 density = 0.)
  {
    crosssections->get_cross_sections(vec, MT, T);
  };

  void get_transfer_cx(std::vector<std::vector<std::vector<real8> > >& vec,
    CX_MT_id MT, real8 T = 300)
  {
    crosssections->get_transfer_cross_sections(vec, MT, T);
  };

  void get_transfer_cx(std::vector<std::vector<std::vector<real8> > >& vec,
    CX_MT_id MT, size_t g_min, size_t g_max, real8 T = 300)
  {
    crosssections->get_transfer_cross_sections(vec, MT, T, g_min, g_max);
  };

  void get_transfer_cx(std::vector<std::vector<std::vector<real8> > >& vec,
    CX_MT_id MT, size_t n_mom, real8 T = 300)
  {
    crosssections->get_transfer_cross_sections(vec, MT, T, n_mom-1);
  };
};

// =============================================================================
// OpacitiesComponent Class
// =============================================================================

class OpacitiesComponent : public BaseComponent
{
public:
  using BaseComponent::opacities;

  OpacitiesComponent(const std::string _id, SpecificHeat _Cv,
                     std::bitset<EDITS_SIZE>* _edits)
    : BaseComponent(_id, 0, _Cv) // 0 scat_order
  {
    opacities = new Opacities(_id);
    this->data_type = BCT_OPAC;
  };

  void set_model(real8 mod_const) { opacities->set_model(mod_const); };
  void set_abc(real8 a, real8 b, real8 c) { opacities->set_abc(a,b,c); };
  void initialize(std::string filename, opac_type _def) {
    opacities->initialize(filename, _def);
    md_units = opacities->get_units();
    md_microscopic = opacities->is_microscopic();
  }
  void set_group_structure(std::vector<real8> g_s) {
    opacities->set_group_structure(g_s);
    this->num_groups = opacities->get_num_groups();
  }

  ~OpacitiesComponent() { if (opacities) delete opacities; };

  // Inherited virtual
  // =================
  void get_sigt(std::vector<real8>& vec, real8 T, real8 m_density) {
    OPAC_STATUS t, d;
    opacities->get_opacities(vec, T, m_density, t, d);
  }

  void get_siga(std::vector<real8>& vec, real8 T, real8 m_density) {
    OPAC_STATUS t, d;
    opacities->get_opacities(vec, T, m_density, t, d);
  }

  void get_opacities(std::vector<real8>& vec, real8 T, real8 m_density, OPAC_STATUS& t_err,
        OPAC_STATUS& d_err, opac_type type = OT_COUNT)
  {
    opacities->get_opacities(vec, T, m_density, t_err, d_err, type);
  };

  void dump_cross_sections(std::ostream& os) {
    opacities->full_write(os);
  };

  void set_gID(int id) {gID=id;}
  int get_gID() {return gID;}

  void set_Anum(int A) {Anum=0;}
  void set_Znum(int Z) {Znum=0;}
  int get_Anum() {return Anum;}
  int get_Znum() {return Znum;}

  bool checkStoreTrue(CX_MT_id process) {return false;}
  void checkDepletionData() {}

  // Invalid for OpacitiesComponent
  // ==============================
  void invalid(const char* name) {
    std::stringstream str;
    str << name << " is not applicable in thermal radiation.";
    throw CommonException( str, CET_TYPE_ERROR );
  }

  void not_yet(const char* name) {
    std::stringstream str;
    str << name << " is not yet implemented for thermal radiation.";
    throw CommonException( str, CET_TYPE_ERROR );
  }

  /// Inherited virtual function
  void get_sigf(std::vector<real8>& vec, real8 T, real8 m_density) {
    invalid("Fission");
  }

  void get_energy_deposition(std::vector<real8>& vec) {
    invalid("Energy deposition");
  }

  void get_nusigf(std::vector<real8>& vec, real8 T) {
    invalid("Nu-sig_f");
  }

  void get_nubar(std::vector<real8>& vec, real8 T) {
    invalid("nubar");
  };

  void get_chi(std::vector<real8>& vec, real8 T) {
    invalid("Chi");
  }

  void get_specific_cx(std::vector<real8>& vec, CX_MT_id mt, real8 T,
    real8 m_density)
  {
    invalid("Acessing cross sections by mt id");
  }

  void get_stopping_power(std::vector<real8>& vec) {
    invalid("Stopping power");
  };

  void get_energy_width(std::vector<real8>& vec) {
    invalid("Energy width");
  };

  void get_transfer_cx(std::vector<std::vector<std::vector<real8> > >& vec,
    CX_MT_id MT, size_t g_min, size_t g_max, real8 T = 300)
  {
    not_yet("Scattering");
  };

  void get_transfer_cx(std::vector<std::vector<std::vector<real8> > >& vec,
    CX_MT_id MT, real8 T = 300)
  {
    not_yet("Scattering");
  };

  void get_transfer_cx(std::vector<std::vector<std::vector<real8> > >& vec,
    CX_MT_id MT, size_t n_mom, real8 T = 300)
  {
    not_yet("Scattering");
  };

  real8 get_lambda() {
    invalid("Access of decay constant");
    return 0.0;
  }

  real8 get_EperFission() {
    invalid("Access of energy per fission");
    return 0.0;
  }

  void set_dontDepleteMe(bool doMe) {dontDepleteMe = doMe;}

  bool get_dontDepleteMe() {return dontDepleteMe;}

  void writeDeplEvents(std::ostream& os) {
    std::stringstream str;
     str << "Opacity Component should not be trying to write depletion events.";
     throw CommonException(str, CET_TYPE_ERROR);
  }

  std::vector<real8>& dQdSCALAR() {
    invalid("Access of dQOI_dSCALAR");
    return dQOI_dSCALAR;
  }

  std::vector<std::vector<std::vector<std::vector<real8> > > >& dQdTRANSFER() {
    invalid("Access of dQOI_dTRANSFER");
    return dQOI_dTRANSFER;
  }

  std::vector<std::vector<real8> >& dQdONED( ){
    invalid("Access of dQOI_dONED");
    return dQOI_dONED;
  }

  std::vector<ReactionEvent>& rxEvent() {
    std::stringstream str;
    str << "Opacity component attempted to access list "
        << "of reaction depletion events";
    throw CommonException(str, CET_TYPE_ERROR);
    return reactionEvents;
  }

  std::vector<DecayEvent>& dkEvent() {
    std::stringstream str;
    str << "Opacity component attempted to access list "
        << "of decay depletion events";
    throw CommonException(str, CET_TYPE_ERROR);
    return decayEvents;
  }

  real8 computeRxProdRate(size_t ind, std::vector<real8>& phi) {
    std::stringstream str;
    str << "Opacity component attempted to call the depletion related "
        << "routine computeRxProdRate(...),";
    throw CommonException(str, CET_TYPE_ERROR);
    return 0.0;
  }

  void retrieveRxXs(size_t ind, std::vector<real8>& sig) {
    std::stringstream str;
    str << "Opacity component attempted to call the depletion related "
        << "routine retrieveRxXs(...),";
    throw CommonException(str, CET_TYPE_ERROR);
  }
};

#endif
