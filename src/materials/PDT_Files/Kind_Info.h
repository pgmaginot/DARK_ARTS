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

#ifndef _KIND_INFO_H
#define _KIND_INFO_H

#include "TD_Info.h"
#include "BP_Info.h"
#include "QOI_Info.h"
#include "KEFF_Info.h"
#include "AdjointController.h"

enum PROBLEM_KIND
{
  NEUTRONICS,
  RADIATIVE_TRANSFER,
  NEUTRON_PHOTON,
  GAMMA,
  ELECTRON,
  K_EIGEN
};

enum OP_TYPE
{
  ARD,
  WGS_PHI,
  AGS_PHI,
  GDA_CFEM,
  GDA_CDFEM
};

enum SPATIAL_METHOD
{
  SDM_WDD,             // Weighted diamond difference
  SDM_CBSTEP,          // Corner balance
  SDM_PWLD,            // Piece-wise linear discontinuous finite element
  SDM_BILD,            // Bilinear discontinuous finite element for XY
  SDM_TRILD,           // Trilinear discontinuous finite element for XYZ
  SDM_UNSTRUCT_BILD,   // Unlumped, bilienar DFEM for unstructured xy geometry
  SDM_UNSTRUCT_BCSZ,   // bilinear DFEM with nonlinear, consistent set to zero  
  SDM_INVALID          // Invalid spatial discretization type
};

enum FEM_TYPE
{
  FEM_UNLUMP,          // Unlumped DFEM
  FEM_LUMP,            // Lumped DFEM
  FEM_INVALID          // Invalid FEM type
};

enum GRID_TYPE
{
  GT_ORTHOGONAL,       // Purely orthogonal grid
  GT_ORTHOGONAL_MOVED, // An orthogonal grid with displaced vertices
  GT_REACTOR,          // Specialized geometry for reactors
  GT_ARBITRARY         // Arbitrary grid
};

typedef enum
{
  GEOMETRY_XY,
  GEOMETRY_XYZ,
  GEOMETRY_INVALID
} GEOMETRY_TYPE;

typedef enum
{
  DO_NOTHING,          // No flux fixup
  S_TO_ZERO,           // Set-to-zero flux fixup
} FIXUP_TYPE;

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

enum REORDERING
{
  none,
  linear,
  logarithmic
};

struct Kind_Info
{
  /*!
      \brief  Indicates whether using a CSZ (Consistent-Set-to-Zero) spatial method
  */
  bool is_csz;
  
  /*!
      \brief  Indicates whether or not the problem is time dependent.
  */
  bool is_time_dependent;

  /*!
     \brief Indicates whether or not the problem uses the k-eigenvalue solver.
  */
  bool has_k_eigenvalue;

  /*!
   \brief Indicates whether or not we are dealing with a depletion problem
  */
  bool depletion_problem;

  /*!
       \brief Indicates whether or not the problem has a QOI
  */
  bool has_QOI;

  /*!
       \brief Is a gradient calculation (adjoint problem)?
  */
  bool gradient_calculation;

  /*!
       \brief Is a steady-state adjoint problem?
  */
  bool adjoint_mode;

  /*!
       \brief Allocate forward vectors?
  */
  bool allocateForwardVectors;

  /*!
       \brief Allocate residual vectors?  (vectors that hold residual estimates at each
       depletion time step
  */
  bool allocateResidualVectors;

  /*!
      \brief This flag defaults to true to indicate the problem contains
      upscattering.  The user can set it to false if no upscattering is present
      in the problem or if the groupsets have been ordered properly.  If set to
      false, no across-groupset convergence checking is done.
  */
  bool _upscattering;

  /*!
      \brief This flag controls whether or not cells will keep a copy of their
      vertex coordinates.  This is very memory-expensive, but it is needed for
      a certain flux fixup.
  */
  bool _cells_keep_coordinates;

  /*!
      \brief This flag controls whether or not we compute pincell powers.
  */
  bool _pin_power;

  /*!
    \brief Problem type.  Valid types are NEUTRONICS, RADIATIVE_TRANSFER,
    NEUTRON_PHOTON, GAMMA, ELECTRON.
  */
  PROBLEM_KIND _kind;

  /*!
    \brief Spatial discretization method.
  */
  SPATIAL_METHOD _sdm;
  FEM_TYPE _fem;

  /*!
    \brief Type of spatial grid.
  */
  GRID_TYPE _grid_type;

  /*!
    \brief Number of dimensions in the problem.
  */
  int _dimensions;

  /*!
    \brief Number of angular moments in the problem.
  */
  int _moments;

  /*!
    \brief Number of energy groups in the problem.
  */
  int _energy_groups;

  /*!
    \brief Number of energy groups for uncollided flux calculation.
  */
  int _uf_energy_groups;

  /*!
    \brief Parameters for time-dependent problems.
  */
  TD_Info _td;

  /*!
    \brief Parameters for time-dependent problems.
  */
  KEFF_Info _k;

  /*!
   \brief Parameters for burnup problems.
  */
  BP_Info _bp;

  /*!
    \brief Parameters that describe a QOI.
  */

  QOI_Info _qoi;

  /*!
    \brief struct containing information about the depletion adjoint controller
  */
  adjointController _adjointController;

  /*!
    \brief Bool to determine whether to allocate for full solution
  */
  bool store_full_solution;

  /*!
    \brief A problem ID, used in checkpoint file naming
  */
  std::string problemID;

  /*!
    \brief Data type for electron problems.
  */
  int bfp;

  /*!
    \brief Group re-order method for electron problems.
  */
  REORDERING reordering;

  bool uncollided;
  bool first_collision_source;

  UNIT_TYPE phi_units;
  UNIT_TYPE ard_units;

  void print(std::ostream& os)
  {
    os << "Problem type is ";
    switch(_kind)
    {
      case NEUTRONICS:
        os << "neutronics";
        break;
      case RADIATIVE_TRANSFER:
        os << "radiative transfer";
        break;
      case NEUTRON_PHOTON:
        os << "coupled neutron/photon";
        break;
      case GAMMA:
        os << "gamma";
        break;
      case ELECTRON:
        os << "electron";
        break;
      case K_EIGEN:
        os << "k eigenvalue";
        break;
    }
    os << std::endl;
  };

  void define_type(stapl::typer& t)
  {
    t.member(is_csz);
    t.member(is_time_dependent);
    t.member(has_k_eigenvalue);
    t.member(depletion_problem);
    t.member(has_QOI);
    t.member(gradient_calculation);
    t.member(adjoint_mode);
    t.member(allocateForwardVectors);
    t.member(allocateResidualVectors);
    t.member(_upscattering);
    t.member(_cells_keep_coordinates);
    t.member(_pin_power);
    t.member(_kind);
    t.member(_dimensions);
    t.member(_td);
    t.member(_k);
    t.member(_bp);
    t.member(_qoi);
    t.member(_adjointController);
    t.member(store_full_solution);
    t.member(problemID);
    t.member(bfp);
    t.member(reordering);
    t.member(uncollided);
    t.member(first_collision_source);
    t.member(phi_units);
    t.member(ard_units);
  }
};
#endif
